import "odis_workflows/tools/utils.wdl" as UTILS
import "odis_workflows/tools/samtools.wdl" as SAM
import "odis_workflows/helper_workflows/collect_large_file_list_wf.wdl" as ZIP

task ImportGVCFs {

    File gvcfs_zip
    String interval
    File ref_fasta
    File ref_fasta_index
    File ref_dict
    String  gvcfs_dir = basename(gvcfs_zip, ".tar.gz")
    String workspace_dir_name
    Int batch_size

    # Using a nightly version of GATK containing fixes for GenomicsDB
    # https://github.com/broadinstitute/gatk/pull/5899
    String gatk_docker = "us.gcr.io/broad-gotc-prod/gatk-nightly:2019-05-07-4.1.2.0-5-g53d015e4f-NIGHTLY-SNAPSHOT"
    Int cpu = 4
    Int mem_gb = 26
    Int max_retries = 1

  command <<<
    set -euo pipefail

    rm -rf ${workspace_dir_name}

    # Unpack gvcfs
    tar -xzvf ${gvcfs_zip} -C ./

    # Create arg file
    ls ${gvcfs_dir}/*.g.vcf.gz | awk 'BEGIN{FS="\t"; OFS="\t"} {$1="-V "$1; print}' > argfile.txt

    echo $(pwd)

    ls -l

    # We've seen some GenomicsDB performance regressions related to intervals, so we're going to pretend we only have a single interval
    # using the --merge-input-intervals arg
    # There's no data in between since we didn't run HaplotypeCaller over those loci so we're not wasting any compute

    # The memory setting here is very important and must be several GiB lower
    # than the total memory allocated to the VM because this tool uses
    # a significant amount of non-heap memory for native libraries.
    # Also, testing has shown that the multithreaded reader initialization
    # does not scale well beyond 5 threads, so don't increase beyond that.
    gatk --java-options -Xms8g \
      GenomicsDBImport \
      --genomicsdb-workspace-path ${workspace_dir_name} \
      --batch-size ${batch_size} \
      --arguments_file argfile.txt \
      -L ${interval} \
      --reader-threads 5 \
      --merge-input-intervals \
      --consolidate

    tar -cf ${workspace_dir_name}.tar ${workspace_dir_name}
  >>>

  runtime {
    memory: mem_gb + " GB"
    cpu: cpu
    docker: gatk_docker
    maxRetries: max_retries
  }

  output {
    File output_genomicsdb = "${workspace_dir_name}.tar"
  }
}


task GenotypeGVCFs {
    File workspace_tar
    String interval

    String output_vcf_filename

    File ref_fasta
    File ref_fasta_index
    File ref_dict

    # This is needed for gVCFs generated with GATK3 HaplotypeCaller
    Boolean allow_old_rms_mapping_quality_annotation_data = false
    String gatk_docker = "us.gcr.io/broad-gatk/gatk:4.1.4.0"
    Int cpu = 2
    Int mem_gb = 26
    Int max_retries = 3

  command <<<
    set -euo pipefail

    tar -xf ${workspace_tar}
    WORKSPACE=$(basename ${workspace_tar} .tar)

    gatk --java-options -Xms8g \
      GenotypeGVCFs \
      -R ${ref_fasta} \
      -O ${output_vcf_filename} \
      -G StandardAnnotation -G AS_StandardAnnotation \
      --only-output-calls-starting-in-intervals \
      --use-new-qual-calculator \
      -V gendb://$WORKSPACE \
      -L ${interval} \
      ${true='--allow-old-rms-mapping-quality-annotation-data' false='' allow_old_rms_mapping_quality_annotation_data} \
      --merge-input-intervals
  >>>

  runtime {
    memory: mem_gb + " GB"
    cpu: cpu
    docker: gatk_docker
    maxRetries: max_retries
  }

  output {
    File output_vcf = "${output_vcf_filename}"
    File output_vcf_index = "${output_vcf_filename}.tbi"
  }
 }

task GatherVcfs {

    Array[File] input_vcfs
    String output_vcf_name

    Int cpu = 1
    Int mem_gb = 7
    Int max_retries = 3
    String gatk_docker = "us.gcr.io/broad-gatk/gatk:4.1.1.0"

  command <<<
    set -euo pipefail

    # --ignore-safety-checks makes a big performance difference so we include it in our invocation.
    # This argument disables expensive checks that the file headers contain the same set of
    # genotyped samples and that files are in order by position of first record.
    gatk --java-options -Xms6g \
      GatherVcfsCloud \
      --ignore-safety-checks \
      --gather-type BLOCK \
      --input ${sep=" --input " input_vcfs} \
      --output ${output_vcf_name}

    tabix ${output_vcf_name}
  >>>

  runtime {
    memory: mem_gb + " GB"
    cpu: cpu
    docker: gatk_docker
    maxRetries: max_retries
  }

  output {
    File output_vcf = "${output_vcf_name}"
    File output_vcf_index = "${output_vcf_name}.tbi"
  }
}

task merge_lots_of_bams{
    File input_bams_dir
    String bam_dir_basename = basename(input_bams_dir, ".tar.gz")
    String output_basename
    String output_filename = "${output_basename}.merged.bam"

    # Runtime environment
    String docker = "rticode/samtools:1.9"
    Int cpu = 8
    Int mem_gb = ceil(cpu*4) + 4
    Int max_retries = 3

    meta {
        description: "Samtools_merge task merges 2 or more sorted bams into single sorted output bam"
    }

    parameter_meta {
        bam: "Input bam file"
        bam_index: "Input bam index file"
        docker: "(optional) the docker image containing the runtime environment for this task"
        mem_gb: "(optional) the amount of memory (GB) to provision for this task"
        cpu: "(optional) the number of cpus to provision for this task"
    }

    command <<<

        tar -xzvf ${input_bams_dir} -C ./

        # merge alignments
        samtools merge \
            -f -c \
            -@ ${cpu} \
            ${output_filename} \
            ${bam_dir_basename}/*.bam

        # Index output bam
        samtools index ${output_filename}

    >>>

    runtime {
        docker: docker
        cpu: cpu
        memory: "${mem_gb} GB"
        maxRetries: max_retries
    }

    output{
        File bam = output_filename
        File bam_index = "${output_filename}.bai"
    }
}

workflow JointGenotyping {

    String analysis_name = "odis_ddrad"
    Array[String] sample_names
    Array[File] gvcfs
    Array[File] gvcf_indexes
    Array[File] bams
    Array[File] bam_indexes

    File ref_fasta
    File ref_fasta_idx
    File ref_fasta_dict

    # Runtime attributes
    String gatk_docker = "broadinstitute/gatk:4.1.4.0"
    String gatk_path = "/gatk/gatk"
    String picard_docker = "us.gcr.io/broad-gotc-prod/gatk4-joint-genotyping:yf_fire_crosscheck_picard_with_nio_fast_fail_fast_sample_map"

    Boolean use_allele_specific_annotations = true
    Boolean cross_check_fingerprints = true
    Boolean scatter_cross_check_fingerprints = false

    # Get names of chromosomes for passing to haplotype caller
    call UTILS.get_chrs_from_faidx as get_ref_chrs{
        input:
            ref_fasta_idx = ref_fasta_idx
    }

    call ZIP.collect_large_file_list_wf as zip_gvcfs{
        input:
            output_dir_name = analysis_name + "_gvcfs",
            input_files = gvcfs
    }

    call ZIP.collect_large_file_list_wf as zip_gvcf_indexes{
        input:
            output_dir_name = analysis_name + "_gvcfs",
            input_files = gvcf_indexes
    }

    call ZIP.collect_chunks as zip_gvcfs_with_index{
        input:
            output_dir_name = analysis_name + "_gvcfs",
            input_files = [zip_gvcfs.output_dir, zip_gvcf_indexes.output_dir]
    }


    call ZIP.collect_large_file_list_wf as zip_bams{
        input:
            output_dir_name = analysis_name + "_bams",
            input_files = bams
    }

    call ZIP.collect_large_file_list_wf as zip_bam_indexes{
        input:
            output_dir_name = analysis_name + "_bams",
            input_files = bam_indexes
    }

    call ZIP.collect_chunks as zip_bams_with_index{
        input:
            output_dir_name = analysis_name + "_bams",
            input_files = [zip_bams.output_dir, zip_bam_indexes.output_dir]
    }

    scatter (chr in get_ref_chrs.chrs) {
        call ImportGVCFs {
            input:
                gvcfs_zip = zip_gvcfs_with_index.output_dir,
                interval = chr,
                ref_fasta = ref_fasta,
                ref_fasta_index = ref_fasta_idx,
                ref_dict = ref_fasta_dict,
                workspace_dir_name = "genomicsdb",
                batch_size = 50,
                gatk_docker = gatk_docker
        }

        call GenotypeGVCFs {
            input:
                workspace_tar = ImportGVCFs.output_genomicsdb,
                interval = chr,
                output_vcf_filename = analysis_name + "." + chr + ".vcf.gz",
                ref_fasta = ref_fasta,
                ref_fasta_index = ref_fasta_idx,
                ref_dict = ref_fasta_dict,
                gatk_docker = gatk_docker
        }
    }

    # Gather genotype VCFs from each chr into single VCF
    call GatherVcfs {
        input:
            input_vcfs = GenotypeGVCFs.output_vcf,
            output_vcf_name = analysis_name + ".merged.vcf.gz",
    }

    call merge_lots_of_bams as merge_bam{
        input:
            input_bams_dir = zip_bams_with_index.output_dir,
            output_basename = analysis_name + ".all_samples"
    }

    output{
        File output_vcf = GatherVcfs.output_vcf
        File output_vcf_index = GatherVcfs.output_vcf_index
        File merged_bam = merge_bam.bam
        File merged_bam_index = merge_bam.bam_index
    }
}


