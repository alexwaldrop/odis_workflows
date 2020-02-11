import "odis_workflows/tools/utils.wdl" as UTILS
import "odis_workflows/helper_workflows/collect_large_file_list_wf.wdl" as ZIP

task ImportGVCFs {

    File gvcfs_zip
    File interval
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
    Int max_retries = 3

  command <<<
    set -euo pipefail

    rm -rf ${workspace_dir_name}

    # Unpack gvcfs
    tar -xzvf ${gvcfs_zip} -C ./

    # Create arg file
    ls ${gvcfs_dir}/*.g.vcf* | awk 'BEGIN{FS="\t"; OFS="\t"} {$1="-V ${gvcfs_dir}/"$1; print}' > argfile.txt

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
    File interval

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


workflow JointGenotyping {

    String analysis_name = "odis_ddrad"
    Array[String] sample_names
    Array[File] gvcfs

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

    # Collect all QC files into tarball for MultiQC input
    call ZIP.collect_large_file_list_wf as zip_gvcfs{
        input:
            output_dir_name = analysis_name + "_gvcfs",
            input_files = gvcfs
    }

    scatter (chr in get_ref_chrs.chrs) {
        call ImportGVCFs {
            input:
                gvcfs_zip = zip_gvcfs.output_dir,
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

    output{
        File output_vcf = GatherVcfs.output_vcf
        File output_vcf_index = GatherVcfs.output_vcf_index
    }
}


