task HaplotypeCaller {
    # Command parameters
    File input_bam
    File input_bam_index
    String interval
    String output_filename
    File ref_dict
    File ref_fasta
    File ref_fasta_index
    Boolean make_gvcf
    String gatk_path
    String? java_options

    # Runtime parameters
    String docker
    Int mem_gb = 8
    Int cpu = 1
    Int max_retries = 3

    String java_opt = select_first([java_options, "-XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10"])
    Int command_mem_gb = mem_gb - 1

  command {
    set -e

    ${gatk_path} --java-options "-Xmx${command_mem_gb}G ${java_opt}" \
      HaplotypeCaller \
      -R ${ref_fasta} \
      -I ${input_bam} \
      -L ${interval} \
      -O ${output_filename} \
      ${true="-ERC GVCF" false="" make_gvcf} \
      -G StandardAnnotation -G AS_StandardAnnotation -G StandardHCAnnotation
  }

  runtime {
    cpu: cpu
    docker: docker
    memory: mem_gb + " GB"
    maxRetries: max_retries
  }

  output {
    File output_gvcf = "${output_filename}"
    File output_gvcf_index = "${output_filename}.tbi"
  }
}


# Merge GVCFs generated per-interval for the same sample
task MergeGVCFs {
    # Command parameters
    Array[File] input_vcfs
    Array[File] input_vcfs_indexes
    String output_filename

    String gatk_path

    # Runtime parameters
    String docker
    Int cpu = 1
    Int mem_gb = 3
    Int command_mem_gb = mem_gb - 1
    Int max_retries = 3

  command {
  set -e

    ${gatk_path} --java-options "-Xmx${command_mem_gb}G"  \
      MergeVcfs \
      --INPUT ${sep=' --INPUT ' input_vcfs} \
      --OUTPUT ${output_filename}
  }

  runtime {
    docker: docker
    cpu: cpu
    memory: mem_gb + " GB"
    maxRetries: max_retries
  }

  output {
    File output_vcf = "${output_filename}"
    File output_vcf_index = "${output_filename}.tbi"
  }
}

workflow HaplotypeCallerGvcf_GATK4 {
    File input_bam
    File input_bam_index
    File ref_dict
    File ref_fasta
    File ref_fasta_index
    Array[String] chrs
    String output_basename

    Boolean make_gvcf = true
    String gatk_docker = "broadinstitute/gatk:4.1.4.0"
    String gatk_path = "/gatk/gatk"
    String gitc_docker = "broadinstitute/genomes-in-the-cloud:2.3.1-1500064817"
    String output_filename = "${output_basename}.g.vcf.gz"

  # Call variants in parallel over grouped calling intervals
  scatter (chr in chrs) {

    # Generate GVCF by interval
    call HaplotypeCaller {
      input:
        input_bam = input_bam,
        input_bam_index = input_bam_index,
        interval = chr,
        output_filename = "${output_basename}.${chr}.g.vcf.gz",
        ref_dict = ref_dict,
        ref_fasta = ref_fasta,
        ref_fasta_index = ref_fasta_index,
        make_gvcf = make_gvcf,
        docker = gatk_docker,
        gatk_path = gatk_path
    }
  }

  # Merge per-interval GVCFs
  call MergeGVCFs {
    input:
      input_vcfs = HaplotypeCaller.output_gvcf,
      input_vcfs_indexes = HaplotypeCaller.output_gvcf_index,
      output_filename = "${output_basename}.merged.g.vcf.gz",
      docker = gatk_docker,
      gatk_path = gatk_path
  }

  # Outputs that will be retained when execution is complete
  output {
    File output_vcf = MergeGVCFs.output_vcf
    File output_vcf_index = MergeGVCFs.output_vcf_index
  }
}