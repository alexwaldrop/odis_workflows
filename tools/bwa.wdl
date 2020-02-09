task bwa_mem_se{
    File fastq
    String rg

    # BWA index files
    File amb
    File ann
    File bwt
    File pac
    File sa
    String bwa_index_name = basename(amb, ".amb")

    String output_basename
    String output_filename = "${output_basename}.bam"

    # Runtime environment
    String docker = "davelabhub/bwa_samtools:20190123"
    Int cpu = 12
    Int mem_gb = ceil(cpu*3)
    Int max_retries = 3

    command <<<
        # Align, convert to bam, and sort
        bwa mem -M -R "${rg}" \
            -t ${cpu} \
            ${bwa_index_name} \
            ${fastq} | samtools view -uS -@ ${cpu} - | samtools sort -@ ${cpu} - -o ${output_filename}

        # Index sorted output bam
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


task bwa_index{
    File fasta
    Boolean use_bwtsw
    String ref_basename = basename(fasta)

    # Runtime environment
    String docker = "davelabhub/bwa_samtools:20190123"
    Int cpu = 16
    Int mem_gb = 32
    Int max_retries = 3

    command <<<
        cp ${fasta} .
        bwa index -a ${true="bwtsw" false="" use_bwtsw} ${ref_basename}
    >>>

    runtime {
        docker: docker
        cpu: cpu
        memory: "${mem_gb} GB"
        maxRetries: max_retries
    }

    output{
        Array[File] bwa_index_files = glob("${ref_basename}.*")
    }
}

