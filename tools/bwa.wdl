task bwa_mem_se{
    File fastq
    String rg
    File bwa_index_tar_file
    String bwa_index_name = basename(bwa_index_tar_file, ".tar.gz")
    String output_basename
    String output_filename = "${output_basename}.bam"


    # Runtime environment
    String docker = "rticode/samtools:1.9"
    Int cpu = 16
    Int mem_gb = ceil(cpu*3)
    Int max_retries = 3

    command <<<
        tar -xvzf ${bwa_index_tar_file}

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

