task headcrop{
    File fastq
    Int crop_len
    Boolean? is_phred33 = true

    # Output file
    String output_basename
    String trimmed_output_filename = "${output_basename}.headcropped.fastq.gz"

    # Runtime environment
    String docker = "rticode/trimmomatic:0.39"
    Int cpu = 4
    Int mem_gb = 8
    Int max_retries = 3

    command{
        java -jar /opt/trimmomatic \
            SE \
            ${true='-phred33' false='-phred64' is_phred33} \
            -threads ${cpu} \
            ${fastq} ${trimmed_output_filename} \
            HEADCROP:${crop_len} \
            > ${output_basename}.trimmomatic.log 2>&1
    }

    runtime{
        docker: docker
        cpu: cpu
        memory: "${mem_gb} GB"
        maxRetries: max_retries
    }

    output{
        File trimmed_fastq = "${trimmed_output_filename}"
        File trim_log = "${output_basename}.trimmomatic.log"
    }
}