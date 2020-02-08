task create_sequence_dictionary{
    File fasta
    String output_basename = basename(fasta, ".fasta")

    # Runtime environment
    String docker = "broadinstitute/picard:latest"
    Int cpu = 4
    Int mem_gb = 16
    Int max_retries = 3

    command <<<
        java -jar /usr/picard/picard.jar CreateSequenceDictionary R=${fasta} O=${output_basename}.dict
    >>>

    runtime {
        docker: docker
        cpu: cpu
        memory: "${mem_gb} GB"
        maxRetries: max_retries
    }

    output{
        File fasta_dict = "${output_basename}.dict"
    }
}