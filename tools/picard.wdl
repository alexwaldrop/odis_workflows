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

task mark_duplicates{
    File input_bam
    File input_bam_index
    Boolean assume_sorted
    Boolean remove_duplicates
    String validation_stringency = "LENIENT"
    String output_basename

    # Runtime environment
    String docker = "broadinstitute/picard:latest"
    Int cpu = 4
    Int mem_gb = 8
    Int max_retries = 3

    command <<<
        # Mark duplicates
        java -jar /usr/picard/picard.jar MarkDuplicates \
            INPUT=${input_bam} \
            OUTPUT=${output_basename}.mrkdup.bam \
            METRICS_FILE=${output_basename}.mrkdup_report.txt \
            ASSUME_SORTED=${true="TRUE" false="FALSE" assume_sorted} \
            REMOVE_DUPLICATES=${true="TRUE" false="FALSE" remove_duplicates} \
            VALIDATION_STRINGENCY=${validation_stringency}

        # Re-build bam index
        java -jar /usr/picard/picard.jar BuildBamIndex INPUT=${output_basename}.mrkdup.bam
    >>>

    runtime {
        docker: docker
        cpu: cpu
        memory: "${mem_gb} GB"
        maxRetries: max_retries
    }

    output{
        File bam = "${output_basename}.mrkdup.bam"
        File bam_index = "${output_basename}.mrkdup.bam.bai"
        File mrkdup_report = "${output_basename}.mrkdup_report.txt"
    }
}