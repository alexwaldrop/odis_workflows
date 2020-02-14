task merge_ref_contigs{
    # Module for creating concatenated artificial contigs from ref file
    # Intended for poorly-scaffolded non-model organism genomes to make downstream processing easier

    File fasta
    Int max_contigs
    String output_basename
    Int? n_spacer_len

    # Runtime environment
    String docker = "alexwaldrop/merge_ref_contigs:07525f3"
    Int cpu = 4
    Int mem_gb = 16
    Int max_retries = 3

    command <<<
        merge_ref_contigs.py -f ${fasta} \
            -c ${max_contigs} \
            -o ${output_basename}.fasta \
            ${'-n-spacer-len ' + n_spacer_len} \
            -vv
    >>>

    runtime {
        docker: docker
        cpu: cpu
        memory: "${mem_gb} GB"
        maxRetries: max_retries
    }

    output{
        File merged_fasta = "${output_basename}.fasta"
        File liftover_bed = "${output_basename}.fasta.liftover.bed"
    }
}