task Samtools_view{
    File input_bam
    File input_bam_index
    String output_basename
    String output_filename = "${output_basename}.bam"
    String? exclude_flag
    String? include_flag
    String? region

    # Runtime environment
    String docker = "rticode/samtools:1.9"
    Int cpu = 16
    Int mem_gb = 16
    Int max_retries = 3

    meta {
        description: "Samtools_view task produces a bam with optinal filtering by region and align flag values"
    }

    parameter_meta {
        bam: "Input bam file"
        bam_index: "Input bam index file"
        exclude_flag: "only include reads with all of the FLAGs in INT present"
        include_flag: "only include reads with none of the FLAGs in INT present"
        region: "only output reads from this region"
        docker: "(optional) the docker image containing the runtime environment for this task"
        mem_gb: "(optional) the amount of memory (GB) to provision for this task"
        cpu: "(optional) the number of cpus to provision for this task"
    }

    command <<<

        # Filter alignments
        samtools view \
            -@ ${cpu} \
            -h \
            -b \
            ${'-F ' + exclude_flag} \
            ${'-f ' + include_flag} \
            ${input_bam} \
            ${region} > ${output_filename}

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


task Samtools_sort{
    File input_bam
    String output_basename
    String output_filename = "${output_basename}.sorted.bam"
    Boolean sort_by_name = false

    # Runtime environment
    String docker = "rticode/samtools:1.9"
    Int cpu = 16
    Int mem_gb_per_cpu = 4
    Int mem_gb = ceil(mem_gb_per_cpu*cpu) + 8
    Int max_retries = 3

    meta {
        description: "Samtools_view task produces a bam with optinal filtering by region and align flag values"
    }

    parameter_meta {
        bam: "Input bam file"
        bam_index: "Input bam index file"
        docker: "(optional) the docker image containing the runtime environment for this task"
        mem_gb: "(optional) the amount of memory (GB) to provision for this task"
        cpu: "(optional) the number of cpus to provision for this task"
    }

    command <<<

        # sort alignments
        samtools sort \
            -@ ${cpu} \
            -m ${mem_gb_per_cpu}G \
            ${true="-n" false="" sort_by_name}\
            ${input_bam} > ${output_filename}

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

task Samtools_merge{
    Array[File] input_bams
    Array[File] input_bam_indices
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

        # merge alignments
        samtools merge \
            -f -c \
            -@ ${cpu} \
            ${output_filename} \
            ${sep=" " input_bams}

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


task Samtools_faidx{
    File fasta
    String output_basename = basename(fasta)

    # Runtime environment
    String docker = "rticode/samtools:1.9"
    Int cpu = 2
    Int mem_gb = 8
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
        samtools faidx ${fasta}
    >>>

    runtime {
        docker: docker
        cpu: cpu
        memory: "${mem_gb} GB"
        maxRetries: max_retries
    }

    output{
        File fasta_index = "${fasta}.fai"
    }
}

