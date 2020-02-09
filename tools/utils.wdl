task append {
    Array[String] a
    String? b

    # Runtime environment
    String docker = "ubuntu:18.04"
    Int cpu = 1
    Int mem_gb = 1

    command {
        cat ${write_lines(a)}
        ${'echo ' + b}
    }

    runtime {
        docker: docker
        cpu: cpu
        memory: "${mem_gb} GB"
    }

    output {
        Array[File] out = read_lines(stdout())
    }
}

task rename_file {
    File in_file
    String new_name

    # Runtime environment
    String docker = "ubuntu:18.04"
    Int cpu = 1
    Int mem_gb = 1

    command {
        cp -R ${in_file} ${new_name}
    }

    runtime {
        docker: docker
        cpu: cpu
        memory: "${mem_gb} GB"
    }

    output {
        File out_file = "${new_name}"
    }
}

task slice{
    # Get 0-based slice of an array
    # end_index is exclusive so it works more or less like python
    Array[String] inputs
    Int start_pos
    Int end_pos
    Int slice_size = end_pos - start_pos

    # Make start pos 1-based because of how tail -N + works
    Int actual_start_pos = start_pos + 1

    # Runtime environment
    String docker = "ubuntu:18.04"
    Int cpu = 1
    Int mem_gb = 1

    command {
        tail -n +${actual_start_pos} ${write_lines(inputs)} | head -${slice_size}
    }

    runtime {
        docker: docker
        cpu: cpu
        memory: "${mem_gb} GB"
    }

    output {
        Array[String] outputs = read_lines(stdout())
    }
}