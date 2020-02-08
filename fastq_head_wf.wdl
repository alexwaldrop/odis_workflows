task head{
    File in_file
    Int n
    String output_filename

    command <<<
        if [[ ${in_file} =~ \.gz$ ]]
        then
            zcat ${in_file} | head -n ${n} > ${output_filename}
        else
           head -n ${n} ${in_file} > ${output_filename}
        fi
    >>>

    runtime {
        docker: "ubuntu:18.04"
        cpu: 1
        memory: "1 GB"
    }

    output {
        File out = output_filename
    }
}

workflow fastq_head_wf{
    Array[File] fastqs
    Int n

    scatter(fastq in fastqs){
        call head{
            input:
                in_file = fastq,
                n = n,
                output_filename = basename(fastq, ".fastq.gz") + ".head.fastq"
        }
    }

    output{
        Array[File] head_fiiles = head.out
    }
}