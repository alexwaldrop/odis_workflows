import "odis_workflows/tools/fastqc.wdl" as FASTQC

workflow combined_reads_fastqc_wf{
    Array[File] fastqs

    scatter(fastq in fastqs){
        call FASTQC.FastQC{
            input:
                fastq = fastq
        }
    }

    output{
        Array[File] html_reports = FastQC.html_report
        Array[File] test_reports = FastQC.text_report
    }
}