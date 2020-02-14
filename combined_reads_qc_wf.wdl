import "odis_workflows/tools/fastqc.wdl" as FASTQC
import "odis_workflows/tools/trimmomatic.wdl" as TRIM

workflow combined_reads_qc_wf{
    Array[File] fastqs
    File adapters

    Int seed_mismatches = 2
    Int pal_clip_thresh = 15
    Int simple_clip_thresh = 10
    Int min_adapter_len = 7
    Int sliding_window_qual = 20
    Int sliding_window_len = 5
    Int min_len = 50

    scatter(fastq in fastqs){
        call FASTQC.FastQC as raw_fastqc{
            input:
                fastq = fastq
        }
        call TRIM.trimmomatic_se as trim{
            input:
                fastq = fastq,
                output_basename = basename(fastq, ".fastq.gz"),
                adapters = adapters,
                seed_mismatches = seed_mismatches,
                pal_clip_thresh = pal_clip_thresh,
                simple_clip_thresh = simple_clip_thresh,
                min_adapter_len = min_adapter_len,
                sliding_window_qual = sliding_window_qual,
                sliding_window_len = sliding_window_len,
                min_len = min_len
        }
        call FASTQC.FastQC as trim_fastqc{
            input:
                fastq = trim.trimmed_fastq
        }
    }

    output{
        Array[File] raw_html_reports = raw_fastqc.html_report
        Array[File] raw_test_reports = raw_fastqc.text_report
        Array[File] trimmed_fastqs = trim.trimmed_fastq
        Array[File] trim_logs = trim.trim_log
        Array[File] trim_html_reports = trim_fastqc.html_report
        Array[File] trim_test_reports = trim_fastqc.text_report
    }
}