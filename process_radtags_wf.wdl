import "odis_workflows/tools/stacks.wdl" as STACKS
import "odis_workflows/tools/fastqc.wdl" as FASTQC
import "odis_workflows/tools/multiqc.wdl" as MULTIQC
import "odis_workflows/tools/utils.wdl" as UTILS
import "odis_workflows/helper_workflows/collect_large_file_list_wf.wdl" as ZIP

workflow process_radtags_wf{

    Array[String] batch_ids
    Array[File] fastqs
    Array[File] barcode_files
    String enzyme_1
    String adapter_seq
    Int adapter_mismatches
    Int min_len

    String barcode_option = "inline-null"
    Boolean discard_low_quality = true
    Boolean rescue_reads = true
    Boolean clean_reads = true
    Boolean is_phred_33 = true
    Float sliding_window_frac = 0.15
    Int sliding_window_score_limit = 10
    Boolean capture_discards = true
    Boolean retain_header = true
    Boolean disable_rad_check = false
    Boolean filter_illumina = true
    Boolean is_bestRAD = false

    # Process Radtags
    scatter(batch_index in range(length(batch_ids))){

        call STACKS.process_radtags_se as process_radtags{
            input:
                input_fastq = fastqs[batch_index],
                barcode_file = barcode_files[batch_index],
                enzyme_1 = enzyme_1,
                is_bestRAD = is_bestRAD,
                adapter_seq = adapter_seq,
                adapter_mismatches = adapter_mismatches,
                barcode_option = barcode_option,
                discard_low_quality = discard_low_quality,
                rescue_reads = rescue_reads,
                clean_reads = clean_reads,
                is_phred_33 = is_phred_33,
                sliding_window_frac = sliding_window_frac,
                sliding_window_score_limit = sliding_window_score_limit,
                capture_discards = capture_discards,
                retain_header = retain_header,
                disable_rad_check = disable_rad_check,
                min_len = min_len,
                filter_illumina = filter_illumina
        }
    }

    # Flatten 2-D array into single list of files
    call UTILS.flatten_string_array as flatten_fastqs{
        input:
            array = process_radtags.demuxed_fastqs
    }

    scatter(demux_fastq in flatten_fastqs.flat_array){
        # Had to include step because negative control files sometimes come back with no reads and fastqc fails
        if(size(demux_fastq, "B") > 100){
            call FASTQC.FastQC{
                input:
                    fastq = demux_fastq
            }
        }
    }

    call ZIP.collect_large_file_list_wf as gather_fastqc{
        input:
            input_files = select_all(FastQC.zip_report),
            output_dir_name  = "odis_ddrad_fastqc"
    }

    call MULTIQC.MultiQC{
        input:
            analysis_dir = gather_fastqc.output_dir,
            force = true,
            dirs_depth = 1,
            cpu = 8,
            mem_gb = 24
    }

    output{
        Array[File] demux_fastqs = flatten_fastqs.flat_array
        File multiqc_input_dir = gather_fastqc.output_dir
        File multiqc_report = MultiQC.multiqc_report
        Array[File] multiqc_data_files = MultiQC.multiqc_data_files
    }
}