import "odis_workflows/helper_workflows/collect_large_file_list_wf.wdl" as ZIP
import "odis_workflows/tools/stacks.wdl" as STACKS

workflow stacks_denovo_wf{
    Array[File] input_fastqs
    File popmap_file

    call ZIP.collect_large_file_list_wf as zip_fastqs{
        input:
            output_dir_name = "stacks_input",
            input_files = input_fastqs
    }

    call STACKS.denovo_map{
        input:
            input_fastq_zip = zip_fastqs.output_dir,
            popmap_file = popmap_file
    }

    output{
        Array[File] stacks_output = denovo_map.stacks_output
    }

}