import "odis_workflows/tools/samtools.wdl" as SAM
import "odis_workflows/tools/picard.wdl" as PIC
import "odis_workflows/tools/bwa.wdl" as BWA
import "odis_workflows/tools/merge_ref_contigs.wdl" as MERGE

workflow build_references_wf{
    File ref_fasta
    Int max_ref_contigs
    String ref_output_basename
    Boolean do_contig_merge

    # Optionally collapse lots of smaller contigs into N larger contigs for easier processesing
    if(do_contig_merge){
        call MERGE.merge_ref_contigs as merge_ref{
            input:
                fasta = ref_fasta,
                max_contigs = max_ref_contigs,
                output_basename = ref_output_basename
        }
    }

    # Create fasta index of genome reference
    call SAM.Samtools_faidx as make_ref_index{
        input:
            fasta = select_first([merge_ref.merged_fasta, ref_fasta])
    }

    # Create sequence dictionary of genome reference for GATK tools
    call PIC.create_sequence_dictionary as make_ref_dict{
        input:
            fasta = select_first([merge_ref.merged_fasta, ref_fasta])
    }

    # Create BWA alignment index according to GATK best practices (bwtsw mode)
    call BWA.bwa_index{
        input:
            fasta = select_first([merge_ref.merged_fasta, ref_fasta]),
            use_bwtsw = true
    }

    output{
        File final_ref = select_first([merge_ref.merged_fasta, ref_fasta])
        Array[File] bwa_index_files = bwa_index.bwa_index_files
        File ref_index = make_ref_index.fasta_index
        File ref_dict  = make_ref_dict.fasta_dict
    }
}