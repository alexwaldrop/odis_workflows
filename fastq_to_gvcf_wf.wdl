import "odis_workflows/tools/multiqc.wdl" as MULTIQC
import "odis_workflows/tools/trimmomatic.wdl" as TRIM
import "odis_workflows/tools/bwa.wdl" as BWA
import "odis_workflows/tools/samtools.wdl" as SAM
import "odis_workflows/tools/picard.wdl" as PICARD
import "odis_workflows/tools/utils.wdl" as UTILS
import "odis_workflows/tools/gatk.wdl" as GATK
import "odis_workflows/helper_workflows/collect_large_file_list_wf.wdl" as ZIP


workflow fastq_to_gvcf_wf{

    String sample_name
    String batch_id
    File fastq
    Int headcrop_len
    Boolean is_phred_33 = true

    # Make read group for bam
    String rg_pu
    String rg_pl
    String rg = "@RG\\tID:${sample_name}\\tPU:${rg_pu}\\tSM:${sample_name}\\tLB:${batch_id}\\tPL:${rg_pl}"

    File ref_fasta
    File ref_fasta_idx
    File ref_fasta_dict
    File ref_bwa_amb
    File ref_bwa_ann
    File ref_bwa_bwt
    File ref_bwa_pac
    File ref_bwa_sa

    # Headrop
    call TRIM.headcrop{
        input:
            fastq = fastq,
            crop_len = headcrop_len,
            output_basename = sample_name
    }

    # BWA align
    call BWA.bwa_mem_se{
        input:
            fastq = headcrop.trimmed_fastq,
            rg = rg,
            ref_fasta = ref_fasta,
            ref_fasta_idx = ref_fasta_idx,
            amb = ref_bwa_amb,
            ann = ref_bwa_ann,
            bwt = ref_bwa_bwt,
            pac = ref_bwa_pac,
            sa  = ref_bwa_sa,
            output_basename = sample_name
    }


    # Flagstat bam
    call SAM.Samtools_flagstat as flagstat{
        input:
            bam = bwa_mem_se.bam,
            bam_index = bwa_mem_se.bam_index
    }

    # Get names of chromosomes for passing to haplotype caller
    call UTILS.get_chrs_from_faidx as get_ref_chrs{
        input:
            ref_fasta_idx = ref_fasta_idx
    }

    # Scatter Haplotype caller across chromosomes
    call GATK.HaplotypeCallerGvcf_GATK4 as hap_caller{
        input:
            input_bam = bwa_mem_se.bam,
            input_bam_index = bwa_mem_se.bam_index,
            ref_dict = ref_fasta_dict,
            ref_fasta = ref_fasta,
            ref_fasta_index = ref_fasta_idx,
            chrs = get_ref_chrs.chrs,
            output_basename = sample_name
    }

    call MULTIQC.Collect_qc_files as gather_qc{
        input:
            input_files = [flagstat.flagstat_log],
            output_dir_name  = "${batch_id}_fastqc"
    }

    call MULTIQC.MultiQC{
        input:
            analysis_dir = gather_qc.output_dir,
            force = true,
            dirs_depth = 1,
            cpu = 2,
            mem_gb = 4
    }

    output{
        File final_bam = bwa_mem_se.bam
        File final_bam_index = bwa_mem_se.bam_index
        File gvcf = hap_caller.output_vcf
        File gvcf_index = hap_caller.output_vcf_index
        File multiqc_input_dir = gather_qc.output_dir
        File multiqc_report = MultiQC.multiqc_report
        Array[File] multiqc_data_files = MultiQC.multiqc_data_files
    }
}