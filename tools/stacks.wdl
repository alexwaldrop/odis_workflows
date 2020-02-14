task process_radtags_se{

    File input_fastq
    File barcode_file
    String barcode_option
    Boolean discard_low_quality
    Boolean rescue_reads
    Boolean clean_reads
    Boolean is_bestRAD
    String adapter_seq
    Int adapter_mismatches
    Boolean is_phred_33
    String enzyme_1

    Float? sliding_window_frac
    Int? sliding_window_score_limit
    Boolean? capture_discards
    Int? truncate_len
    String? output_type
    String? enzyme_2
    Boolean? retain_header
    Boolean? disable_rad_check
    Int? min_len
    Int? barcode_dist_1
    Int? barcode_dist_2
    Boolean? filter_illumina

    String fastq_basename = basename(input_fastq)

    # Runtime environment
    String docker = "alexwaldrop/stacks:2.5"
    Int cpu = 32
    Int mem_gb = 64
    Int max_retries = 3

    command <<<
        # Make output directory
        mkdir stacks_output

        # Run process radtags
        process_radtags -f ${input_fastq} \
            -b ${barcode_file} \
            -o stacks_output \
            --${barcode_option} \
            ${true='--clean' false='' clean_reads} \
            ${true='--rescue' false='' rescue_reads} \
            ${true='--quality' false='' discard_low_quality} \
            ${true='--bestrad' false='' is_bestRAD} \
            --adapter-1 "${adapter_seq}" \
            --adapter-mm ${adapter_mismatches} \
            --renz-1 ${enzyme_1} \
            ${true='-E phred33' false='-E phred64' is_phred_33} \
            ${'-w ' + sliding_window_frac} \
            ${'-s ' + sliding_window_score_limit} \
            ${'-y ' + output_type} \
            ${true='-D' false='' capture_discards} \
            ${true='--filter-illumina' false='' filter_illumina} \
            ${true='--disable-rad-check' false='' disable_rad_check} \
            ${'--len-limit ' + min_len} \
            ${'--barcode-dist-1 ' + barcode_dist_1} \
            ${'--barcode-dist-2 ' + barcode_dist_2} \
            ${true='--retain-header' false='' retain_header} \
            ${'--renz-2 ' + enzyme_2} \
            ${'-t ' + truncate_len}
    >>>

    runtime {
        docker: docker
        cpu: cpu
        memory: "${mem_gb} GB"
        maxRetries: max_retries
    }

    output{
        Array[File] demuxed_fastqs = glob("stacks_output/*.fq*")
        File discard_fastq = "stacks_output/${fastq_basename}.discards"
    }
}

task denovo_map{

    File input_fastq_zip
    File popmap_file
    Boolean? paired_end
    Boolean? remove_pcr_dups
    String? additional_options
    Int? mismatch_within_individuals
    Int? mismatch_between_individuals
    Float? var_alpha
    Float? gt_alpha
    Int? min_samples_per_pop
    Int? min_pops

    String input_fastq_dir = basename(input_fastq_zip, ".tar.gz")

    # Runtime environment
    String docker = "alexwaldrop/stacks:2.5"
    Int cpu = 48
    Int mem_gb = 72
    Int max_retries = 3

    command <<<
        # Make output directory
        mkdir stacks_output

        tar xzvf ${input_fastq_zip} -C .

        # Run denovo_map
        denovo_map.pl --samples ${input_fastq_dir} \
            --popmap ${popmap_file} \
            -o stacks_output \
            -T ${cpu} \
            ${true='--paired' false='' paired_end} \
            ${true='--rm-pcr-duplicates' false='' remove_pcr_dups} \
            ${'-X "' + additional_options + '"'} \
            ${'-M ' + mismatch_within_individuals} \
            ${'-n ' + mismatch_between_individuals} \
            ${'--var-alpha ' + var_alpha} \
            ${'--gt-alpha ' + gt_alpha} \
            ${'-r ' + min_samples_per_pop} \
            ${'-p ' + min_pops}
    >>>

    runtime {
        docker: docker
        cpu: cpu
        memory: "${mem_gb} GB"
        maxRetries: max_retries
    }

    output{
        Array[File] stacks_output = glob("stacks_output/*")
    }
}
