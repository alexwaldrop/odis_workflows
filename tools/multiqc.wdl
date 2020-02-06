task MultiQC{
    File analysis_dir # Expects a tarball because WDL doesn't do directories
    Boolean force = false
    Boolean dirs = false
    Int? dirs_depth
    Boolean full_names = false
    String? title
    String? comment
    String? file_name
    String? template
    String? tag
    String? ignore
    String? ignore_samples
    Boolean ignore_symlinks = false
    File? sample_names
    File? file_list
    Array[String]+? exclude
    Array[String]+? module
    Boolean data_dir = false
    Boolean no_data_dir = false
    String? data_format
    Boolean zip_data_dir = false
    Boolean export = false
    Boolean flat = false
    Boolean interactive = true
    Boolean lint = false
    Boolean pdf = false
    Boolean megaQC_upload = false
    File? config
    String? cl_config
    Boolean verbose  = false
    Boolean quiet = false
    String analysis_dir_name = basename(analysis_dir, ".tar.gz")

    # Name of report filename
    String report_filename = if (defined(file_name))
        then sub(select_first([file_name]), "\\.html$", "")
        else "multiqc"

    # Runtime options
    String docker = "rticode/multiqc:1.7"
    Int cpu = 2
    String mem_gb = 8
    Int max_retries = 3

    meta {
        description: "MultiQC task runs MultiQC on tar-zipped analysis directory"
    }

    parameter_meta {
        analysis_dir: "Tarballed (.tar.gz) directory containing analysis logs used as input for MultiQC"
        docker: "(optional) the docker image containing the runtime environment for this task"
        mem_gb: "(optional) the amount of memory (GB) to provision for this task"
        cpu: "(optional) the number of cpus to provision for this task"
    }

    command {
        set -e

        # Unzip analysis directory to working directory
        tar xzvf ${analysis_dir} -C .

        multiqc \
            ${true="--force" false="" force} \
            ${true="--dirs" false="" dirs} \
            ${"--dirs-depth " + dirs_depth} \
            ${true="--fullnames" false="" full_names} \
            ${"--title " + title} \
            ${"--comment " + comment} \
            ${"--filename " + file_name} \
            ${"--template " + template} \
            ${"--tag " + tag} \
            ${"--ignore " + ignore} \
            ${"--ignore-samples" + ignore_samples} \
            ${true="--ignore-symlinks" false="" ignore_symlinks} \
            ${"--sample-names " + sample_names} \
            ${"--file-list " + file_list} \
            ${true="--exclude " false="" defined(exclude)}${sep=" --exclude " exclude} \
            ${true="--module " false="" defined(module)}${sep=" --module " module} \
            ${true="--data-dir" false="" data_dir} \
            ${true="--no-data-dir" false="" no_data_dir} \
            ${"--data-format " + data_format} \
            ${true="--zip-data-dir" false="" zip_data_dir} \
            ${true="--export" false="" export} \
            ${true="--flat" false="" flat} \
            ${true="--interactive" false="" interactive} \
            ${true="--lint" false="" lint} \
            ${true="--pdf" false="" pdf} \
            ${false="--no-megaqc-upload" true="" megaQC_upload} \
            ${"--config " + config} \
            ${"--cl-config " + cl_config } \
            ${analysis_dir_name}
    }

    output {
        File multiqc_report = "${report_filename}_report.html"
        Array[File] multiqc_data_files = glob("${report_filename}_data/*")
    }

    runtime {
        docker: docker
        memory: "${mem_gb} GiB"
        cpu: cpu
        maxRetries: max_retries
    }
}

task Collect_qc_files{
    # Combine multiple files or directories into a single directory for MultiQC
    # Can be used to combine MultiQC directories across samples or aggregate files from a single sample
    # Right now this basically has to exist because WDL doesn't have a Directory type and MultiQC is based almost entirely on directory structure
    # NOTE: This task always assumes .tar.gz files should be unzipped into the parent directory
    # Allows for "recursively" combining directories to create more complex directory structures required for Multi-sample analysis
    # Was also intended for cases like Salmon output files that need to be in certain directories in order to be found (e.g. flenDist.txt)
    Array[File] input_files
    String output_dir_name

    # Runtime environment
    String docker = "ubuntu:18.04"
    Int cpu = 2
    Int mem_gb = 4
    Int max_retries = 3

    meta {
        description: "Gather multiple files into a single gzipped tarball that can be unzipped and directly input to MultiQC"
    }

    parameter_meta {
        input_files: "Files to zip"
        output_dir_name: "Name of directory that will be created and tarballed. Usually this should be a sample name"
        docker: "(optional) the docker image containing the runtime environment for this task"
        mem_gb: "(optional) the amount of memory (GB) to provision for this task"
        cpu: "(optional) the number of cpus to provision for this task"
    }

    command <<<
        set -e

        # Create directory and copy files into directory
        mkdir -p ${output_dir_name}

        # Loop through files in input file and copy/decompress them to output dir
        for input_file in ${sep=" " input_files}; do

            if [[ $input_file == *.tar.gz ]]; then
                # Untar directory into MultiQC dir
                # This is how we can combine compressed directories from multiple samples or modules into a single MultiQC directory
                tar -xvzf "$input_file" -C ${output_dir_name}
            else
                # Just copy flat files to MultiQC dir
                cp "$input_file" ${output_dir_name}
            fi

        done

        # Make a list of files in directory
        find ${output_dir_name}/* -type f > ${output_dir_name}.contents.txt

        # Compress directory so it can be placed inside higher-level MultiQC directories
        tar -cvzf ${output_dir_name}.tar.gz ${output_dir_name}

    >>>

    runtime {
        docker: docker
        cpu: cpu
        memory: "${mem_gb} GB"
        maxRetries: max_retries
    }

    output{
        File output_dir_file_list = "${output_dir_name}.contents.txt"
        File output_dir = "${output_dir_name}.tar.gz"
    }
}