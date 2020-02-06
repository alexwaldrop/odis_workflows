task FastQC {
   File fastq
   Boolean? nogroup = true
   Boolean? casava
   Boolean? nano
   Boolean? nofilter
   Int? min_length
   String? format
   File? contaminants
   File? adapters
   File? limits
   Int? kmers
   Boolean? quiet

   # Get the basename of the FastQC output files
   String output_basename = sub(basename(sub(fastq, "\\.gz$","")), "\\.[^\\.]*$", "_fastqc")

   # Runtime environment
   String docker = "rticode/fastqc:0.11.8"
   Int cpu = 2
   Int mem_gb = 8
   Int max_retries = 3

   meta {
    description: "FastQC task will compute basic QC stats for single Fastq file"
   }

   parameter_meta {
    fastq: "input fastq file"
    casava: "Files come from raw casava output"
    nano: "Files come from nanopore sequences and are in fast5 format"
    nogroup: "(optional) Disable grouping of bases for reads >50bp"
    nofilter: "(optional)  If running with --casava then don't remove read flagged by casava as poor quality when performing the QC analysis"
    min_length: "(optional) specify minimum read length of summarized reads"
    format: "(optional) override auto-detected format"
    contaminants: "(optional) file of contaminant reads to check for"
    adapters: "(optional) file of adapter sequences to check for"
    limits: "(optional) specify non-default file for warn/fail/pass limits for specific tests"
    kmers: "(optional) Specifies the length of Kmer to look for in the Kmer content module. Defaults to 7"
    docker: "(optional) the docker image containing the runtime environment for this task"
    mem_gb: "(optional) the amount of memory (GB) to provision for this task"
    cpu: "(optional) the number of cpus to provision for this task"
    }

   command <<<
      #FastQC takes a FastQ file and runs a series of tests on it to generate a comprehensive QC report
      fastqc --outdir ./ --extract \
            ${true='--casava' false='' casava} \
            ${true='--nano' false='' nano} \
            ${true='--nofilter' false='' nofilter} \
            ${true='--nogroup' false='' nogroup} \
            ${'--min_length ' + min_length} \
            ${'--format ' + format} \
            ${'--contaminants ' + contaminants} \
            ${'--adapters ' + adapters} \
            ${'--limits ' + limits} \
            ${'--kmers ' + kmers} \
            ${true='--quiet' false='' quiet} \
            ${fastq}
   >>>

   runtime {
        docker: docker
        cpu: cpu
        memory: "${mem_gb} GB"
        maxRetries: max_retries
   }

   output{
        File html_report = "${output_basename}.html"
        File text_report = "${output_basename}/fastqc_data.txt"
        File zip_report = "${output_basename}.zip"
   }
}
