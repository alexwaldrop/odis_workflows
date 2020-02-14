#!/usr/bin/env python3

import sys
import os
import argparse
import logging
import math

import pyfastx
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def configure_logging(verbosity):
    # Setting the format of the logs
    FORMAT = "[%(asctime)s] %(levelname)s: %(message)s"

    # Configuring the logging system to the lowest level
    logging.basicConfig(level=logging.DEBUG, format=FORMAT, stream=sys.stderr)

    # Defining the ANSI Escape characters
    BOLD = '\033[1m'
    DEBUG = '\033[92m'
    INFO = '\033[94m'
    WARNING = '\033[93m'
    ERROR = '\033[91m'
    END = '\033[0m'

    # Coloring the log levels
    if sys.stderr.isatty():
        logging.addLevelName(logging.ERROR, "%s%s%s%s%s" % (BOLD, ERROR, "ERROR", END, END))
        logging.addLevelName(logging.WARNING, "%s%s%s%s%s" % (BOLD, WARNING, "WARNING", END, END))
        logging.addLevelName(logging.INFO, "%s%s%s%s%s" % (BOLD, INFO, "INFO", END, END))
        logging.addLevelName(logging.DEBUG, "%s%s%s%s%s" % (BOLD, DEBUG, "DEBUG", END, END))
    else:
        logging.addLevelName(logging.ERROR, "ERROR")
        logging.addLevelName(logging.WARNING, "WARNING")
        logging.addLevelName(logging.INFO, "INFO")
        logging.addLevelName(logging.DEBUG, "DEBUG")

    # Setting the level of the logs
    level = [logging.ERROR, logging.WARNING, logging.INFO, logging.DEBUG][verbosity]
    logging.getLogger().setLevel(level)


def get_argparser():
    # Configure and return argparser object for reading command line arguments
    argparser_obj = argparse.ArgumentParser(prog="make_ref_contigs")

    def file_type(arg_string):
        if not os.path.exists(arg_string):
            err_msg = "%s does not exist! " \
                      "Please provide a valid file!" % arg_string
            raise argparse.ArgumentTypeError(err_msg)
        return arg_string

    # Path to input reference fasta
    argparser_obj.add_argument("--fasta",
                               "-f",
                               action="store",
                               type=file_type,
                               dest="in_fasta",
                               required=True,
                               help="Path to fasta formatted input reference contigs you want to merge")

    # Path to merged output fasta
    argparser_obj.add_argument("-output",
                               "-o",
                               action="store",
                               type=str,
                               dest="out_fasta",
                               required=True,
                               help="Merged fasta-formatted output file")


    # Path to input reference fasta
    argparser_obj.add_argument("--num-contigs",
                               "-c",
                               action="store",
                               type=int,
                               dest="num_contigs",
                               required=True,
                               help="Number of contigs desired in output fasta")

    # Length of N spacer region
    argparser_obj.add_argument("--n-spacer-len",
                               action="store",
                               type=int,
                               dest="spacer_len",
                               required=False,
                               default=150,
                               help="Length (bp) of N-spacer used to combine contigs")

    # Verbosity level
    argparser_obj.add_argument("-v",
                               action='count',
                               dest='verbosity_level',
                               required=False,
                               default=0,
                               help="Increase verbosity of the program."
                                    "Multiple -v's increase the verbosity level:\n"
                                    "0 = Errors\n"
                                    "1 = Errors + Warnings\n"
                                    "2 = Errors + Warnings + Info\n"
                                    "3 = Errors + Warnings + Info + Debug")

    return argparser_obj


def merge_contigs(fasta, num_contigs, min_contig_size, output_file, bed_output_file, spacer_len=150):
    seqs = []
    spacer = 'N' * spacer_len
    length = 0
    merged_contig_num = 1
    out_fh = open(output_file, "w")
    bed_out_fh = open(bed_output_file, "w")
    bed_out_fh.write("chrom\tchromStart\tchromEnd\torigContig\n")
    contig_pos = 0
    for seq in fasta:
        seq_name = seq.name
        seq = seq.seq
        length += len(seq)
        if length < min_contig_size or merged_contig_num >= num_contigs:
            # Keep adding to current merged contig until max length reached
            # OR this is the last merged contig (last merged contig can exceed max size to make sure you have the right number of chroms)
            seqs.append(str(seq))

        else:
            # Create new merged contig from sequences if current merged contig is long enough
            logging.info("Writing merged contig {0}...".format(merged_contig_num))
            # Use biopython for pretty printing
            merged_contig = SeqRecord(Seq(spacer.join(seqs)),
                                      id="merged_contig_{0}".format(merged_contig_num),
                                      description="")
            # Output to file
            out_fh.write("{0}".format(merged_contig.format("fasta")))

            # Reset counters
            seqs = []
            length = len(seq)
            seqs.append(str(seq))
            merged_contig_num += 1
            contig_pos = 0

        # Write bed record and update current contig position
        bed_out_fh.write("{0}\t{1}\t{2}\t{3}\n".format("merged_contig_{0}".format(merged_contig_num),
                                                       contig_pos,
                                                       contig_pos + len(seq),
                                                       seq_name,))
        contig_pos = contig_pos + len(seq) + spacer_len

    # Join and output remaining contigs
    if seqs:
        logging.info("Writing merged contig {0}...".format(merged_contig_num))
        merged_contig = SeqRecord(Seq(spacer.join(seqs)),
                                  id="merged_contig_{0}".format(merged_contig_num),
                                  description="")
        out_fh.write("{0}".format(merged_contig.format("fasta")))

    # Close output files
    out_fh.close()
    bed_out_fh.close()


def main():
    # Configure argparser
    argparser = get_argparser()

    # Parse the arguments
    args = argparser.parse_args()

    # Input files: json input file to be used as template and
    in_fasta = args.in_fasta
    num_contigs = args.num_contigs
    spacer_len = args.spacer_len
    out_fasta = args.out_fasta
    out_bed = "{0}.liftover.bed".format(out_fasta)

    # Configure logging appropriate for verbosity
    configure_logging(args.verbosity_level)

    # Count total length of bases in fasta
    logging.info("Reading input fasta...")
    fa = pyfastx.Fasta(in_fasta)
    total_contigs = len(fa)
    total_bases = fa.size
    logging.info("Total input contigs: {0}".format(total_contigs))
    logging.info("Total input bases: {0}".format(total_bases))

    # Determine minimum length of merged contigs
    min_contig_size = math.ceil(total_bases/float(num_contigs))

    # Merge contigs
    logging.info("Merging contigs")
    merge_contigs(fa, num_contigs, min_contig_size, out_fasta, out_bed, spacer_len)

    # Remove previous indices if they exist
    if os.path.exists("{0}.fxi".format(out_fasta)):
        os.remove("{0}.fxi".format(out_fasta))

    # Print info of new file
    logging.info("Summarizing output fasta...")
    fa_out = pyfastx.Fasta(out_fasta)
    out_total_contigs = len(fa_out)
    out_total_bases = fa_out.size
    logging.info("Total output contigs: {0}".format(out_total_contigs))
    logging.info("Total output bases: {0}".format(out_total_bases))

    # Check to make sure everything went well
    assert out_total_bases - (spacer_len*total_contigs) + (out_total_contigs*spacer_len) == total_bases

    # Log that output fasta has been fully validated
    logging.info("Output fasta validated! Output fasta size (bp) == "
                 "Input fasta size (bp) after accounting for N-spacers!")


if __name__ == "__main__":
    sys.exit(main())
