#!/usr/bin/env python3

import sys
import os
import argparse
import logging

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


def filter_fasta(fasta,
                 output_fasta,
                 gap_chars=True,
                 space_chars=True,
                 convert_uracil=True,
                 cap_seqs=True):
    """ filters fasta file(s), writes filtered output file(s)

    fasta_seqs: fasta filepath(s).  Multiple filepaths are separated by a
     colon.
    output_dir: output directory.
    gap_chars: if True, remove '.' and '-' chars.
    space_chars: if True, remove spaces from sequence.
    convert_uracil: if True, convert all "U" to "T".
    """

    # Gap characters to remove
    gap_chars_to_remove = []
    if gap_chars:
        gap_chars_to_remove = ['.', '-']
    if space_chars:
        gap_chars_to_remove.append(' ')

    # DNA characters to replace with the given value
    bad_dna_chars = {'U': 'T'}
    out_f = open(output_fasta, "w")

    count = 0
    for seq in fasta:
        name = seq.name
        seq = seq.seq

        if count % 1000 == 0:
            logging.info("Processed {0} records...".format(count))
        # Remove gaps and/or spaces
        for curr_gap_char in gap_chars_to_remove:
            seq = seq.replace(curr_gap_char, "")

        if cap_seqs:
            seq = seq.upper()

        # Replace bad dna characters
        if convert_uracil:
            for bad_dna_char, good_dna_char in bad_dna_chars.items():
                seq = seq.replace(bad_dna_char, good_dna_char)

        filtered = SeqRecord(Seq(seq),
                             id=name,
                             description="")
        # Output to file
        out_f.write("{0}".format(filtered.format("fasta")))
        count += 1

    out_f.close()


def main():
    # Configure argparser
    argparser = get_argparser()

    # Parse the arguments
    args = argparser.parse_args()

    # Input files: json input file to be used as template and
    in_fasta = args.in_fasta
    out_fasta = args.out_fasta

    # Configure logging appropriate for verbosity
    configure_logging(args.verbosity_level)

    # Count total length of bases in fasta
    logging.info("Reading input fasta...")
    fa = pyfastx.Fasta(in_fasta)
    total_contigs = len(fa)
    total_bases = fa.size
    logging.info("Total input contigs: {0}".format(total_contigs))
    logging.info("Total input bases: {0}".format(total_bases))

    # Merge contigs
    logging.info("Cleaning fastq records")
    filter_fasta(fa, out_fasta)


if __name__ == "__main__":
    sys.exit(main())
