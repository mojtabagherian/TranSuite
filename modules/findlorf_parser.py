"""
Created on Wed May 15 14:25:00 2019
@author: Juan C Entizne
@email: e.entizne[at]dundee.ac.uk
"""
import os
import sys
import time
import traceback
from lib.logger.logger import logger, clean_log
from argparse import ArgumentParser, RawTextHelpFormatter
from modules.findlorf_main import findlorf_main
from lib.parsing.gtf_parsing_tools import add_features_to_gtf
from lib.tools.input_tools import check_input


description = \
    "Description:\n" + \
    "FindLORF find and annotates the longest ORF of a transcript as its putative CDS.\n"

parser = ArgumentParser(description=description, formatter_class=RawTextHelpFormatter,
                        add_help=False)

parser.add_argument('--gtf',
                    dest="gtf",
                    help="Transcriptome annotation file in GTF format.")

parser.add_argument('--fasta',
                    dest="fasta",
                    help='Transcripts fasta file (nucleotide sequence of exonic regions).')

parser.add_argument("--cds",
                    dest="cds_th", type=int, default=30,
                    help="Minimum number of amino-acids an ORF must have to be considered as a potential CDS. "
                         "Default: 30 AA.")

parser.add_argument('--outpath',
                    dest="outpath", default=None,
                    help="Path of the output folder.")

parser.add_argument('--outname',
                    dest="outname", default=None,
                    help="Prefix for the output files.")

parser.add_argument('--chimeric',
                    dest="chimeric", default=None,
                    help="Table indicating chimeric genes in the annotation.")


def main():
    args = parser.parse_args()

    # Check and sanitize user input
    args = check_input(args)

    # Create logfile to track the analysis, overwrite it if it exist ("w+" mode)
    time_stamp = time.strftime("%Y%m%d-%H%M%S")
    logfile = os.path.join(args.outpath, f"{time_stamp}_{args.outname}_logfile_temp.out")
    logger(logfile, w_mode="w+")

    # Record executed command
    command = " ".join(sys.argv)
    print(f"\n{command}\n", flush=True)

    # Run analysis
    try:
        transfind_gtf, orf_indes_file = findlorf_main(args.gtf, args.fasta, args.outpath, args.outname, args.cds_th)

        # Annotate additional features (UTR, start/stop codons)
        _ = add_features_to_gtf(transfind_gtf)

    except SystemExit as err:
        # Valid for python 3.5+
        print("".join(traceback.TracebackException.from_exception(err).format()))
    except Exception as err:
        print("".join(traceback.TracebackException.from_exception(err).format()))
        sys.exit(f"{err}")

    clean_log(logfile)
