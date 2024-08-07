import argparse
import sys
import os 
import pandas as pd
import polars as pl
from glob import glob
from simple_colors import green
from .__init__ import __version__
from . import ngs


def add_counter_parser(parent_subparsers, parent):
    name = "guidecounter"
    desc = "Process FASTQ files to count sgRNA sequences."
    help = """
    Example usage:

        `screenpro guidecounter --single-guide-design -l <path-to-library> -s <path-to-sample-sheet>`
        `screenpro guidecounter --dual-guide-design -l <path-to-library> -s <path-to-sample-sheet>`
    
    """

    sub_parser = parent_subparsers.add_parser(
        name, parents=[parent], description=desc, help=help, add_help=True
    )

    sub_parser.add_argument(
        '-c',
        '--cas-type',
        type=str,
        default='cas9',
        help='Type of Cas protein used for the screen. Default: cas9'
    )

    sub_parser.add_argument(
        "--single-guide-design",
        action="store_true",
        help="Single guide design library.",
    )

    sub_parser.add_argument(
        "--dual-guide-design",
        action="store_true",
        help="Dual guide design library.",
    )

    sub_parser.add_argument(
        "-l",
        "--library",
        type=str,
        required=True,
        help="Path to library file.",
    )

    sub_parser.add_argument(
        "-p",
        "--path",
        type=str,
        required=True,
        help="Path to directory containing FASTQ files."
    )
    
    sub_parser.add_argument(
        "-s",
        "--samples",
        required=True,
        type=str,
        help="Sample ID(s) to process. Regex can be used to match multiple samples."
    )

    sub_parser.add_argument(
        "--write-count-matrix",
        action="store_true",
        help="Write count matrix to file."
    )

    sub_parser.add_argument(
        '-o',
        '--output',
        type=str,
        help='Path to output directory.'
    )

    return sub_parser


def main():
    """
    Function containing argparse parsers and arguments to allow the use of screenpro from the terminal.
    """

    # Define parent parser
    parent_parser = argparse.ArgumentParser(
        description=f"screenpro v{__version__}", add_help=False
    )
    # Initiate subparsers and define parent
    parent_subparsers = parent_parser.add_subparsers(dest="command")
    parent = argparse.ArgumentParser(add_help=False)

    # Add custom help argument to parent parser
    parent_parser.add_argument(
        "-h", "--help", action="store_true", help="Print manual."
    )
    # Add custom version argument to parent parser
    parent_parser.add_argument(
        "-v", "--version", action="store_true", help="Print version."
    )

    ## screenpro commands

    # counter subcommand
    counter_parser = add_counter_parser(parent_subparsers, parent)

    ## Define return values
    args = parent_parser.parse_args()

    # Help return
    if args.help:
        # Retrieve all subparsers from the parent parser
        subparsers_actions = [
            action
            for action in parent_parser._actions
            if isinstance(action, argparse._SubParsersAction)
        ]
        for subparsers_action in subparsers_actions:
            # Get all subparsers and print help
            for choice, subparser in subparsers_action.choices.items():
                print("Subparser '{}'".format(choice))
                print(subparser.format_help())
        sys.exit(1)

    # Version return
    if args.version:
        print(f"screenpro version: {__version__}")
        sys.exit(1)

    # Show help when no arguments are given
    if len(sys.argv) == 1:
        parent_parser.print_help(sys.stderr)
        sys.exit(1)

    # Show module specific help if only module but no further arguments are given
    command_to_parser = {
        "guidecounter": counter_parser,
    }

    if len(sys.argv) == 2:
        if sys.argv[1] in command_to_parser:
            command_to_parser[sys.argv[1]].print_help(sys.stderr)
        else:
            parent_parser.print_help(sys.stderr)
        sys.exit(1)
    
    ## counter return
    if args.command == "guidecounter":
        ### Perform checks on input arguments
        ## Check if library platform is provided by user
        if args.single_guide_design:
            args.library_type = "single_guide_design"
        elif args.dual_guide_design:
            args.library_type = "dual_guide_design"
        else:
            print("Library type not provided. Available options are '--single-guide-design' or '--dual-guide-design'. Exiting...")
            sys.exit(1)

        if args.output:
            pass
        else:
            print("No output directory provided. Exiting...")
            sys.exit(1)
        
        # Create saving directory
        directory = "/".join(args.output.split("/")[:-1])
        if directory != "":
            os.makedirs(directory, exist_ok=True)   
        
        ### Parse input files: library, samples, and fastq files
        counter = ngs.GuideCounter(args.cas_type, args.library_type)
        
        ## 1. Load library table and check if required columns are present
        counter.load_library(
            args.library, sep = '\t', index_col=False, 
            verbose=True
        )

        counter.library.to_pandas().to_csv(
            f"{args.output}/library.reformatted.tsv", sep = '\t',
        )
        print(f"Library table saved to {args.output}/library.reformatted.tsv")

        ## 2. Get list of samples to process
        #TODO: Implement regex to match multiple samples
        samples = args.samples.split(",")

        if len(samples) == 0:
            print(f"No samples found in {args.path}/{args.samples}. Exiting...")
            sys.exit(1)
        
        print(f"Samples to process: {','.join(samples)}")

        ## 3. Load FASTQ files and count sgRNA sequences
        counter.get_counts_matrix(
            fastq_dir = args.path,
            samples = samples,
            # get_recombinant=True,
            # write='force',
            verbose = True
        )
        print("Finished FASTQ processing.")

        # Save counts for each sample
        for sample in samples:
            cnt = counter.counts_mat[[sample]].rename(columns={sample: 'count'}).astype({'count': int})
            cnt.to_csv(f"{args.output}/{sample}.counts.tsv", sep='\t')
            print(f"Counts for sample {sample} saved to {args.output}/{sample}.counts.tsv")

        # Save count matrix
        if args.write_count_matrix:
            counter.counts_mat.to_csv(
                f"{args.output}/count_matrix.tsv",
                sep='\t'
            )
            print(f"Count matrix saved to {args.output}/count_matrix.tsv")

        # Pipeline finished
        print(green("Your run is finished successfully."))
        sys.exit(0)
