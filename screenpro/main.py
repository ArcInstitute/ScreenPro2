import argparse
import sys
import os 
import pandas as pd
import polars as pl

from .__init__ import __version__
from . import ngs

def add_fq2cnt_parser(parent_subparsers, parent):
    name = "fq2cnt",
    desc = "Process FASTQ files to count sgRNA sequences."
    help = """
    Example usage:

        `screenpro fq2cnt --single-guide-design -l <path-to-library> -s <path-to-sample-sheet>`
        `screenpro fq2cnt --dual-guide-design -l <path-to-library> -s <path-to-sample-sheet>`
    
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
        "-s",
        "--sample_sheet",
        type=str,
        help="Path to sample sheet.",
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

    # fq2cnt subcommand
    fq2cnt_parser = add_fq2cnt_parser(parent_subparsers, parent)

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
        "fq2cnt": fq2cnt_parser,
    }

    if len(sys.argv) == 2:
        if sys.argv[1] in command_to_parser:
            command_to_parser[sys.argv[1]].print_help(sys.stderr)
        else:
            parent_parser.print_help(sys.stderr)
        sys.exit(1)
    