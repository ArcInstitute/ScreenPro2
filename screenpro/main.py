import argparse
import sys

# Custom functions
from .__init__ import __version__


def add_fq2count_parser(parent_subparsers, parent):
    name = "fq2count",
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

    # fq2count subcommand
    fq2count_parser = add_fq2count_parser(parent_subparsers, parent)
    