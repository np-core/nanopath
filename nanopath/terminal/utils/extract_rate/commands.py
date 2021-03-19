import click

from pathlib import Path
from nanopath.pathfinder import phybeast_extract_rate


@click.command()
@click.option(
    "--file", "-f", default="lsd.out", type=Path,
    help="Input meta data file, tab-delimited, includes: name, date columns.",
)
@click.option(
    "--prep", "-p", default="lsd2", help="Prepare output file from: lsd2, treetime.", type=str,
)
@click.option(
    "--output", "-o", default="rate.txt", type=Path,
    help="Rate and tMRCA estimates written to first line of output file for workflow collection.",
)
def extract_rate(file, prep, output):

    """ Extract substitution rate from LSD2 or TimeTree output """

    phybeast_extract_rate(result_file=file, prep=prep, output_file=output)


