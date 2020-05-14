import click

from pathlib import Path
from nanopath.pathfinder import phybeast_randomise_date_file


@click.command()
@click.option(
    "--date_file", "-d", default="input_alignment.fasta", help="Input alignment.", type=Path,
)
@click.option(
    "--output", "-o", default="output_alignment.fasta", help="Output alignment.", type=Path,
)
def randomise_dates(date_file, output):

    """ Randomise the dates in a meta data file with columns: name, date """

    phybeast_randomise_date_file(date_file=date_file, output_file=output)


