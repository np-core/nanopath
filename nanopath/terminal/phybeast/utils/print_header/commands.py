import click

from pathlib import Path
from nanopath.pathfinder import print_fastx_header


@click.command()
@click.option(
    "--fasta", "-f", default="contigs.fasta", help="Input FASTA", type=Path,
)
def print_header(fasta):

    """ Rename FASTA headers """

    print_fastx_header(fastx=fasta)


