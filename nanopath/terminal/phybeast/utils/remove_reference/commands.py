import click

from pathlib import Path
from nanopath.pathfinder import remove_sample


@click.command()
@click.option(
    "--alignment", "-a", default="input_alignment.fasta", help="Input alignment.", type=Path,
)
@click.option(
    "--output", "-o", default="output_alignment.fasta", help="Output alignment.", type=Path,
)
def remove_reference(alignment, output):

    """ Remove 'Reference' from Snippy alignment output file """

    remove_sample(alignment=alignment, outfile=output, remove='Reference')


