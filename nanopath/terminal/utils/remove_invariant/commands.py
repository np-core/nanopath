import click

from pathlib import Path
from nanopath.processors import MegalodonCore


@click.command()
@click.option(
    "--alignment", "-a", default="input_alignment.fasta", help="Input alignment", type=Path,
)
@click.option(
    "--output", "-o", default="output_alignment.fasta", help="Output alignment", type=Path,
)
def remove_invariant(alignment, output):

    """ Remove invariant sites from an alignment """

    ml = MegalodonCore() # use method in MegalodonCore

    ml.filter_invariant_sites(alignment=alignment, fasta=output)


