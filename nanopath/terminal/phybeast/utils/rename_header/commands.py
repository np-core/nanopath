import click

from pathlib import Path
from nanopath.pathfinder import rename_fasta_header


@click.command()
@click.option(
    "--fasta", "-f", default="contigs.fasta", help="Input FASTA", type=Path,
)
@click.option(
    "--prefix", "-p", default="contig_", help="Prefix for headers", type=str,
)
@click.option(
    "--output", "-o", default="output_alignment.fasta", help="Output FASTA", type=Path,
)
def rename_header(fasta, prefix, output):

    """ Rename FASTA headers """

    rename_fasta_header(fasta=fasta, outfile=output, prefix=prefix)


