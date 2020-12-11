import click

from pathlib import Path
from nanopath.pathfinder import subset_alignment


@click.command()
@click.option(
    "--alignment", "-a", help="Input alignment", type=Path,
)
@click.option(
    "--data", "-d", help="Data file with columns name and and to subset", type=Path,
)
@click.option(
    "--column", "-c", help="Column with subset identifiers", type=str, default='outbreak'
)
@click.option(
    "--outdir", "-o", default="aln_subsets", help="Output path for subset alignments", type=Path,
)
@click.option(
    "--meta", "-m", default=None, help="Provide a meta data file to apply subset to as well", type=Path,
)
@click.option(
    "--prefix", "-p", default="", help="Prefix for outputs", type=str,
)
def subset_aln(alignment, outdir, data, column, meta, prefix):

    """ Subset an alignment by values in a data table column """

    subset_alignment(alignment=alignment, outdir=outdir, data=data, column=column, meta=meta, prefix=prefix)


