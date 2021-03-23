import click

from pathlib import Path
from nanopath.pathfinder import phybeast_prepare_metadata_file


@click.command()
@click.option(
    "--meta_data", "-m", default="meta.tsv", type=Path,
    help="Input meta data file, tab-delimited, includes: name, date columns.",
)
@click.option(
    "--alignment", "-a", default="input_alignment.fasta", type=Path,
    help="Input alignment to check dates against sequences",
)
@click.option(
    "--prep", "-p", default="lsd2", help="Prepare metadata file for one of: lsd2, treetime.", type=str,
)
@click.option(
    "--output", "-o", default="lsd2.meta", help="Output path of prepared meta data file.", type=Path,
)
def prepare_metadata(meta_data, alignment, prep, output):

    """ Randomise the dates in a meta data file with columns: name, date """

    phybeast_prepare_metadata_file(meta_file=meta_data, alignment=alignment, prep=prep, output_file=output)


