from pathlib import Path
import click

from nanopath.surveillance import BioSampler


@click.command()
@click.option(
    '--file', '-f', type=Path, required=False,
    help='Data table to extract accessions from'
)
@click.option(
    '--outfile', '-o', type=Path, required=False, default="[queries.json]",
    help='Output response from queries by query key JSON format [queries.json]'
)
@click.option(
    '--query_column',  type=str, required=False, default="sample",
    help='Column in data containing accessions to query'
)
def query(file, outfile, query_column):

    """ Query and process metadata from the BioSample DB """

    bs = BioSampler()

    bs.process_list(
        file=file,
        sep='\t',
        query_column=query_column,
        outfile=outfile
    )

