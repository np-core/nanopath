from pathlib import Path
import click

from nanopath.surveillance import BioSampler


@click.command()
@click.option(
    '--file', '-f', type=Path, required=False,
    help='Data table to extract accessions from'
)
@click.option(
    '--output', '-o', type=Path, required=False, default="run_query.tsv",
    help='Output response from run queries [run_query.tsv]'
)
@click.option(
    '--query_column',  type=str, required=False, default="sample_accession",
    help='Column in data containing accessions to query'
)
def query(file, output, query_column):

    """ Query metadata by run accession from the BioSample DB """

    bs = BioSampler()

    bs.process_list(
        file=file,
        sep='\t',
        query_column=query_column,
        output=output
    )

