from pathlib import Path
import click
import pandas

from nanopath.surveillance import Survey


@click.command()
@click.option(
    '--data', '-d', type=Path, default=None, required=False,
    help='Path to input data file tab-delineated with column headers [none]'
)
@click.option(
    '--accession', '-a', type=str, default=None,
    help='Column name in data file to download accessions [none]'
)
@click.option(
    '--project', '-p', type=str, default=None,
    help='Name of ENA project accession to compile [none]'
)
@click.option(
    '--species', '-s', type=str, default=None,
    help='Name of species to search and compile accessions for [none]'
)
@click.option(
    '--outdir', '-o', type=Path, default="query_out",
    help='Output directory for query [query_out]'
)
@click.option(
    '--tech', '-t', type=str, default=None,
    help='Combination of WGS sequencing techs, one of: illumina, ont, none [none]'
)
@click.option(
    '--sub', is_flag=True,
    help='Get the submitter read files if available for the search terms'
)
@click.option(
    '--term', type=str, default=None,
    help='Free-flow search tearm to ENA to return query form from [none]'
)
def find(data, accession, project, species, outdir, tech, sub, term):

    """ Find a collection short read genomes in ENA """

    outdir.mkdir(exist_ok=True, parents=True)

    su = Survey()

    sample_list = None
    if data is not None:
        su.read_query_file(file=data)
        if accession in su.query.columns:
            sample_list = su.query[accession].tolist()

    query_results = su.query_ena(
        study=project,
        sample=sample_list,
        species=species,
        scheme=tech,
        term=term,
        submitted_fastq=sub,
        allow_missing=False,
        allow_merged_fastq=False
    )

    for search, data in query_results.items():
        data['results'].to_csv(
            outdir / f'query_{search}.tsv', sep='\t', index=False
        )