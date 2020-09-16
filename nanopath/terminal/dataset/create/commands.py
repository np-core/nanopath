import click

from pathlib import Path
from nanopath.dataset import Dataset


@click.command()
@click.option(
    '--name',
    '-n',
    type=str,
    required=True,
    help='Name of the dataset (short)'
)
@click.option(
    '--collection',
    '-c',
    required=True,
    type=str,
    help='Dataset collection'
)
@click.option(
    '--version',
    '-v',
    required=True,
    type=int,
    help='Dataset version'
)
@click.option(
    '--author',
    '-a',
    required=True,
    type=str,
    help='Dataset author'
)
@click.option(
    '--data',
    '-d',
    type=Path,
    required=True,
    help='Path to directory containing data for dataset'
)
def create(name, data, collection, author, version):

    """ Create a new dataset from a data directory """

    ds = Dataset(
        name=name,
        collection=collection,
        version=version,
        data=data,
        author=author
    )

    ds.create_template_directory()
    ds.create_metadata()