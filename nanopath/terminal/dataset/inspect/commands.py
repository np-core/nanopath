import click
import os

from pathlib import Path
from nanopath.dataset import DataStorage


@click.command()
@click.option(
    '--dataset',
    '-d',
    type=Path,
    required=True,
    help='Exact name of dataset in storage on node, format: {collection}-{name}-v{version}'
)
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
    '--config',
    '-c',
    default=None,
    type=Path,
    help='Netflow node configuration file with entry: storage'
)
def inspect(dataset, name, version, collection, config):

    """ Inspect a dataset in storage """

    if config is None:
        try:
            config = Path(os.environ['NP_STORAGE_CONFIG'])
        except KeyError:
            print('Could not determine configuration file path from environment variable: NP_STORAGE_CONFIG')
            exit(1)

    if name and version and collection:
        dataset = f"{collection}-{name}-v{version}"

    ds = DataStorage(data=None, config=config)

    for node, node_client in ds.netflow.clients.items():
        node_client.inspect_dataset(dataset=dataset)