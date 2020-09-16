import click
import os

from pathlib import Path
from nanopath.dataset import DataStorage


@click.command()
@click.option(
    '--data',
    '-d',
    type=Path,
    required=True,
    help='Path to directory containing data for dataset'
)
@click.option(
    '--config',
    '-c',
    default=None,
    type=Path,
    help='Netflow node configuration file with entry: storage'
)
def store(data, config):

    """ Create a new dataset from a data directory """

    if config is None:
        try:
            config = Path(os.environ['NP_STORAGE_CONFIG'])
        except KeyError:
            print('Could not determine configuration file path from environment variable: NP_STORAGE_CONFIG')
            exit(1)

    if not data.exists():
        print("File does not exist")
        exit(1)

    ds = DataStorage(data=data, config=config)

    for node, node_client in ds.netflow.clients.items():
        ds.logger.info(f'Uploading to node {node} to {node_client.storage_path}')
        node_client.upload_dataset_directory(local_path=data)
        node_client.disconnect()