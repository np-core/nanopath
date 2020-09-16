import click
import os

from pathlib import Path
from nanopath.dataset import DataStorage


@click.command()
@click.option(
    '--dataset',
    '-d',
    type=str,
    required=True,
    help='Exact name of dataset directory in storage on node'
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
def delete(name, version, collection, dataset, config):

    """ Delete a dataset from storage """

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
        dat = f"{node_client.storage_path}/{dataset}"
        ds.logger.info(f"Checking for dataset to delete: {dat}")
        # Check if dataset exists on node
        out = node_client.execute_cmd(f'[ -d "{dat}" ] && echo "1"')
        if not out:
            ds.logger.info(
                f'Could not detect storage directory: {dat}'
            )
            exit(1)

        click.confirm(f"Do you want to delete the dataset in storage at: {dat} ({node_client.host})")
        out = node_client.execute_cmd(f'rm -r {dat} && echo "1"')
        if not out:
            ds.logger.info(
                f'Could not delete storage directory: {dat}'
            )
            exit(1)

        node_client.disconnect()