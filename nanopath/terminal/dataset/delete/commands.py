import click
import os

from pathlib import Path
from nanopath.dataset import DataStorage


@click.command()
@click.option(
    '--name',
    '-n',
    type=Path,
    required=True,
    help='Exact name of dataset directory in storage on node'
)
@click.option(
    '--config',
    '-c',
    default=None,
    type=Path,
    help='Netflow node configuration file with entry: storage'
)
def delete(name, config):

    """ Delete a dataset from storage """

    if config is None:
        try:
            config = Path(os.environ['NP_STORAGE_CONFIG'])
        except KeyError:
            print('Could not determine configuration file path from environment variable: NP_STORAGE_CONFIG')
            exit(1)

    ds = DataStorage(data=None, config=config)

    for node, node_client in ds.netflow.clients.items():
        dataset = f"{node_client.storage_path}/{name}"
        ds.logger.info(f"Checking for dataset to delete: {dataset}")
        # Check if dataset exists on node
        out = node_client.execute_cmd(f'[ -d "{dataset}" ] && echo "1"')
        if not out:
            ds.logger.info(
                f'Could not detect storage directory: {dataset}'
            )
            exit(1)

        click.confirm(f"Do you want to delete the dataset in storage at: {dataset} ({node_client.host})")
        out = node_client.execute_cmd(f'rm -r {dataset} && echo "1"')
        if not out:
            ds.logger.info(
                f'Could not delete storage directory: {dataset}'
            )
            exit(1)

        node_client.disconnect()