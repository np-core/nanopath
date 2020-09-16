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
    help='Exact name of dataset in storage on node, format: {collection}-{name}-v{version}'
)
@click.option(
    '--outdir',
    '-o',
    type=Path,
    help='Output directory to retrieve data to, default dataset'
)
@click.option(
    '--config',
    '-c',
    default=None,
    type=Path,
    help='Netflow node configuration file with entry: storage'
)
def retrieve(dataset, outdir, config):

    """ Retrieve a dataset from storage """

    if config is None:
        try:
            config = Path(os.environ['NP_STORAGE_CONFIG'])
        except KeyError:
            print('Could not determine configuration file path from environment variable: NP_STORAGE_CONFIG')
            exit(1)

    if outdir is None:
        outdir = Path(dataset)

    ds = DataStorage(data=None, config=config)

    for node, node_client in ds.netflow.clients.items():
        ds.logger.info(f'Retrieve {dataset} from {node} to {outdir}')
        node_client.download_dataset_directory(dataset=dataset, outdir=outdir)
        node_client.disconnect()