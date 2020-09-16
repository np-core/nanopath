import click
import os

from pathlib import Path
from nanopath.dataset import DataStorage


@click.command()
@click.option(
    '--config',
    '-c',
    default=None,
    type=Path,
    help='DataStorage netflow configuration file with entry: storage'
)
def list(config):

    """ List the datasets present in storage """

    if config is None:
        try:
            config = Path(os.environ['NP_STORAGE_CONFIG'])
        except KeyError:
            print('Could not determine configuration file path from environment variable: NP_STORAGE_CONFIG')
            exit(1)

    ds = DataStorage(data=None, config=config)

    ds.list()
