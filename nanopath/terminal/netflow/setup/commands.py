import click
from pathlib import Path


@click.command()
@click.option(
    '--yaml',
    type=Path,
    help='Path to server and pipeline configuration file'
)
def setup(yaml):

    """ Setup server pipelines and resources """

    pass