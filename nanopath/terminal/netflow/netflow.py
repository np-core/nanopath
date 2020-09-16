import click

from .setup import setup
from .run import run

VERSION = "0.1"


@click.group()
@click.version_option(version=VERSION)
def netflow():
    """ Run the dashboard and pipeline netflow for real-time analysis """
    pass


netflow.add_command(setup)
netflow.add_command(run)
