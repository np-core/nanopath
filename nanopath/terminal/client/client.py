import click

from .setup import setup
from .run import run

VERSION = "0.1"


@click.group()
@click.version_option(version=VERSION)
def client():
    """ Run the dashboard and pipeline client for real-time analysis """
    pass


client.add_command(setup)
client.add_command(run)
