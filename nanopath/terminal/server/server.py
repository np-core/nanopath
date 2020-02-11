import click

from .process_kraken import process_kraken

VERSION = "0.1"


@click.group()
@click.version_option(version=VERSION)
def server():
    """ Tasks for: supplementary utility functions """
    pass


server.add_command(process_kraken)
