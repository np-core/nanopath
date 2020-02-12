import click

from .process_kraken import process_kraken

VERSION = "0.1"


@click.group()
@click.version_option(version=VERSION)
def server():
    """ Functions for server-side compute and pipelines """
    pass


server.add_command(process_kraken)
