import click

from .run import run

VERSION = "0.1"


@click.group()
@click.version_option(version=VERSION)
def server():
    """ Live monitoring of sequencing runs """
    pass


server.add_command(sepsis)
