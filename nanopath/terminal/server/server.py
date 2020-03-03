import click

from .sepsis import sepsis

VERSION = "0.1"


@click.group()
@click.version_option(version=VERSION)
def server():
    """ Functions for server-side compute and pipelines """
    pass


server.add_command(sepsis)
