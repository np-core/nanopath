import click
from .server import server

VERSION = '0.1'


@click.group()
@click.version_option(version=VERSION)
def terminal_client():
    """ NanoPath Server: Python CLI """
    pass


terminal_client.add_command(server)