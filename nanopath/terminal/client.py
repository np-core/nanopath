import click

from .test import test
from .server import server
from .utils import utils

VERSION = '0.1'


@click.group()
@click.version_option(version=VERSION)
def terminal_client():
    """ NanoPath Server: Python CLI """
    pass


terminal_client.add_command(test)
terminal_client.add_command(utils)
terminal_client.add_command(server)
