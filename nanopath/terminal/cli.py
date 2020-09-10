import click

from .pathogen import pathogen
from .signal import signal
from .client import client
from .utils import utils
from .phybeast import phybeast
from .assembly import assembly
from .survey import survey

VERSION = '0.1'


@click.group()
@click.version_option(version=VERSION)
def terminal_client():
    """ NanoPath: Python CLI """
    pass


terminal_client.add_command(pathogen)
terminal_client.add_command(utils)
terminal_client.add_command(signal)
terminal_client.add_command(survey)
terminal_client.add_command(client)
terminal_client.add_command(phybeast)
terminal_client.add_command(assembly)