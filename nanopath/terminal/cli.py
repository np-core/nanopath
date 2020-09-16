import click

from .dataset import dataset
from .beastling import beastling
from .pathogen import pathogen
from .signal import signal
from .netflow import netflow
from .utils import utils
from .phybeast import phybeast
from .assembly import assembly
from .survey import survey


VERSION = '1'


@click.group()
@click.version_option(version=VERSION)
def terminal_client():
    """ NanoPath: Command Line Interface """
    pass


terminal_client.add_command(dataset)
terminal_client.add_command(beastling)
terminal_client.add_command(pathogen)
terminal_client.add_command(utils)
terminal_client.add_command(signal)
terminal_client.add_command(survey)
terminal_client.add_command(netflow)
terminal_client.add_command(phybeast)
terminal_client.add_command(assembly)