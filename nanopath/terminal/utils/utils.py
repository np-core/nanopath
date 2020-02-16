import click

from .assess import assess
from .compose import compose
from .decontaminate import decontaminate

VERSION = "0.1"


@click.group()
@click.version_option(version=VERSION)
def utils():
    """ Utility functions for NanoPath """
    pass


utils.add_command(assess)
utils.add_command(compose)
utils.add_command(decontaminate)