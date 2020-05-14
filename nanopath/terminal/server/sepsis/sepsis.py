import click

from .collect import collect
from .filter_host import filter_host

VERSION = "0.1"


@click.group()
@click.version_option(version=VERSION)
def sepsis():
    """ Sepsis command line tasks for pipeline support """
    pass


sepsis.add_command(collect)
sepsis.add_command(filter_host)