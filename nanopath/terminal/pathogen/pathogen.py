import click

from .utils import utils
from .collect import collect
from .raven_mag import raven_mag
from .remove_human import remove_human

VERSION = "0.1"


@click.group()
@click.version_option(version=VERSION)
def pathogen():

    """ ONT pathogen detection in Nextflow """

    pass


pathogen.add_command(remove_human)
pathogen.add_command(raven_mag)
pathogen.add_command(utils)
pathogen.add_command(collect)
