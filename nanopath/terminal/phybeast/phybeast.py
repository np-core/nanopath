import click

from .utils import utils
from .call_clair import call_clair
from .call_medaka import call_medaka

VERSION = "1"


@click.group()
@click.version_option(version=VERSION)
def phybeast():

    """ Outbreak phylogenetics and -dynamics pipeline in Nextflow """

    pass

phybeast.add_command(call_medaka)
phybeast.add_command(call_clair)
phybeast.add_command(utils)