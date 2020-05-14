import click

from .utils import client as utils_client
from .call_megalodon import call_megalodon

VERSION = "0.1"


@click.group()
@click.version_option(version=VERSION)
def phybeast():

    """ Outbreak phylogenetics and -dynamics pipeline in Nextflow """

    pass


phybeast.add_command(call_megalodon)
phybeast.add_command(utils_client)