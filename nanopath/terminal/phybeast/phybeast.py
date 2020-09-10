import click

from .utils import utils
from .beastling import beastling
from .call_medaka import call_medaka
from .call_megalodon import call_megalodon
from .snp_distance import snp_distance
from .snp_difference import snp_difference

VERSION = "0.1"


@click.group()
@click.version_option(version=VERSION)
def phybeast():

    """ Outbreak phylogenetics and -dynamics pipeline in Nextflow """

    pass


phybeast.add_command(snp_difference)
phybeast.add_command(snp_distance)
phybeast.add_command(call_megalodon)
phybeast.add_command(call_medaka)
phybeast.add_command(utils)
phybeast.add_command(beastling)