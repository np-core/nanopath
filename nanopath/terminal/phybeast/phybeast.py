import click

from .utils import utils
from .call_clair import call_clair
from .call_medaka import call_medaka
from .call_megalodon import call_megalodon
from .snp_distance import snp_distance
from .forest_filter import forest_filter
from .forest_train import forest_train
from .forest_evaluate import forest_evaluate

VERSION = "1"


@click.group()
@click.version_option(version=VERSION)
def phybeast():

    """ Outbreak phylogenetics and -dynamics pipeline in Nextflow """

    pass


phybeast.add_command(forest_evaluate)
phybeast.add_command(forest_train)
phybeast.add_command(forest_filter)
phybeast.add_command(snp_distance)
phybeast.add_command(call_megalodon)
phybeast.add_command(call_medaka)
phybeast.add_command(call_clair)
phybeast.add_command(utils)