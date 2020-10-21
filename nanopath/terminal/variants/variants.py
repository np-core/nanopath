import click

from .hybrid_megalodon import hybrid_megalodon
from .hybrid_denovo import hybrid_denovo
from .snp_distance import snp_distance
from .forest_filter import forest_filter
from .forest_train import forest_train
from .forest_evaluate import forest_evaluate

VERSION = "1"


@click.group()
@click.version_option(version=VERSION)
def variants():

    """ Outbreak phylogenetics and -dynamics variant calling in Nextflow """

    pass


variants.add_command(hybrid_megalodon)
variants.add_command(hybrid_denovo)
variants.add_command(forest_evaluate)
variants.add_command(forest_train)
variants.add_command(forest_filter)
variants.add_command(snp_distance)