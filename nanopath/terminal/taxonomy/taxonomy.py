import click

from .list import list

VERSION = "0.1"


@click.group()
@click.version_option(version=VERSION)
def dataset():

    """ Hybrid genome assembly and pangenome pipeline in Nextflow """

    pass


dataset.add_command(list)