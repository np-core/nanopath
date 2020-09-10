import click

from .collect import collect

VERSION = "0.1"


@click.group()
@click.version_option(version=VERSION)
def assembly():

    """ Hybrid genome assembly and pangenome pipeline in Nextflow """

    pass


assembly.add_command(collect)