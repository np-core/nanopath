import click

from .rename_barcodes import rename_barcodes

VERSION = "0.1"


@click.group()
@click.version_option(version=VERSION)
def signal():

    """ Signal-level analysis including Nextflow """

    pass


signal.add_command(rename_barcodes)