import click

from .merge_kraken import merge_kraken

VERSION = '0.1'


@click.group()
@click.version_option(version=VERSION)
def utils():

    """ Pathhogen: utility tasks for pipeline execution """

    pass


utils.add_command(merge_kraken)

