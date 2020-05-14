import click

from .extract_rate import extract_rate
from .remove_reference import remove_reference
from .randomise_dates import randomise_dates
from .prepare_metadata import prepare_metadata
from .plot_date_randomisation import plot_date_randomisation

VERSION = '0.1'


@click.group()
@click.version_option(version=VERSION)
def utils():

    """ Phybeast: utility tasks for pipeline execution """

    pass


utils.add_command(extract_rate)
utils.add_command(remove_reference)
utils.add_command(randomise_dates)
utils.add_command(prepare_metadata)
utils.add_command(plot_date_randomisation)
