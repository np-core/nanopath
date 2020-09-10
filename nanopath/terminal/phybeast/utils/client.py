import click

from .get_fast5 import get_fast5
from .show_dates import show_dates
from .subset_aln import subset_aln
from .sort_traits import sort_traits
from .join_traits import join_traits
from .extract_rate import extract_rate
from .remove_reference import remove_reference
from .date_random_test import date_random_test
from .prepare_metadata import prepare_metadata
from .remove_invariant import remove_invariant
from .plot_date_randomisation import plot_date_randomisation

VERSION = '0.1'


@click.group()
@click.version_option(version=VERSION)
def utils():

    """ Phybeast: utility tasks for pipeline execution """

    pass


utils.add_command(get_fast5)
utils.add_command(show_dates)
utils.add_command(subset_aln)
utils.add_command(sort_traits)
utils.add_command(join_traits)
utils.add_command(extract_rate)
utils.add_command(remove_reference)
utils.add_command(date_random_test)
utils.add_command(prepare_metadata)
utils.add_command(plot_date_randomisation)
utils.add_command(remove_invariant)
