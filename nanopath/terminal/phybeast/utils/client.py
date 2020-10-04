import click

from .show_dates import show_dates
from .subset_aln import subset_aln
from .extract_rate import extract_rate
from .remove_reference import remove_reference
from .date_random_test import date_random_test
from .prepare_metadata import prepare_metadata
from .remove_invariant import remove_invariant
from .rename_header import rename_header
from .print_header import print_header

VERSION = '0.1'


@click.group()
@click.version_option(version=VERSION)
def utils():

    """ Phybeast: utility tasks for pipeline execution """

    pass


utils.add_command(show_dates)
utils.add_command(subset_aln)
utils.add_command(extract_rate)
utils.add_command(remove_reference)
utils.add_command(date_random_test)
utils.add_command(prepare_metadata)
utils.add_command(remove_invariant)
utils.add_command(rename_header)
utils.add_command(print_header)