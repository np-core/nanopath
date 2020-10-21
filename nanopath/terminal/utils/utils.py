import click

from .rename_barcodes import rename_barcodes
from .snippy_core_to_clair import snippy_core_to_clair
from .get_nctc import get_nctc
from .assess import assess
from .compose import compose
from .combine_df import combine_df
from .sort_df import sort_df
from .merge_df import merge_df
from .plot_variant_summary import plot_variant_summary

VERSION = "1"


@click.group()
@click.version_option(version=VERSION)
def utils():
    """ Utility functions for NanoPath """
    pass


utils.add_command(rename_barcodes)
utils.add_command(assess)
utils.add_command(compose)
utils.add_command(combine_df)
utils.add_command(sort_df)
utils.add_command(merge_df)
utils.add_command(get_nctc)
utils.add_command(plot_variant_summary)
utils.add_command(snippy_core_to_clair)