import click

from .plot_assembly_dnadiff import plot_assembly_dnadiff
from .get_fast5 import get_fast5
from .plot_coverage_summary import plot_coverage_summary
from .rename_barcodes import rename_barcodes
from .snippy_core_to_clair import snippy_core_to_clair
from .get_nctc import get_nctc
from .assess import assess
from .compose import compose
from .combine_df import combine_df
from .sort_df import sort_df
from .merge_df import merge_df
from .plot_variant_violins import plot_variant_violins
from .plot_mlst_heatmap import plot_mlst_heatmap


VERSION = "1"


@click.group()
@click.version_option(version=VERSION)
def utils():
    """ Utility functions for NanoPath """
    pass


utils.add_command(plot_assembly_dnadiff)
utils.add_command(get_fast5)
utils.add_command(rename_barcodes)
utils.add_command(assess)
utils.add_command(compose)
utils.add_command(combine_df)
utils.add_command(sort_df)
utils.add_command(merge_df)
utils.add_command(get_nctc)
utils.add_command(plot_variant_violins)
utils.add_command(plot_mlst_heatmap)
utils.add_command(snippy_core_to_clair)
utils.add_command(plot_coverage_summary)
