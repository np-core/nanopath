import click
import pandas
from nanopath.pipelines import PathogenPipeline
from pathlib import Path


@click.command()
@click.option(
    '--path',
    '-p',
    type=Path,
    help='Path to pipeline output directory containing results'
)
@click.option(
    '--barcode',
    '-b',
    is_flag=True,
    help='Group files by barcode regex and merge before collecting [none]'
)
@click.option(
    '--groupby',
    '-g',
    type=str,
    default=None,
    help='Group files by filename regex and merge before collecting [none]'
)
@click.option(
    '--outdir',
    '-o',
    default='collected',
    type=Path,
    help='Path to server collection output directory'
)
def collect(path, barcode, groupby, outdir):

    """ Collect output from pathogen pipeline """

    pp = PathogenPipeline()

    if barcode:
        groupby = r"barcode\d\d"  # EXCLUDING unclassified

    pp.collect_results(path=path, groupby=groupby)

    pp.clean()