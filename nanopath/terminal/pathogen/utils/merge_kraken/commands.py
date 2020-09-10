import click

from pathlib import Path
from nanopath.tools.kraken import merge_reports


@click.command()
@click.option(
    "--reports", "-r", help="Report files", multiple=True, type=Path
)
@click.option(
    "--output", "-o", help="Merged report file", type=Path,
)
def merge_kraken(reports, output):

    """ Merge Kraken-style reports"""

    merge_reports(r_files=reports, output=str(output))


