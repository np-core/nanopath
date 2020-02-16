import click
import json

from pathlib import Path
from nanopath.processors import KrakenProcessor


@click.command()
@click.option(
    '--report_file',
    '-r',
    type=Path,
    help='Path to tax report files'
)
@click.option(
    '--read_file',
    '-f',
    type=Path,
    help='Path to read classification file'
)
@click.option(
    '--output',
    '-o',
    type=str,
    default='assess',
    help='Output prefix for composition assessment files [assess]'
)
@click.option(
    '--compose',
    '-c',
    type=Path,
    default=None,
    help='Composition JSON to run mixture assessment'
)
@click.option(
    '--palette',
    '-p',
    type=str,
    default='Greens',
    help='Color palette for central donut plot'
)
@click.option(
    '--verbose',
    '-v',
    is_flag=True,
    help='Verbose terminal output'
)
def assess(report_file, read_file, output, compose, palette, verbose):

    """ Compose artificial mixtures by sampling from read files """

    kp = KrakenProcessor(
        report=report_file,
        reads=read_file,
        verbose=verbose
    )

    summary, summary_props, negatives = \
        kp.assess_composition(compose_file=compose)

    kp.create_summary_plot(
        summary, summary_props, negatives, palette=palette, prefix=output
    )
