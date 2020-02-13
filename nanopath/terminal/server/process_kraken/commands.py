import click
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
    type=Path,
    default='plot.png',
    help='Output plot file'
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
    help='Color palette for central donut plot.'
)
@click.option(
    '--verbose',
    '-v',
    is_flag=True,
    help='Verbose terminal output'
)
def process_kraken(report_file, read_file, output, compose, palette, verbose):

    """ Process results and generate server data in the metagenome pipeline  """

    kp = KrakenProcessor(report=report_file, reads=read_file, verbose=verbose)

    if compose is not None:
        summary, summary_props, negatives = \
            kp.assess_composition(compose_file=compose)
        kp.create_summary_plot(summary, summary_props, negatives)