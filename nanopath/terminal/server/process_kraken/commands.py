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
    type=Path,
    default='server_data.json',
    help='Output server data JSON'
)
@click.option(
    '--verbose',
    '-v',
    is_flag=True,
    help='Verbose terminal output'
)
def process_kraken(
    report_file,
    read_file,
    output,
    verbose
):

    """ Process results and generate server data in the metagenome pipeline  """

    kp = KrakenProcessor(
        report=report_file,
        reads=read_file,
        verbose=verbose
    )

    sub_reports, sub_data = kp.summarise_report()
    server_data = kp.get_server_data(*sub_reports, *sub_data)

    with output.open('w') as server_out:
        json.dump(server_data, server_out)