import click

from nanopath.processors import KrakenProcessor
from nanopath.utils import get_output_handle
from pysam import FastxFile
from pathlib import Path


@click.command()
@click.option(
    '--fastq',
    '-f',
    type=Path,
    help='Path to Fastq file to filter'
)
@click.option(
    '--report_file',
    '-r',
    type=Path,
    help='Path to tax report file from Kraken'
)
@click.option(
    '--read_file',
    '-d',
    type=Path,
    help='Path to read classification file from Kraken'
)
@click.option(
    '--host',
    '-h',
    default=None,
    type=str,
    help='Host name to filter out'
)
@click.option(
    '--output',
    '-o',
    type=Path,
    default='filtered.fq',
    help='Path to filtered output Fastq'
)
@click.option(
    '--verbose',
    '-v',
    is_flag=True,
    help='Verbose terminal output'
)
def filter_host(
    fastq,
    report_file,
    read_file,
    host,
    output,
    verbose
):

    """ Filter host reads using read classification from Kraken """

    kp = KrakenProcessor(
        report=report_file,
        reads=read_file,
        verbose=verbose,
    )

    if host is not None:
        kp.host = host

    percent, total, reads = kp.get_host(server_json=False)

    if not reads:
        kp.logger.info(
            f'No host reads detected. writing all reads to {output}'
        )
    else:
        reads = reads.read.tolist()
        kp.logger.info(
            f'Detected {len(reads)} reads classified as {kp.host}'
        )
        kp.logger.info(
            f'Writing {kp.total_reads-len(reads)} non-host reads to {output}'
        )

    with FastxFile(fastq) as fin, \
            get_output_handle(output) as fout:
        for read in fin:
            if read.name not in reads:
                fout.write(str(read)+"\n")


