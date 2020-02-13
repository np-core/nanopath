import click

from pathlib import Path
from nanopath.utils import HostDecontaminator


@click.command()
@click.option(
    '--fastq', '-f', type=Path, help='Input Fastq file, mixed sample (host, microbial)'
)
def decon(fastq):

    """ Decontaminate a mixed sample by removing host reads """

    hd = HostDecontaminator(fastq=fastq)
    hd.decontaminate(method='minikraken2', db=Path.home() / 'resources' / 'minikraken2')
