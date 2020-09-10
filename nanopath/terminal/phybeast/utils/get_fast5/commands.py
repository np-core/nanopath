import click
import pysam

from pathlib import Path
from ont_fast5_api.multi_fast5 import MultiFast5File

@click.command()
@click.option(
    "--fastq", "-fq", type=Path, help="List of Fastq files to extract Fast5 files for"
)
@click.option(
    "--fast5", "-f5", type=Path, help="Directory of Fast5 files to extract reads from"
)
@click.option(
    "--outdir", "-o", type=Path, help="Output directory"
)
def get_fast5(fastq, fast5, outdir):

    """ Get the Fast5 reads for reads in a Fastq """

    outdir.mkdir(parents=True, exist_ok=True)

    with pysam.FastxFile(fastq) as fin:
        with (outdir / 'seqids.txt').open('w') as fout:
            for entry in fin:
                fout.write(str(entry.name) + '\n')