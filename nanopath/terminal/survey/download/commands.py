import click
import pandas

from nanopath.surveillance import Survey
from nanopath.surveillance import MiniAspera

from pathlib import Path


@click.command()
@click.option(
    '--file', '-f', type=Path, default=None, required=True,
    help='Output file from ENA query [none]'
)
@click.option(
    '--outdir', '-o', type=Path, default="pf-download",
    help='Output directory for read files'
)
@click.option(
    '--batch', '-b', type=int, default=0,
    help='Batch large search results into their own output directories'
)
@click.option(
    '--filter', '-fi', type=str, default=None,
    help='Pandas filter string on the return fields of the ENA query '
         'for example: "coverage < 700 & coverage > 50" '
)
@click.option(
    '--ftp', '-f', is_flag=True,
    help='Force download from FTP instead of Aspera (slow)'
)
@click.option(
    '--raw', '-r', is_flag=True,
    help='Read raw table query file [false]'
)
@click.option(
    '--sub', '-s', is_flag=True,
    help='Submitted files instead of fastq from ENA [false]'
)
def download(outdir, batch, file, filter, ftp, raw, sub):
    """ Download sequence read data from ENA """

    survey = Survey(outdir=outdir)
    survey.read_query_file(file=file)

    if not survey.query.empty:
        print(f"Total download size: {sum(survey.query['size']) / 1024:.2f} GB [n = {len(survey.query)}]")

    if raw:
        ascp = MiniAspera()
        ascp.download_raw(
            data=survey.query,
            outdir=outdir,
            submitted=sub,
            ftp=ftp
        )
    else:
        if filter is not None:
            survey.filter_query(filter)

        if batch > 0:
            batches = survey.batch(query_result=survey.query, batch_size=batch)
            batches = survey.batch_output(batches, outdir=outdir)
        else:
            batches = [(
                outdir, survey.query
            )]

        for batch_path, batch_data in batches:
            ascp = MiniAspera()
            ascp.download_batch(
                data=batch_data,
                outdir=batch_path,
                ftp=ftp
            )
