import click

from pathlib import Path


@click.command()
@click.option(
    '--path',
    '-p',
    type=Path,
    help='Path to pipeline output directory containing results'
)
@click.option(
    '--outdir',
    '-o',
    default='collected',
    type=Path,
    help='Path to server collection output directory'
)
def run(path, outdir):
    """ Collect output from sepsis pipeline into server data """

    mp = MetagenomePipeline(path=path, outdir=outdir)
    mp.collect()

