import click
import pandas
import shutil
import logging
from pathlib import Path
from nanopath.utils import PoreLogger


@click.command()
@click.option(
    '--directory',
    '-d',
    type=Path,
    help='Path to directory containing barcode read files'
)
@click.option(
    '--file',
    '-f',
    type=Path,
    help='Path to tab-delimited file with columns: barcode, name, (panel)'
)
@click.option(
    '--panel',
    '-p',
    type=str,
    default=None,
    help='Use the column `panel` to select a subset panel from master file'
)
@click.option(
    '--outdir',
    '-o',
    default='renamed',
    type=Path,
    help='Path to output directory for renamed barcode read files'
)
@click.option(
    '--extension',
    '-e',
    default='.fastq',
    type=str,
    help='Fastq barcode file extension [.fastq]'
)
def rename_barcodes(directory, file, panel, outdir, extension):

    """ Rename barcode read files after basecalling"""

    pl = PoreLogger(level=logging.INFO, name="Utils").logger

    pl.info(f'Creating output directory: {outdir}')
    outdir.mkdir(parents=True, exist_ok=True)

    pl.info(f'Reading panel data from file: {file}')
    df = pandas.read_csv(file, sep='\t', header=0)

    if panel:
        pl.info(f'Subsetting panel data to: {panel}')
        df = df.loc[df['panel'].isin([str(panel)]), :]

    for p in directory.rglob("*"):
        if p.stem in df.barcode.tolist():
            row = df.loc[df.barcode == p.stem, :]
            name = row.name.values[0] + f"{extension}"
            pl.info(f'Rename and copy {p.name} to:  {outdir / name}')
            shutil.copy(
                str(p), str(outdir / name)
            )
