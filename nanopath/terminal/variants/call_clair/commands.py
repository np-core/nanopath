import click

from pathlib import Path
from nanopath.processors import ClairSample


@click.command()
@click.option(
    '--clair_directory',
    type=Path,
    help='Path to directory containing output directories with VCFs from Clair'
)
@click.option(
    '--snippy_directory',
    type=Path,
    help='Path to directory containing matching output VCFs from Snippy'
)
@click.option(
    '--reference',
    type=Path,
    help='Reference genome used in Snippy and Clair [none]'
)
@click.option(
    '--prefix',
    type=str,
    default='clair_',
    help='Prefix for main output files [megalodon]'
)
@click.option(
    '--outdir',
    type=Path,
    default='clair_core',
    help='Path to output directory [medaka_core]'
)
@click.option(
    '--include_reference',
    is_flag=True,
    help='Include the reference sequence in variant alignments [false]'
)
def call_clair(
    clair_directory,
    snippy_directory,
    reference,
    include_reference,
    outdir,
    prefix,
):

    """ Detect core SNPs from Medaka and Snippy """

    for vcf in clair_directory.glob("*.vcf"):
        print(f"Processing: {vcf.name}")
        cs = ClairSample(vcf=vcf)
        cs.filter(min_quality=30)
        cs.write_vcf(Path(f"{vcf.name}.out.vcf"))
