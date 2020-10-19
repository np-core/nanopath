import click

from pathlib import Path
from nanopath.variants import CoreGenome


@click.command()
@click.option(
    '--vcf_path',
    default=Path.cwd(),
    type=Path,
    help='Path to VCFs to construct core genome from [$CWD]'
)
@click.option(
    '--vcf_glob',
    default="*.vcf",
    type=str,
    help='Glob of VCFs in --vcf_path to construct core genome from  ["*.vcf"]'
)
@click.option(
    '--allow_missing',
    default=0.1,
    type=float,
    help='Allow for missing SNP rate across isolates in the core SNP alignment, include missing as --missing_char'
)
@click.option(
    '--missing_char',
    default="N",
    type=str,
    help='Missing character to replace for gaps in alignment due to --allow_missing ["N"]'
)
@click.option(
    '--snp_limit',
    default=0,
    type=int,
    help='Set a SNP limit filter as parsed from the VCF [0]'
)
@click.option(
    '--prefix',
    type=str,
    default="core",
    help='Prefix for output alignment and VCF'
)
def snp_core(vcf_path, vcf_glob, prefix, allow_missing, missing_char, snp_limit):

    cg = CoreGenome()
    cg.parse_snp_vcf(path=vcf_path, vcf_glob=vcf_glob)
    cg.core_genome(allow_missing=allow_missing, missing_char=missing_char, snp_limit=snp_limit)
