import click

from pathlib import Path
from nanopath.variants import HybridCoreGenome


@click.command()
@click.option(
    '--vcf_snippy',
    default=Path.cwd() / 'snippy',
    type=Path,
    help='Path to .vcf / .aligned.fa from background population Snippy calls [snippy/]'
)
@click.option(
    '--vcf_ont',
    default=Path.cwd() / 'ont',
    type=Path,
    help='Path to .vcf from RandomForest filtered ONT calls [ont/]'
)
@click.option(
    '--vcf_glob',
    default="*.vcf",
    type=str,
    help='Glob of VCFs in --vcf_path to construct core genome from  ["*.vcf"]'
)
@click.option(
    '--prefix',
    type=str,
    default="core",
    help='Prefix for output alignment and VCF'
)
def snp_core(vcf_snippy, vcf_ont, vcf_glob, prefix):

    cg = HybridCoreGenome()

    # cg.parse_snippy_vcf(path=vcf_snippy, vcf_glob=vcf_glob, break_complex=False)  # only consider snps
    cg.parse_ont_vcf(path=vcf_ont, vcf_glob=vcf_glob, min_cov=10)
    #
    # cg.find_core_genome_sites()
