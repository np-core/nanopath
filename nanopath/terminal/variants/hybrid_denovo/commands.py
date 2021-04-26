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
    '--reference',
    default=None,
    type=Path,
    help='Reference genome for full SNP alignment [none]'
)
@click.option(
    '--vcf_glob',
    default="*.vcf",
    type=str,
    help='Glob of VCFs in --vcf_snippy / --vcf_ont to construct core genome from  ["*.vcf"]'
)
@click.option(
    '--min_cov',
    default=2,
    type=int,
    help='Minimum coverage in ONT reference mapping to allow a core SNP site across all isolates; disable with -1 [2]'
)
@click.option(
    '--allow_missing',
    default=0,
    type=float,
    help='Allow sites to be included if they do not exceed this proportional threshold [0]'
)
@click.option(
    '--break_complex',
    is_flag=True,
    help='Break complex variant types in Snippy (MNP, COMPLEX) and include in core alignment'
)
@click.option(
    '--keep_all',
    is_flag=True,
    help='Do not exclude any SNPs based on gaps in reference regions from .aligned.fa (non-core sites) [false]'
)
@click.option(
    '--include_reference',
    is_flag=True,
    help='Include reference in final alignment [false]'
)
@click.option(
    '--prefix',
    type=str,
    default="core",
    help='Prefix for output alignment and VCF'
)
def hybrid_denovo(vcf_snippy, vcf_ont, vcf_glob, reference, prefix, min_cov, break_complex, allow_missing, include_reference, keep_all):

    cg = HybridCoreGenome(prefix=prefix, reference=reference)

    cg.parse_snippy_vcf(path=vcf_snippy, vcf_glob=vcf_glob, break_complex=break_complex)
    cg.parse_ont_vcf(path=vcf_ont, vcf_glob=vcf_glob, min_cov=min_cov)

    cg.call_hybrid_core(include_reference=include_reference, keep_all=keep_all)
