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
    help='Glob of VCFs in --vcf_path to construct core genome from  ["*.vcf"]'
)
@click.option(
    '--min_cov',
    default=5,
    type=int,
    help='Minimum coverage in ONT reference mapping to allow a core SNP site across all isolates; disable with -1 [5]'
)
@click.option(
    '--prefix',
    type=str,
    default="core",
    help='Prefix for output alignment and VCF'
)
def hybrid_denovo(vcf_snippy, vcf_ont, vcf_glob, reference, prefix):

    cg = HybridCoreGenome(prefix=prefix, reference=reference)

    cg.parse_snippy_vcf(path=vcf_snippy, vcf_glob=vcf_glob, break_complex=False)  # only consider snps
    cg.parse_ont_vcf(path=vcf_ont, vcf_glob=vcf_glob, min_cov=5)  # min_cov -1 do not assess low cov / gaps in ONT

    cg.call_hybrid_core()
