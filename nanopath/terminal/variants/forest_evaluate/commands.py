import click

from pathlib import Path
from nanopath.variants import RandomForestFilter


@click.command()
@click.option(
    '--model',
    type=Path,
    help='Path to trained Random Forest model [.sav] (from: forest-train)'
)
@click.option(
    '--vcf_snippy',
    type=Path,
    help='Path to reference Snippy variant file (from: np-core/np-variants --workflow core)'
)
@click.option(
    '--vcf_ont',
    type=Path,
    help='Path to ONT (Medaka or Clair) variant file (from: np-core/np-variants --workflow denovo)'
)
@click.option(
    '--stats_ont',
    type=Path,
    help='Path to ONT (Medaka or Clair) Pysamstats file (from: np-core/np-variants --workflow denovo)'
)
@click.option(
    '--dir_snippy',
    type=Path,
    help='Path to directory with VCFs from Snippy reference variant calls (from: np-core/np-variants)'
)
@click.option(
    '--dir_ont',
    type=Path,
    help='Path to directory with VCFs from Medaka or Clair (from: np-core/np-variants)'
)
@click.option(
    '--outdir',
    type=Path,
    help='Path to output directory'
)
@click.option(
    '--caller',
    type=str,
    default="clair",
    help='ONT variant caller, one of: medaka, clair'
)
@click.option(
    '--prefix',
    type=str,
    default="prefix",
    help='Output evaluation prefix [prefix]'
)
@click.option(
    '--mask_weak',
    type=float,
    default=0.,
    help='Mask weak position with < base support than <0. - 1.> [0.8]'
)
@click.option(
    '--break_complex',
    is_flag=True,
    help='Break complex SNP variants from Freebayes/Snippy in reference [false]'
)
def forest_evaluate(dir_snippy, dir_ont, outdir, caller, model, prefix, mask_weak, break_complex, vcf_snippy, vcf_ont, stats_ont):

    ft = RandomForestFilter(outdir=outdir)

    ft.evaluate_model(
        dir_snippy=dir_snippy,
        dir_ont=dir_ont,
        vcf_snippy=vcf_snippy,
        vcf_ont=vcf_ont,
        stats_ont=stats_ont,
        model_file=model,
        caller=caller,
        break_complex=break_complex,
        prefix=prefix,
        mask_weak=mask_weak
    )



