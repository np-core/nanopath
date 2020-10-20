import click

from pathlib import Path
from nanopath.variants import RandomForestFilter

@click.command()
@click.option(
    '--vcf',
    type=Path,
    help=''
)
@click.option(
    '--vcf_dir',
    type=Path,
    help=''
)
@click.option(
    '--model',
    type=Path,
    help=''
)
@click.option(
    '--outdir',
    type=Path,
    default='filtered_calls',
    help=''
)
@click.option(
    '--caller',
    type=str,
    default='clair',
    help=''
)
@click.option(
    '--mask_weak',
    type=float,
    default=0.8,
    help=''
)
def forest_filter(
    vcf,
    vcf_dir,
    outdir,
    model,
    caller,
    mask_weak
):

    if vcf_dir:
        ont_vcf = list(vcf_dir.glob("*.vcf"))
    elif vcf:
        ont_vcf = [vcf]
    else:
        raise ValueError("One of vcf, vcf_dir must be passed as arguments.")

    rff = RandomForestFilter(outdir=outdir)


    # TODO: add multi-threading
    rff.filter_vcf(ont_vcf=ont_vcf, model_file=model, caller=caller, mask_weak=mask_weak)


