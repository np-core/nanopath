import click

from pathlib import Path
from nanopath.variants import RandomForestFilter


@click.command()
@click.option(
    '--dir_snippy',
    type=Path,
    help='Path to directory with VCFs from Snippy reference variant calls .ref.vcf (np-core)'
)
@click.option(
    '--dir_ont',
    type=Path,
    help='Path to directory with VCFs from coverage subset Medaka or Clair (np-core)'
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
    default="model",
    help='Model name prefix for output to {prefix}_{composite,qual}.sav [model]'
)
@click.option(
    '--break_complex',
    is_flag=True,
    help='Break complex SNP variants from Freebayes/Snippy in reference [false]'
)
@click.option(
    '--test_size',
    is_flag=True,
    type=float,
    help='Test set size in test-train split [0.3]'
)
def forest_train(dir_snippy, dir_ont, outdir, caller, test_size, break_complex, prefix):

    ft = RandomForestFilter(outdir=outdir)

    ft.prepare_training_data(
        dir_snippy=dir_snippy,
        dir_ont=dir_ont,
        caller=caller,
        break_complex=break_complex
    )
    ft.train_models(
        test_size=test_size,
        model_prefix=prefix
    )



