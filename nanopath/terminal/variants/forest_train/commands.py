import click

from pathlib import Path
from nanopath.variants import RandomForestFilter


@click.command()
@click.option(
    '--dir_snippy',
    type=Path,
    help='Path to directory with VCFs from Snippy reference variant calls (np-core)'
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
def forest_train(dir_snippy, dir_ont, outdir, caller):

    ft = RandomForestFilter(outdir=outdir)

    ft.prepare_training_data(
        dir_snippy=dir_snippy,
        dir_ont=dir_ont,
        caller=caller,
        break_complex=True
    )
    ft.train_models(
        test_size=0.3,
        model_prefix="test_model"
    )



