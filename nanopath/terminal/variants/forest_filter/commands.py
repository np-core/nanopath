import click

from pathlib import Path
from nanopath.variants import ForestClassifier

@click.command()
@click.option(
    '--vcf',
    type=Path,
    help=''
)
@click.option(
    '--stats',
    type=Path,
    help=''
)
@click.option(
    '--model',
    type=Path,
    help=''
)
@click.option(
    '--output',
    type=Path,
    default='filtered.vcf',
    help=''
)
@click.option(
    '--caller',
    type=str,
    default='clair',
    help=''
)
@click.option(
    '--probability',
    type=float,
    default=0,
    help=''
)
@click.option(
    '--mask_weak',
    is_flag=True,
    help=''
)
def forest_filter(
    vcf,
    stats,
    output,
    model,
    caller,
    mask_weak,
    probability
):
    c = ForestClassifier(vcf, stats, output, model, caller, mask_weak, probability)

    pos = c.read_vcf()
    if pos > 0:
        c.classify()
    c.filter()




