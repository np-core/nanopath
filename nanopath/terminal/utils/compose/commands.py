import click
import json

from pathlib import Path
from nanopath.utils import ArtificialMixture


@click.command()
@click.option(
    '--composition', '-c', type=Path, help='JSON file, composition configuration'
)
@click.option(
    '--output', '-o', type=Path, default='artificial.fq', help='Output reads file path'
)
@click.option(
    '--reads', '-r', type=int, default=1000, help='Total reads to sample for the mixture'
)
@click.option(
    '--shuffle', '-s', is_flag=True, help='Shuffle output reads'
)
def compose(composition, reads, output, shuffle):

    """ Compose artificial mixtures by sampling from read files """

    with composition.open() as config_file:
        config = json.load(config_file)

    am = ArtificialMixture(
        composition=config,
        reads=reads
    )

    am.compose(fout=output, shuffle=shuffle)
