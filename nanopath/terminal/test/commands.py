import click

from pathlib import Path
from nanopath.processors import AssemblyProcessor


@click.command()
@click.option(
    '--fasta',
    '-f',
    type=Path,
    help='Path to assembled sequences in Fasta'
)
@click.option(
    '--info',
    '-i',
    type=Path,
    default=None,
    help='Path to assembly info file (Flye)'
)
@click.option(
    '--verbose',
    '-v',
    is_flag=True,
    help='Verbose terminal output'
)
def test(fasta, info, verbose):

    """ Task for testing command line functions """

    ap = AssemblyProcessor(
        fasta=fasta,
        info_file=info,
        verbose=verbose
    )

    ap.get_assembly_data()