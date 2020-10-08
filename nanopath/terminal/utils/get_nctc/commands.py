import click
from pathlib import Path
from nanopath.tools.nctc import NCTC3000


@click.command()
@click.option(
    "--outdir", "-d", default=".", help="Output directory of NCTC3000 DB", type=Path,
)
@click.option(
    "--species", "-s", default="Staphylococcus aureus", help="Species to obtain reference genomes for", type=str,
)
@click.option(
    "--chromosomes", "-c", is_flag=True, help=""
)
@click.option(
    "--force", "-f", is_flag=True, help=""
)
@click.option(
    "--verbose", "-v", is_flag=True, help=""
)
def get_nctc(outdir, species, chromosomes, force, verbose):

    """ Get a NCTC collection """

    nctc = NCTC3000(path=outdir, species=species, verbose=verbose)
    nctc.make(strict=chromosomes, force=force)
