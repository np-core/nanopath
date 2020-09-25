import click

from pathlib import Path


@click.command()
@click.option(
    "--model", "-m", required=False, type=str, default="bdss",
    help="Model to select configuration template: bdss, bdsc, mtbd, cosky",
)
@click.option(
    "--output", "-o", required=False, type=Path, default=Path("bdss.yaml"),
    help="Output configuration file (YAML)",
)
def config(model, output):

    """ Get a Beastling configuration template YAML """

