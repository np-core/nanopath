import click
import dendropy

from pathlib import Path
from nanopath.pathfinder import phybeast_plot_date_randomisation


@click.command()
@click.option(
    "--file", "-f", type=Path,
    help="Replicate rate file from process DateRandomisationPlot",
)
@click.option(
    "--rate", "-r", type=Path,
    help="True rate file from process MolecularClock; one line.",
)
@click.option(
    "--tree", "-t", type=Path,
    help="Rooted tree file from process MolecularClock;",
)
@click.option(
    "--dates", "-d", type=Path,
    help="Date input file for the workflow (tab-delimited, includes columns: name, date) ",
)
@click.option(
    "--tree_format", default='nexus',
    help="Tree format scheme: nexus (LSD) or newick (TimeTree) ",
)
@click.option(
    "--output", "-o", default="date_randomisation.png",
    help="Output plot file; set extension for format.", type=Path,
)
def plot_date_randomisation(file, output, rate, tree, dates, tree_format):

    """ Plot date randomisation test by Duchene et al. """

    phybeast_plot_date_randomisation(
        replicate_file=file,
        rate_file=rate,
        tree_file=tree,
        date_file=dates,
        tree_format=tree_format,
        output_file=output,
    )


