import click
import pandas
from pathlib import Path
from matplotlib import pyplot as plt
import seaborn as sns


@click.command()
@click.option(
    "--cov_path", "-c", type=Path, help="Locus specific coverage files from BED file loop directory of .cov files"
)
@click.option(
    "--plot_file", "-p", type=Path, default="multi_locus_coverage.pdf",
    help="Plot output file [multi_locus_coverage.pdf]"
)
def plot_loci_cov(
    cov_path, plot_file
):

    """ Plot a multiple enriched loci from an adaptive smpling run """

    cov_files = list(cov_path.glob("*.cov"))

    nrow, ncol = len(cov_files)//2, 2

    fig, ax = plt.subplots(
        nrows=nrow, ncols=ncol, figsize=(
            nrow*7, ncol*4.5
        )
    )

    for locus in cov_files:
        coverage = pandas.read_csv(locus, sep="\t", header=None, names=["locus", 'position', 'coverage'])
        print(coverage)


    sns.scatterplot(data=data, x="gdpPercap", y="lifeExp", size="pop", legend=False, sizes=(20, 2000), ax=ax)

    plt.tight_layout()
    plt.ylabel("Count\n")
    plt.xlabel("\nAssembly")

    fig.savefig(f'{plot_file}')

