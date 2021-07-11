import click
import pandas
from pathlib import Path
from matplotlib import pyplot as plt
import seaborn as sns


@click.command()
@click.option(
    "--bracken_combined", "-b", type=Path, help="Multi sample Bracken abundance estimates"
)
@click.option(
    "--plot_file", "-p", type=Path, default="multi_abundance.pdf", help="Plot output file [multi_abundance.pdf]"
)
def plot_multi_abundance(
    bracken_combined, plot_file
):

    """ Plot a multiple sample Bracken output"""

    nrow, ncol = 1, 1

    fig, ax = plt.subplots(
        nrows=nrow, ncols=ncol, figsize=(
            nrow*7, ncol*4.5
        )
    )

    data = pandas.read_csv(bracken_combined, sep='\t', index_col='name', header=0)

    print(data.name.tolist())


    sns.scatterplot(data=data, x="gdpPercap", y="lifeExp", size="pop", legend=False, sizes=(20, 2000), ax=ax)

    # plot grid behind markers
    plt.grid(ls="--", zorder=1)
    # take care of long labels
    fig.autofmt_xdate()

    plt.tight_layout()
    plt.ylabel("Count\n")
    plt.xlabel("\nAssembly")

    fig.savefig(f'{plot_file}')

