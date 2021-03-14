import click
import pandas
from pathlib import Path
from matplotlib import pyplot as plt
import seaborn as sns
from collections import Counter


@click.command()
@click.option(
    "--data", "-d", type=Path, help="Data file long format with header: name,location,mlst"
)
@click.option(
    "--name", "-n", type=str, help="Plot file output name [mlst]", default="heatmap"
)
@click.option(
    "--palette", "-p", type=str, help="Heatmap palette [Blues]", default="Blues"
)
@click.option(
    "--column", "-c", type=str, help="Column to subset the data", default="cluster"
)
def plot_coverage_summary(
    data, name, palette, column
):

    """ Join two dataframes to add selected trait columns """

    nrow, ncol = 1, 1

    fig, ax = plt.subplots(
        nrows=nrow, ncols=ncol, figsize=(
            nrow*7, ncol*4.5
        )
    )

    df = pandas.read_csv(data, sep="\t")

    sns.stripplot(y=column, x="mean_coverage", data=df, ax=ax, palette=palette)

    plt.tight_layout()

    fig.savefig(f'{name}.pdf')
    fig.savefig(f'{name}.svg')
    fig.savefig(f'{name}.png')

