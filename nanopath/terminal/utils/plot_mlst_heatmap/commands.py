import click
import pandas
from pathlib import Path
from matplotlib import pyplot as plt
import seaborn as sns
from collections import Counter
import re
import numpy as np


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
def plot_mlst_heatmap(
    data, name, palette
):

    """ Join two dataframes to add selected trait columns """

    nrow, ncol = 1, 1

    fig, ax = plt.subplots(
        nrows=nrow, ncols=ncol, figsize=(
            nrow*7, ncol*4.5
        )
    )

    df = pandas.read_csv(data, sep="\t")
    if 'name' in df.columns:
        df.drop(columns='name', inplace=True)

    counters = df.apply(Counter, axis='columns').value_counts()

    rows = []
    for i, c in enumerate(counters):
        meta = list(
            counters.index[i].keys()
        )
        row = [meta[0], meta[1], c]
        rows.append(row)

    df = pandas.DataFrame(rows, columns=["location", 'mlst', 'count'])
    df = df.dropna()

    df = df[df['mlst'] != 'ST-']  # exclude MLST SLVs
    df['mlst'] = [int(v.strip("ST")) for v in df['mlst']]
    df = df.sort_values('mlst')

    sns.heatmap(
        df.pivot("mlst", "location", "count"), cmap=palette, ax=ax, linewidths=1, annot=True, cbar=False
    )
    ax.set_facecolor('#eeeeee')
    ax.yaxis.label.set_color('gray')
    ax.tick_params(axis='x', colors='gray')
    ax.tick_params(axis='y', colors='gray')

    plt.ylabel("Sequence type (ST)\n")
    plt.xlabel("")
    plt.yticks(rotation=0)

    plt.tight_layout()
    fig.savefig(f'{name}.pdf')
    fig.savefig(f'{name}.svg')
    fig.savefig(f'{name}.png')

