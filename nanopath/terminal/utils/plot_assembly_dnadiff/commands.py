import click
import pandas
from pathlib import Path
from matplotlib import pyplot as plt
import seaborn as sns
from collections import Counter


@click.command()
@click.option(
    "--dir", "-d", type=Path, help="Data directory to collected results from assembly pipeline"
)
@click.option(
    "--name", "-n", type=str, help="Plot file output name [dnadiff]", default="dnadiff"
)
@click.option(
    "--palette", "-p", type=str, help="Heatmap palette [Blues]", default="Greens"
)
@click.option(
    "--mincov", "-m", type=int, help="Minimum coverage threshold", default=3
)
def plot_assembly_dnadiff(
    dir, name, palette, mincov
):

    """ Join two dataframes to add selected trait columns """

    nrow, ncol = 1, 1

    fig, ax = plt.subplots(
        nrows=nrow, ncols=ncol, figsize=(
            nrow*7, ncol*4.5
        )
    )

    qc = pandas.read_csv(dir / 'read_qc_all.tsv', sep="\t")

    print(qc)

    dnadiff = pandas.concat(
        [pandas.read_csv(dir / f, sep='\t') for f in ('ont_vs_ref.tsv', 'hybrid_vs_ref.tsv', 'unicycler_vs_ref.tsv')]
    )

    dnadiff_cov = dnadiff.merge(qc, on='name')

    print(plot_assembly_dnadiff)

    if mincov > 0:
        dnadiff_cov = dnadiff_cov[dnadiff_cov['mean_coverage'] > mincov]

    print(dnadiff_cov)

    snps = dnadiff['snps']
    indels = dnadiff['indels']
    branch = dnadiff['branch']
    names = dnadiff['name']

    data = pandas.DataFrame(
        {
            'variants': ['snp' for _ in snps] + ['indel' for _ in indels],
            'count': snps + indels,
            'branch': branch + branch,
            'name': names + names
         }
    )

    sns.despine()
    sns.stripplot(y="count", x="branch", hue="variant", data=data, ax=ax, palette=palette, edgecolor='gray')

    plt.tight_layout()
    plt.ylabel("Outbreak\n")
    plt.xlabel("\nMean coverage")

    fig.savefig(f'{name}.pdf')
    fig.savefig(f'{name}.svg')
    fig.savefig(f'{name}.png')

