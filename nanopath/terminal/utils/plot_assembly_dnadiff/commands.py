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
@click.option(
    "--count_limit", "-l", type=int, help="Variant count limit to exclude tails of low coverage isoaltes", default=0
)
def plot_assembly_dnadiff(
    dir, name, palette, mincov, count_limit
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

    print(dnadiff_cov)

    if mincov > 0:
        dnadiff_cov = dnadiff_cov[dnadiff_cov['mean_coverage_x'] > mincov]


    print(dnadiff_cov)

    snps = dnadiff['snps'].tolist()
    indels = dnadiff['indels'].tolist()
    branch = dnadiff['branch'].tolist()

    data = pandas.DataFrame(
        {
            'variant': ['snp' for _ in snps] + ['indel' for _ in indels],
            'count': snps + indels,
            'branch': branch + branch
         }
    )

    if count_limit > 0:
        data = data[data['count'] <= count_limit]

    sns.despine()
    sns.stripplot(y="count", x="branch", hue="variant", data=data, ax=ax, palette=palette, edgecolor='gray', dodge=True)

    plt.tight_layout()
    plt.ylabel("Count\n")
    plt.xlabel("\nAssembly")

    fig.savefig(f'{name}.pdf')
    fig.savefig(f'{name}.svg')
    fig.savefig(f'{name}.png')

