import click
import pandas
from pathlib import Path
from matplotlib import pyplot as plt
import seaborn as sns
import fnmatch


PATHOGENS = [
    'Burkholderia pseudomallei',
    'Escherichia coli',
    'Klebsiella pneumoniae',
    'Pseudomonas aeruginosa',
    'Salmonella enterica',
    'Staphylococcus aureus'
    'Staphylococcus epidermidis',
    'Streptococcus pneumoniae'
]

CONTAM = [
    'Enterococcus faecium',
    'Prevotella dentalis',
    'Pseudomonas azotoformans',
    'Pseudomonas putida',
    'Pseudomonas poae',
    'Streptococcus oralis'
]


@click.command()
@click.option(
    "--bracken_combined", "-b", type=Path, help="Multi sample Bracken abundance estimates"
)
@click.option(
    "--min_percent", "-m", type=float, default=0.001,
    help="Minimum percent abundance of species (S) in sample to include [0.001]"
)
@click.option(
    "--plot_file", "-p", type=Path, default="multi_abundance.pdf", help="Plot output file [multi_abundance.pdf]"
)
def plot_multi_abundance(
    bracken_combined, plot_file, min_percent
):

    """ Plot a multiple sample Bracken output"""

    nrow, ncol = 1, 1

    fig, ax = plt.subplots(
        nrows=nrow, ncols=ncol, figsize=(
            nrow*7, ncol*4.5
        )
    )

    data = pandas.read_csv(bracken_combined, sep='\t', index_col='name', header=0)

    # Use percentage rather than total reads across samples
    data = data[[c for c in data.columns if 'frac' in c]]
    # Separate viruses
    viruses = []
    for name in data.index.tolist():
        if 'virus' in name.lower():
            viruses.append(data[data.index == name])
    viruses = pandas.concat(viruses)
    data = data.drop(viruses.index.tolist())
    print(viruses)

    human = data[data.index == 'Homo sapiens']
    data = data.drop('Homo sapiens')
    print(human)

    if min_percent > 0:
        keep_idx = []
        for i, row in data.iterrows():
            keep_col = [True for v in row if v >= min_percent]
            if any(keep_col):
                keep_idx.append(row.name)
        data = data[data.index.isin(keep_idx)]

    print(data)

    groups = collapse_taxa(viruses, suffix="virus")

    print(groups)

    groups = collapse_taxa(data, genus=True)

    data.reset_index(level=0, inplace=True)

    data_melt = data.melt(id_vars=['name'], value_name="abundance", var_name="sample")
    data_melt['sample'] = data_melt['sample'].str.replace(".bracken_frac", "")

    sns.scatterplot(data=data, x="gdpPercap", y="lifeExp", size="pop", legend=False, sizes=(20, 2000), ax=ax)

    # plot grid behind markers
    plt.grid(ls="--", zorder=1)
    # take care of long labels
    fig.autofmt_xdate()

    plt.tight_layout()
    plt.ylabel("Count\n")
    plt.xlabel("\nAssembly")

    fig.savefig(f'{plot_file}')


def collapse_taxa(df, genus: bool = False, suffix: str = None) -> list:
    """ Taxa names to collapse are in index of DataFrame """

    if suffix is not None:
        groups = []
        for name in df.index.tolist():
            if suffix in name:
                groups.append(name.split(suffix)[0] + 'virus')
        groups = list(set(groups))
    elif genus:
        groups = []
        for name in df.index.tolist():
            try:
                g = name.split(' ')[0]
            except IndexError:
                continue
            groups.append(g)
    else:
        raise ValueError('Genus or suffix parameters must be set')

    return groups