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
    'Staphylococcus aureus',
    'Staphylococcus epidermidis',
    'Streptococcus pneumoniae',
]

CONTAM = [
    'Enterococcus faecium',
    'Prevotella dentalis',
    'Pseudomonas azotoformans',
    'Pseudomonas putida',
    'Pseudomonas poae',
    'Streptococcus oralis',
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
    data = data[[c for c in data.columns if 'bracken_frac' in c]]
    data.columns = [d.replace(".bracken_frac", "") for d in data.columns]

    # Separate viruses
    viruses = []
    for name in data.index.tolist():
        if 'virus' in name.lower():
            viruses.append(data[data.index == name])
    viruses = pandas.concat(viruses)
    data = data.drop(viruses.index.tolist())

    human = data[data.index == 'Homo sapiens']
    data = data.drop('Homo sapiens')
    print(human)

    pathogens = data[data.index.isin(PATHOGENS)].sort_index()
    data = data.drop(PATHOGENS)

    contams = data[data.index.isin(CONTAM)].sort_index()
    data = data.drop(CONTAM)

    print(pathogens)
    print(contams)

    if min_percent > 0:
        keep_idx = []
        for i, row in data.iterrows():
            keep_col = [True for v in row if v >= min_percent]
            if any(keep_col):
                keep_idx.append(row.name)
        data = data[data.index.isin(keep_idx)]

    viruses_collapsed = collapse_taxa(viruses, suffix="virus")

    print(viruses_collapsed)

    other_collapsed = collapse_taxa(data, genus=True)

    print(other_collapsed)

    # data.reset_index(level=0, inplace=True)
    #
    # data_melt = data.melt(id_vars=['name'], value_name="abundance", var_name="sample")
    #
    # sns.scatterplot(data=data, x="gdpPercap", y="lifeExp", size="pop", legend=False, sizes=(20, 2000), ax=ax)
    #
    # # plot grid behind markers
    # plt.grid(ls="--", zorder=1)
    # # take care of long labels
    # fig.autofmt_xdate()
    #
    # plt.tight_layout()
    # plt.ylabel("Count\n")
    # plt.xlabel("\nAssembly")
    #
    # fig.savefig(f'{plot_file}')


def collapse_taxa(df: pandas.DataFrame, genus: bool = False, suffix: str = None) -> pandas.DataFrame:
    """ Taxa names to collapse are in index of DataFrame """

    if suffix is not None:
        group_names = []
        for name in df.index.tolist():
            if suffix in name:
                group_names.append(name.split(suffix)[0] + 'virus')
    elif genus:
        group_names = []
        for name in df.index.tolist():
            try:
                g = name.split(' ')[0]
            except IndexError:
                continue
            group_names.append(g)
    else:
        raise ValueError('Genus or suffix parameters must be set')

    grouped = []
    df.index = group_names
    for group, gdata in df.groupby(df.index):
        # print(f'Collapsing species in genus: {group}')
        d = gdata.apply(sum, axis=0)
        d.name = f"{group} spp." if genus else f"{group}"
        grouped.append(d)

    grouped = pandas.DataFrame(grouped)

    return grouped