import click
import pandas
from pathlib import Path
from matplotlib import ticker as mticker
from matplotlib import pyplot as plt
import seaborn as sns
from numpy import log10, linspace


@click.command()
@click.option(
    "--dfs", "-d", type=str, help="Data file with header, merged assembly and variant caller evaluations"
)
@click.option(
    "--output", "-o", type=Path, help="Plot file output [pub_mlst_variants.png]"
)
@click.option(
    "--column", "-c", type=str, default="illumina mlst", help="Get header name tp plot a binary violin by categorical data as rows over all unique values [mlst]"
)
@click.option(
    "--binary", "-b", type=str, default="ST93", help="Get one value and plot a binary violin by categorical data vs others from column [none]"
)
@click.option(
    "--variables", "-v", type=str, default="apr", help="One or multiple y-axis variables to plot [apr=accuracy,precision,recall]"
)
@click.option(
    "--sep", "-si", default="\t", help="Input delimiter [\\t]", type=str,
)
@click.option(
    "--title_ref", "-tr", default="ST93", help="Title of reference used for variant calling [ST93]", type=str,
)
@click.option(
    "--title_caller", "-tc", default="Clair SNP calls raw,Clair", help="Variant caller used for RandomForest training and variant calling", type=str,
)
@click.option(
    "--title_rff", "-tf", default=" , + ST93-MRSA RandomForest", help="RandomForest classifier model used for SNP validation", type=str,
)
@click.option(
    "--legend_size", "-ls", default=8, help="Legend size", type=int,
)
@click.option(
    "--palettes", "-p", default="YlOrRd,Greens", help="Color palettes for panels in each row [YlOrRd,Greens]", type=str,
)
@click.option(
    "--locs", "-l", default="upper_left,lower_left", help="Position of legends for plots", type=str,
)
@click.option(
    "--min_column", "-mc", default=None, help="Column for lower threshold subset of data", type=str,
)
@click.option(
    "--min_value", "-mv", default=0, help="Minimum column value for lower threshold of data", type=float,
)
def plot_variant_summary(
        dfs, output, column, binary, sep, variables, palettes, title_rff,
        title_ref, title_caller, locs, min_column, min_value, legend_size
):

    """ Join two dataframes to add selected trait columns """

    attr = dict(
        title_caller=[s.strip() for s in title_caller.split(',')],
        title_rff=[s.strip() for s in title_rff.split(',')],
        palettes=[s.strip() for s in palettes.split(',')],
        locs=[lo.replace("_", " ").strip() for lo in locs.split(",")]
    )

    if variables == "apr":
        var = ["accuracy", "precision", "recall"]
    else:
        var = ["true_positives", "true_negatives", "false_positives", "false_negatives"]

    dfs = dfs.split(",")

    new_attr = {}
    for a, attr_data in attr.items():
        if len(attr_data) > len(dfs):
            new_attr[a] = attr_data[:len(dfs)]
        elif len(attr_data) == len(dfs):
            new_attr[a] = attr_data
        else:
            raise ValueError(f"Too few values provided for {a}")

    print(dfs)

    if len(dfs) == 2:
        nrow = 1
    elif len(dfs) % 2 == 0:
        nrow = len(dfs) // 2
    else:
        nrow = (len(dfs) // 2)+1

    fig, axes = plt.subplots(
        nrows=nrow, ncols=2, sharey='row', figsize=(
            2 * 7, nrow * 4.5
        )
    )

    def chunks(iterable, size):
        from itertools import chain, islice
        iterator = iter(iterable)
        for first in iterator:
            yield list(chain([first], islice(iterator, size - 1)))

    dfs = [b for b in chunks(dfs, 2)]
    total = 0
    for r, batch in enumerate(dfs):
        for c, df in enumerate(batch):
            df = pandas.read_csv(df, sep=sep, header=0)

            if min_column:
                df = df[df[min_column] > min_value]

            if binary:
                df[column] = [binary if _ == binary else "Other" for _ in df[column]]
                vc = df[column].value_counts()
                df[column] = [f"{v} (n = {vc[v]})" for v in df[column]]

            dfm = df.melt(id_vars=['name', column], value_vars=var).sort_values(['name', 'variable'])

            if variables == "apr":
                dfm['value'] = dfm['value']*100
                print(dfm['value'])
            else:
                dfm['value'] = [log10(d) for d in dfm['value'].tolist()]

            fig.subplots_adjust(hspace=1.5)

            ax = axes[c] if nrow == 1 else axes[r][c]

            sns.violinplot(y='value', x=column, hue='variable', data=dfm, color="0.8", ax=ax,
                           palette=new_attr["palettes"][total])
            sns.stripplot(y='value', x=column, hue='variable', data=dfm, jitter=True,
                          zorder=1, palette="Greys", linewidth=1, ax=ax, dodge=True)

            if c == 0:
                ax.set_ylabel(
                    "Illumina reference SNPs (%)\n"
                    if variables == "apr" else "ONT vs Illumina reference SNPs (n)"
                )
            else:
                ax.set_ylabel("\n")

            if r == len(dfs)-1:
                ax.set_xlabel(f"\n{column.upper()}")
            else:
                ax.set_xlabel(f"\n")

            handles, labels = ax.get_legend_handles_labels()

            ax.legend(
                handles[0:len(var)],
                labels[0:len(var)],
                title="", loc=attr["locs"][total], prop={'size': legend_size}
            )
            ax.set_title(f'\n{title_ref} reference: {new_attr["title_caller"][total]} {new_attr["title_rff"][total]}')

            if variables != "apr":
                ax.yaxis.set_major_formatter(mticker.StrMethodFormatter("$10^{{{x:.0f}}}$"))
                ax.yaxis.set_ticks(
                    [log10(x) for p in range(1, 6) for x in linspace(10 ** p, 10 ** (p + 1), 10)], minor=True
                )
            total += 1

    plt.tight_layout()
    fig.savefig(f'{output}')