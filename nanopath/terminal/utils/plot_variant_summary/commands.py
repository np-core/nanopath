import click
import pandas
from pathlib import Path
from matplotlib import ticker as mticker
from matplotlib import pyplot as plt
import seaborn as sns
from numpy import log10, linspace


@click.command()
@click.option(
    "--data", "-d", type=str, help="Data file with header, merged assembly and variant caller evaluations"
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
    "--panels", "-p", type=str, default="model", help="Column in dfs to use as paneling for rows [model]"
)
@click.option(
    "--variables", "-v", type=str, default="default", help="One or multiple y-axis variables to plot [apr=accuracy,precision,recall]"
)
@click.option(
    "--sep", "-si", default="\t", help="Input delimiter [\\t]", type=str,
)
@click.option(
    "--title_ref", "-tr", default="ST93", help="Title of reference used for variant calling [ST93]", type=str,
)
@click.option(
    "--legend_size", "-ls", default=8, help="Legend size", type=int,
)
@click.option(
    "--palettes", "-pa", default="Blues,Purples", help="Color palettes for panels in each row [Blues,Purples]", type=str,
)
@click.option(
    "--locs", "-l", default="best,upper_left,best,upper_left,upper_left,upper_left", help="Position of legends for plots", type=str,
)
@click.option(
    "--min_column", "-mc", default=None, help="Column for lower threshold subset of data", type=str,
)
@click.option(
    "--min_value", "-mv", default=0, help="Minimum column value for lower threshold of data", type=float,
)
def plot_variant_summary(
        data, output, column, binary, sep, variables, palettes,
        title_ref, locs, min_column, min_value, legend_size, panels
):

    """ Join two dataframes to add selected trait columns """

    palettes = [s.strip() for s in palettes.split(',')]

    attr = dict(
        locs=[lo.replace("_", " ").strip() for lo in locs.split(",")]
    )

    if variables == "default":
        varis = [
            ["accuracy", "precision", "recall"],
            ["true_positives", "true_negatives", "false_positives", "false_negatives"]
        ]
        varis_order = ["apr", "raw"]
    else:
        raise ValueError

    data_frame = pandas.read_csv(data, sep=sep, header=0)

    panel_rows = sorted(list(
        set(data_frame[panels].unique().tolist())
    ))

    new_attr = {}
    for a, attr_data in attr.items():
        if len(attr_data) > len(panel_rows)*len(varis):
            new_attr[a] = attr_data[:len(panel_rows)]
        elif len(attr_data) == len(panel_rows)*len(varis):
            new_attr[a] = attr_data
        else:
            new_attr[a] = [attr_data[0] for _ in range(len(panel_rows)*len(varis))]

    nrow = len(panel_rows)

    fig, axes = plt.subplots(
        nrows=nrow, ncols=len(varis), figsize=(
            len(varis) * 7, nrow * 4.5
        )
    )

    t = 0
    for r, _ in enumerate(panel_rows):
        for c, var in enumerate(varis):

            # subset data to correct panel row value (e.g. model)
            df = data_frame.loc[data_frame[panels] == panel_rows[r], :]

            if min_column:
                df = df[df[min_column] > min_value]

            if binary:
                df[column] = [binary if _ == binary else "Other" for _ in df[column]]
                vc = df[column].value_counts()
                df[column] = [f"{v} (n = {vc[v]})" for v in df[column]]

            dfm = df.melt(id_vars=['name', column], value_vars=var).sort_values(['name', 'variable'])

            if varis_order[c] == "apr":
                dfm['value'] = dfm['value']*100
            else:
                dfm['value'] = [log10(d) for d in dfm['value'].tolist()]

            fig.subplots_adjust(hspace=1.5)

            ax = axes[c] if nrow == 1 else axes[r][c]
            pal = palettes[c]

            sns.violinplot(y='value', x=column, hue='variable', data=dfm, color="0.8", ax=ax, palette=pal)
            sns.stripplot(y='value', x=column, hue='variable', data=dfm, jitter=True,
                          zorder=1, palette="Greys", linewidth=1, ax=ax, dodge=True)

            if c == 0:
                ax.set_ylabel(
                    "Illumina reference SNPs (%)\n"
                    if varis_order[c] == "apr" else "ONT vs Illumina reference SNPs (n)"
                )
            else:
                ax.set_ylabel("\n")

            if r == len(panel_rows)-1:
                ax.set_xlabel(f"\n{column.upper()}")
            else:
                ax.set_xlabel(f"\n")

            handles, labels = ax.get_legend_handles_labels()

            ax.legend(
                handles[0:len(var)],
                labels[0:len(var)],
                title="", loc=new_attr["locs"][t], prop={'size': legend_size}
            )
            ax.set_title(f'\n{title_ref} reference: {panel_rows[r]}')

            if varis_order[c] != "apr":
                ax.yaxis.set_major_formatter(mticker.StrMethodFormatter("$10^{{{x:.0f}}}$"))
                ax.yaxis.set_ticks(
                    [log10(x) for p in range(1, 6) for x in linspace(10 ** p, 10 ** (p + 1), 10)], minor=True
                )
            t += 1

    plt.tight_layout()
    fig.savefig(f'{output}')