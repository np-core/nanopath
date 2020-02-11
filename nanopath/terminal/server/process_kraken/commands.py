import click
import pandas
import json

from datetime import datetime
from pathlib import Path
from matplotlib import pyplot as plt
import numpy as np
import seaborn as sns
from nanopath.pipelines import Metagenome

@click.command()
@click.option(
    '--report', '-r', type=str, help='Path or file glob to tax report files'
)
@click.option(
    '--prefix', '-p', type=str, help='Output prefix for plot file.'
)
@click.option(
    '--top', '-t', type=int, default=10,
    help='Show top taxonomic levels in plots [10]'
)
@click.option(
    '--color', '-c', type=str, default='Greens',
    help='Color palette for central donut plot.'
)
@click.option(
    '--title', '-t', type=str, default=None,
    help='Row titles for center plot, comma separated string.'
)
@click.option(
    '--sub', '-s', is_flag=True,
    help='Add subplot titles for each column.'
)
def process_kraken(report, prefix, top, color, title, sub):

    """ Process results and generate server data in the metagenome pipeline  """

    mg = Metagenome()
    mg.set_meta(
        uuid=None,
        user=None,
        name=None,
        submitted=None,
        host_removal=''
    )
    mg.process_kraken(report=report)
    mg.get_server_data(
        fout=Path(f'{prefix}.json')
    )
    mg.plot_kraken_summary(
        plot_file=Path(f'{prefix}.png'),
        palette=color,
        top_minor=top,
        title=title,
        subtitles=sub
    )


def plot_overview(df, ax, level, color=None):

    target = sns.color_palette(color, 2)[-1]

    overview_data, overview_labels = [], []

    human = df[df['taxonomy'] == 'Homo sapiens']

    if not human.empty:
        overview_data.append(
            float(human['percent'])
        )
        overview_labels.append(f'Human [{float(human["percent"])}%]')

    unclassified = df[df['taxonomy'] == 'unclassified']

    if not unclassified.empty:
        overview_data.append(
            float(unclassified['percent'])
        )
        overview_labels.append(
            f'Unclassified [{float(unclassified["percent"])}%]'
        )

    microbial = df[(~df.taxonomy.isin(
        ['Homo sapiens', 'Homo', 'unclassified']
    )) & (df.level == level)]

    if not microbial.empty:

        print(microbial)
        print(human)

        microbial_species = round(sum(microbial['percent']), 2)
        mreads = sum(microbial['reads'])

        overview_data.append(microbial_species)
        overview_labels.append(f'Microbial [{microbial_species}%]')
    else:
        microbial_species = 0.
        mreads = 0

    if len(overview_labels) == 2:
        # No human
        colors = ['#A9A9A9', target]
    else:
        # With human
        colors = ['#fec44f', '#A9A9A9', target]

    plot_annotated_donut(
        labels=overview_labels, data=overview_data, ax=ax, colors=colors
    )

    if human.empty:
        human_percent = 0
        human_reads = 0
    else:
        human_percent = human['percent']
        human_reads = human['reads']

    if unclassified.empty:
        unclassified_percent = 0
    else:
        unclassified_percent = unclassified['percent']

    return float(human_percent), int(human_reads), float(unclassified_percent),\
        int(unclassified.reads), float(microbial_species), int(mreads)


def plot_major_composition(df, ax, level, other=0., color=None) -> pandas.DataFrame:

    # Major

    df2 = df[(~df.taxonomy.isin(
        ['Homo sapiens', 'Homo']
    )) & (df.percent >= 10.)]

    tax = df2[
        df2['level'] == level
    ].sort_values(by='percent', ascending=False)

    data = tax.percent.tolist()
    labels = tax.taxonomy.tolist()

    # Minor

    df3 = df[(~df.taxonomy.isin(
        ['Homo sapiens', 'Homo']
    )) & (df.percent < 10.)]

    tax2 = df3[
        df3['level'] == level
    ].sort_values(by='percent', ascending=False)

    minor = round(sum(tax2.percent.values), ndigits=2)

    print(minor, tax2)

    data.insert(0, minor)
    labels.insert(0, 'Minor')

    labels = [f'{l}: [{data[i]}%]' for i, l in enumerate(labels)]

    color = sns.color_palette(color, len(data))
    color[0] = '#808080'

    plot_annotated_donut(
        labels=labels, data=data, ax=ax, colors=color
    )

    return tax


def plot_minor_composition(df, ax, level, top) -> pandas.DataFrame:

    # Remove human and filter percent < 10%

    df2 = df[(~df.taxonomy.isin(
        ['Homo sapiens', 'Homo']
    )) & (df.percent < 10.)]

    tax = df2[
        df2['level'] == level
    ].sort_values(by='percent', ascending=False)

    df = tax[:top]

    df = df.sort_values(by='percent', ascending=True)

    df.plot(
        y='percent', x='taxonomy', kind='barh', ax=ax, legend=False, color='#808080'
    )

    values = df.percent.tolist()
    ax.set_xlabel('Percent of reads classified (%)')
    ax.set_ylabel('')
    ax.set_xlim(0, 10)

    # set individual bar lables using above list
    for i, patch in enumerate(ax.patches):
        # get_width pulls left or right; get_y pushes up or down
        ax.text(
            patch.get_width() + .3, patch.get_y() + 0.1,
            str(values[i]) + "%",
            fontsize=10,
            color='dimgrey'
        )

    return tax


def plot_annotated_donut(labels: [str], data: [float], ax: None, colors=None):

    """ Donut chart with labels example
    """

    if len(labels) != len(data):
        raise ValueError('Data and label lists most be of the same length.')

    wedges, texts = ax.pie(
        data, wedgeprops=dict(width=0.5), startangle=-40, colors=colors
    )

    bbox_props = dict(
        boxstyle="square,pad=0.3", fc="w", ec="k", lw=0.72
    )

    kw = dict(
        arrowprops=dict(arrowstyle="-"),
        bbox=bbox_props, zorder=0, va="center")

    for i, p in enumerate(wedges):
        ang = (p.theta2 - p.theta1) / 2. + p.theta1
        y = np.sin(np.deg2rad(ang))
        x = np.cos(np.deg2rad(ang))
        horizontalalignment = {-1: "right", 1: "left"}[int(np.sign(x))]
        connectionstyle = "angle,angleA=0,angleB={}".format(ang)
        kw["arrowprops"].update({"connectionstyle": connectionstyle})
        ax.annotate(labels[i], xy=(x, y), xytext=(1.35 * np.sign(x), 1.4 * y),
                    horizontalalignment=horizontalalignment, **kw, fontsize=16)

