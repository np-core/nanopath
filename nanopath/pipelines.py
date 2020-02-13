""" Module to supplement and process server-side analysis pipelines """

import json
import pandas
from datetime import datetime
from pathlib import Path
import seaborn as sns
import numpy as np
from matplotlib import pyplot as plt


class Pipeline:

    def __init__(
        self,
        uuid: str = None,
        user: str = None,
        name: str = None,
        submitted: datetime = None,
        host_removal: str = 'minimap2'
    ):

        # UI client side
        self.uuid = uuid
        self.user = user
        self.name = name
        self.submitted = submitted
        self.host_removal = host_removal

    def set_meta(
        self,
        uuid: str = None,
        user: str = None,
        name: str = None,
        submitted: datetime = None,
        host_removal: str = 'minimap2'
    ):

        """ Meta data from request submission of pipeline run """

        self.uuid = uuid
        self.user = user
        self.name = name
        self.submitted = submitted
        self.host_removal = host_removal


class Metagenome(Pipeline):

    """ Pipeline processing for NP-SEPSIS currently excludes Eukaryots / PARASITES / FUNGI """

    def __init__(self):

        Pipeline.__init__(self)

        self.pipeline = 'sepsis'

        self.kraken = {}
        self.level = "S"
        self.host = "Homo sapiens"

        self.overview = None
        self.microbial = None

        # Pipeline server telemetry
        self.completed = datetime.now().strftime("%d-%m-%Y-%H:%M:%S")
        self.database = 'minikraken2'

    def get_server_data(self, fout: Path = None):

        server_data = {
            'meta': dict(
                uuid=self.uuid,
                submitted=self.submitted,
                completed=self.completed,
                host=self.host,
            ),
            'summary': self.kraken['summary'],
            'summary_plot': [
                {
                    'percent': float(
                        self.kraken['summary']['unclassified_percent']
                    ),
                    'reads': int(
                        self.kraken['summary']['unclassified_reads']
                    ),
                    'species': "Unclassified"
                },
                {
                    'percent': float(
                        self.kraken['summary']['microbial_percent']
                    ),
                    'reads': int(
                        self.kraken['summary']['microbial_reads']
                    ),
                    'species': "Microbial"
                },
                {
                    'percent': float(
                        self.kraken['summary']['host_percent']
                    ),
                    'reads': int(
                        self.kraken['summary']['host_reads']
                    ),
                    'species': "Host"
                },
            ],
            'major_species_plot': [
                 {
                     'percent': float(
                         self.kraken['summary']['minor_species_percent']
                     ),
                     'reads': int(
                         self.kraken['summary']['minor_species_reads']
                     ),
                     'species': "Minor"
                 }
            ] + [
                {'percent': float(row[0]), 'reads': int(row[1]), 'species': row[-1]} for row in [
                    row.tolist() for _, row in self.kraken['major_species'].iterrows()
                ]
             ],
            'minor_species_plot': [
                {'percent': float(row[0]), 'reads': int(row[1]), 'species': row[-1]} for row in [
                    row.tolist() for _, row in self.kraken['minor_species'].iterrows()
                ]
            ],
        }

        if fout:
            with Path(f'{fout}').open('w') as json_out:
                json.dump(server_data, json_out)

    def process_kraken(self, report: Path):

        df = pandas.read_csv(
            report, header=None, sep="\t", names=[
                "percent", "reads", "direct", "level", "taxid", "taxonomy",
            ]
        )
        df.taxonomy = df.taxonomy.str.strip()

        report_summary = self.get_report_data(df)

        minor_species = self.microbial[
            (self.microbial.percent < 0.1) & (self.microbial.level == self.level)
        ].sort_values('percent', ascending=False)

        major_species = self.microbial[
            (self.microbial.percent >= 0.1) & (self.microbial.level == self.level)
        ].sort_values('percent', ascending=False)

        report_summary['minor_species_percent'] = float(
            sum(minor_species.percent.values)
        )
        report_summary['minor_species_reads'] = int(
            sum(minor_species.reads.values)
        )

        report_summary['major_species_percent'] = float(
            sum(major_species.percent.values)
        )
        report_summary['major_species_reads'] = int(
            sum(major_species.reads.values)
        )

        self.kraken = {
            'summary': report_summary,
            'minor_species': minor_species,
            'major_species': major_species
        }

        return self.kraken

    def plot_overview(self, overview_data, overview_labels, palette: str, ax=None):

        color = sns.color_palette(palette, 2)[-1]

        if len(overview_labels) == 2:
            # No human
            colors = ['#A9A9A9', color]
        else:
            # With human
            colors = ['#fec44f', '#A9A9A9', color]

        self.plot_annotated_donut(
            labels=overview_labels, data=overview_data, ax=ax, colors=colors
        )

    def get_report_data(self, df: pandas.DataFrame) -> (dict, tuple):

        df.percent = df.percent/100

        # Overview data for plotting:

        data, labels = [], []

        classified_percent = float(
            df.loc[df['taxonomy'] == "root", "percent"]
        )
        classified_reads = int(
            df.loc[df['taxonomy'] == "root", "reads"]
        )

        # Human classification

        host = df[df['taxonomy'] == self.host]
        if not host.empty:
            host_percent = float(host['percent'])
            host_reads = int(host['reads'])

            data.append(host_percent)
            labels.append(f'Host [{round(host_percent*100, 4)}%]')
        else:
            host_percent = 0.
            host_reads = 0

        # Unclassified

        unclassified = df[df['taxonomy'] == 'unclassified']
        if not unclassified.empty:
            try:
                unclassified_percent = float(unclassified['percent'])
                unclassified_reads = int(unclassified['reads'])
            except TypeError:
                raise

            data.append(unclassified_percent)
            labels.append(f'Unclassified [{round(unclassified_percent*100, 4)}%]')
        else:
            unclassified_percent = 0.
            unclassified_reads = 0

        # Taxon level microbial classification - EXCLUDES METAZOAN PARASITES / FUNGI

        bacteria = self.get_subtaxa(df, 'Bacteria')

        try:
            bacteria_percent, bacteria_reads = \
                bacteria.percent.values[0], bacteria.reads.values[0]
        except IndexError:
            bacteria_percent, bacteria_reads = 0., 0
        viruses = self.get_subtaxa(df, 'Viruses')

        try:
            viruses_percent, viruses_reads = \
                viruses.percent.values[0], viruses.reads.values[0]
        except IndexError:
            viruses_percent, viruses_reads = 0., 0

        archaea = self.get_subtaxa(df, 'Archaea')

        try:
            archaea_percent, archaea_reads = \
                archaea.percent.values[0], archaea.reads.values[0]
        except IndexError:
            archaea_percent, archaea_reads = 0., 0

        microbial_percent = bacteria_percent + viruses_percent + archaea_percent
        microbial_reads = bacteria_reads + viruses_reads + archaea_reads

        data.append(microbial_percent)
        labels.append(f'Microbial [{round(microbial_percent*100, 4)}%]')

        self.overview = (data, labels)
        self.microbial = pandas.concat((bacteria, viruses, archaea))

        return {
            "microbial_percent": float(microbial_percent),
            "microbial_reads": int(microbial_reads),
            "host_percent": float(host_percent),
            "host_reads": int(host_reads),
            "unclassified_percent": float(unclassified_percent),
            "unclassified_reads": int(unclassified_reads),
            "bacteria_percent": float(bacteria_percent),
            "bacteria_reads": int(bacteria_reads),
            "viruses_percent": float(viruses_percent),
            "viruses_reads": int(viruses_reads),
            "archaea_percent": float(archaea_percent),
            "archaea_reads": int(archaea_reads),
            "classified_percent": float(classified_percent),
            "classified_reads": int(classified_reads),
            "total_reads": int(classified_reads+unclassified_reads)
        }

    def plot_kraken_summary(
        self,
        plot_file: Path = Path('kraken.png'),
        palette: str = "Greens",
        top_minor: int = 10,
        title: str = "",
        subtitles: bool = False
    ):

        """ Matplotlib summary plot """

        fig, (ax1, ax2, ax3) = plt.subplots(
            nrows=1, ncols=3, figsize=(27.0, 9),
            subplot_kw=dict(aspect="equal")
        )

        self.plot_overview(
            self.overview[0], self.overview[1], palette=palette, ax=ax1
        )

        self.plot_major_species(palette=palette, ax=ax2)
        self.plot_minor_species(top=top_minor, ax=ax3)

        total = self.kraken['summary']["total_reads"]
        if subtitles:
            ax1.title.set_text(f'Total reads ({total})')
            ax2.title.set_text(f'Major Taxa (> 10%)')
            ax3.title.set_text(f'Minor Taxa (<= 10%)')

        plt.tight_layout()
        fig.subplots_adjust(top=1)

        fig.suptitle(title, fontsize=24)

        fig.savefig(
            f'{plot_file}', pad_inches=0.5
        )

    @staticmethod
    def get_subtaxa(
        report: pandas.DataFrame, taxon: str
    ) -> pandas.DataFrame:

        """ Get a domain-specific subset of the report dataframe"""

        df = report.reset_index()

        taxon_index = int(df[df['taxonomy'] == taxon].index.values)

        domain_index = 0
        for i, level in enumerate(
            df['level'][taxon_index:]
        ):
            if level == 'D' and i > 0:
                domain_index = taxon_index+i
                break

        if domain_index == 0:
            return report.iloc[taxon_index:, :]
        else:
            return report.iloc[taxon_index:domain_index, :]

    @staticmethod
    def plot_annotated_donut(
            labels: [str], data: [float], ax: None, colors=None
    ):

        """ Donut chart with labels """

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

    def plot_major_species(self, palette: str = 'Greens', ax=None):

        data = self.kraken['major_species'].percent.tolist()
        labels = self.kraken['major_species'].taxonomy.tolist()

        data = [round(p*100, 2) for p in data]

        minor = self.kraken['summary']['minor_species_percent']
        data.insert(0, round(minor*100, 2))
        labels.insert(0, 'Minor')

        color = sns.color_palette(palette, len(data))
        color[0] = '#808080'

        labels = [f'{l}: [{data[i]}%]' for i, l in enumerate(labels)]
        self.plot_annotated_donut(
            labels=labels, data=data, ax=ax, colors=color
        )

    def plot_minor_species(self, top: int = 10, ax=None):

        df = self.kraken['minor_species'][:top]

        df = df.sort_values(by='percent', ascending=True)

        df['data'] = [round(p*100, 2) for p in df['percent']]

        df.plot(
            y='data',
            x='taxonomy',
            kind='barh',
            ax=ax,
            legend=False,
            color='#808080'
        )

        values = df.data.tolist()
        ax.set_xlabel('Percent of reads classified (%)')
        ax.set_ylabel('')
        ax.set_xlim(0, 10)

        # set individual bar labels using above list
        for i, patch in enumerate(ax.patches):
            # get_width pulls left or right; get_y pushes up or down
            ax.text(
                patch.get_width() + .3, patch.get_y() + 0.1,
                str(values[i]) + "%",
                fontsize=10,
                color='dimgrey'
            )