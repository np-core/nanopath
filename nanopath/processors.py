""" Modules to process program output from the pipelines """

import json
import pandas
import logging
import numpy as np
import seaborn as sns

from matplotlib import pyplot as plt
from pySankey import sankey

from pathlib import Path
from nanopath.utils import PoreLogger


class KrakenProcessor(PoreLogger):

    def __init__(
        self,
        report: Path,
        reads: Path,
        level: str = 'S',
        host: str = 'Homo sapiens',
        verbose: bool = True
    ):

        PoreLogger.__init__(
            self, level=logging.INFO if verbose else logging.ERROR
        )

        self.report = report
        self.reads = reads

        self.level = level
        self.host = host

        self.report_data = self.read_report()
        self.read_data = self.read_classifications()

    def get_unclassified(self):

        unclassified = self.report_data[
            self.report_data['taxonomy'] == 'unclassified'
        ]

        percent = float(unclassified['percent'])
        total = int(unclassified['reads'])
        reads = self.read_data[self.read_data['classified'] == 'U']

        return reads, percent, total

    def get_host(self):

        host = self.report_data[
            self.report_data['taxonomy'] == self.host
        ]

        percent = float(host['percent'])
        total = int(host['reads'])

        read_idx = []
        for i, row in self.read_data.iterrows():
            name, taxon = row['read'], row['taxon']
            if self.host in taxon:
                read_idx.append(i)
        reads = self.read_data.iloc[read_idx]

        return reads, percent, total

    def check_classifications(self, host_prefix='human'):

        # Get unclassified, check how many host, how many microbial

        ureads, upercent, utotal = self.get_unclassified()
        uhost, umicrobial = self._assess(ureads, host_prefix)

        self.logger.info(f'Unclassified reads: {uhost} host, {umicrobial} microbial')

        # Get reads classified as host, check how many host, microbial

        hreads, hpercent, htotal = self.get_host()
        hhost, hother = self._assess(hreads, host_prefix)

        htotal = hhost+hother

        self.logger.info(
            f'Reads classified as host: {hhost}/{htotal} true positives (host)'
        )
        self.logger.info(
            f'Reads classified as host: {hother}/{htotal} false positives (other)'
        )

    def get_microbial(self):

        pass

    def assess_composition(self, compose_file: Path, host_prefix='human'):

        """ Assess a composed mixture for host / microbial classification

         Given a composition file matching the reads from the input
         classification, investigate how each species in the mixture
         is classified.

        :param compose_file: composition file JSON format from `np utils compose`
        :param host_prefix: the prefix usedi n the composition file as host reads

        """

        # Read the compose file for the report and

        with compose_file.open() as jin:
            composition = json.load(jin)

        if host_prefix not in composition.keys():
            raise ValueError(f'Host prefix {host_prefix} not in compose file.')

        return self._assess_composition(
            composition=composition, host_prefix=host_prefix
        )

    def _assess_composition(self, composition: dict, host_prefix: str):

        true_positives, negatives, unclassified, nonsense = {}, {}, {}, {}
        negative_data = []

        for i, row in self.read_data.iterrows():
            name, taxon = row['read'], row['taxon']

            try:
                prefix = name.split('_')[0]
            except IndexError:
                self.logger.error(
                    'Something is wrong with the compose file prefixes.'
                )
                raise

            try:
                prefix_data = composition[prefix]
            except KeyError:
                self.logger.error(
                    f'Could not find {prefix} in composition.'
                )
                raise

            taxon_name = prefix_data['taxon_name']

            if taxon_name in taxon:
                if prefix not in true_positives.keys():
                    true_positives[prefix] = 1
                else:
                    true_positives[prefix] += 1
            elif 'unclassified' in taxon:
                if prefix not in unclassified.keys():
                    unclassified[prefix] = 1
                else:
                    unclassified[prefix] += 1
            elif 'root' in taxon or 'cellular organism' in taxon:
                if prefix not in nonsense.keys():
                    nonsense[prefix] = 1
                else:
                    nonsense[prefix] += 1
            else:
                if prefix not in negatives.keys():
                    negatives[prefix] = 1
                else:
                    negatives[prefix] += 1

                negative_data.append((name, taxon))

        totals = {}  # Fill in null counts and get total count
        for prefix in composition.keys():
            for count in (true_positives, negatives, unclassified, nonsense):
                if prefix not in count.keys():
                    count[prefix] = 0

                if prefix not in totals.keys():
                    totals[prefix] = count[prefix]
                else:
                    totals[prefix] += count[prefix]

        negative_data = pandas.DataFrame(negative_data)

        summary_data = (true_positives, negatives, unclassified, nonsense, totals)

        summary = pandas.DataFrame(
            summary_data,
            index=['positives', 'negatives', 'unclassified', 'nonsense', 'totals']
        )

        summary_proportions = (summary / summary.loc['totals', :])
        summary_proportions = summary_proportions.drop('totals', axis=0)

        # Assess negative data

        microbial, conflict, domain = [], [], []

        negative_data.columns = ['prefix', 'predicted']
        for i, row in negative_data.iterrows():
            prefix = row['prefix'].split('_')[0]

            if prefix == host_prefix:
                true_domain = 'Host'
            else:
                true_domain = 'Microbial'  # Bacteria, Archaea, Viruses !

            predicted_domain = self._get_super_domain(
                taxon_name=row['predicted']
            )

            microbial.append(true_domain)
            domain.append(predicted_domain)

            if predicted_domain in ('Bacteria', 'Archaea', 'Viruses') \
                    and true_domain == 'Microbial':
                conflict.append('no')
            elif predicted_domain == 'Viruses' and 'Human gammaherpesvirus 4' in row['predicted']:
                conflict.append('retro')
            elif predicted_domain == 'Viruses' and true_domain == 'Host':
                conflict.append('yes')
            elif predicted_domain in ('Bacteria', 'Archaea', 'Viruses') \
                    and true_domain == 'Host':
                conflict.append('yes')
            elif predicted_domain == 'Host' and true_domain == 'Microbial':
                conflict.append('yes')
            else:
                conflict.append('unknown')

        negative_data['predicted_domain'] = domain
        negative_data['prefix_class'] = microbial
        negative_data['conflict'] = conflict

        negative_data.index.name = 'negative'
        negative_data.sort_values(
            ['conflict'], inplace=True, ascending=False
        )

        return summary, summary_proportions, negative_data

    def create_summary_plot(
        self,
        summary,
        summary_proportions,
        negative_data,
        prefix: str = 'summary',
        palette: str = 'Blues'
    ):

        fig, axes = plt.subplots(
            nrows=3, ncols=1, figsize=(1 * 12, 3 * 9)
        )

        fig.subplots_adjust(hspace=0.8)

        if len(axes) != len(summary_proportions.columns):
            raise ValueError('Axes must be of same length as columns.')

        summary_proportions.loc['negatives', :] = \
            summary_proportions.loc['nonsense', :] + \
            summary_proportions.loc['negatives', :]

        summary_proportions.drop('nonsense', axis=0, inplace=True)

        for i, column in enumerate(summary_proportions):
            c = summary_proportions[column]
            labels = [
                name + f" [{round(c[i]*100, 2)}%]"
                for i, name in enumerate(c.index)
            ]
            self.plot_annotated_donut(
                labels=labels, data=c.tolist(),
                palette=palette, ax=axes[i], title=column
            )

        plt.tight_layout()
        fig.savefig(f'{prefix}.png')

        self.plot_sankey(data=negative_data, prefix=prefix)

        negative_data.to_csv(f'{prefix}.negatives.tsv', sep='\t', index=False)
        summary_proportions.to_csv(f'{prefix}.props.tsv', sep='\t', index=False)
        summary.to_csv(f'{prefix}.reads.tsv', sep='\t', index=False)

    def _get_super_domain(self, taxon_name: str):

        """ Get a domain-specific subset of the report dataframe"""

        df = self.report_data.reset_index()

        try:
            taxon_index = [
                i for i, taxon in enumerate(df['taxonomy']) if taxon in taxon_name
            ][0]
        except IndexError:
            self.logger.error('Something went wrong in _get_super_domain')
            raise

        super_levels = df['level'][0:taxon_index].tolist()
        super_levels.reverse()

        domain_index = 0
        for i, level in enumerate(super_levels):
            if level == 'D' and i > 0:
                domain_index = taxon_index - i - 1
                break

        domain_row = df.iloc[domain_index, :]

        return str(domain_row['taxonomy'])

    @staticmethod
    def _assess(reads: pandas.DataFrame, host_prefix: str):

        host = 0
        microbial = 0
        for name in reads['read']:
            if name.startswith(host_prefix):
                host += 1
            else:
                microbial += 1

        return host, microbial

    def read_report(self):

        df = pandas.read_csv(
            self.report, header=None, sep="\t", names=[
                "percent", "reads", "direct", "level", "taxid", "taxonomy",
            ]
        )
        df.taxonomy = df.taxonomy.str.strip()  # strip indentation

        return df

    def read_classifications(self):
        return pandas.read_csv(
            self.reads, header=None, sep="\t", names=[
                'classified', 'read', 'taxon', 'node', 'path'
            ]
        )

    @staticmethod
    def plot_sankey(
        data: pandas.DataFrame,
        true: str = 'prefix_class',
        predicted: str = 'predicted_domain',
        palette: str = 'Blues',
        prefix: str = None,
    ):

        # if prefix is not None:
        #     data['prefix'] = [p.split('_')[0] for p in data['prefix']]
        #     data = data[data['prefix'] == prefix]

        categories = data[true].unique().tolist() + data[predicted].unique().tolist()

        colors = {
            categories[i]: c for i, c in enumerate(
                sns.color_palette(palette, len(categories))
            )
        }

        sankey.sankey(
            data[true],
            data[predicted],
            aspect=40,
            colorDict=colors,
            fontsize=12,
            figure_name=prefix if prefix is not None else 'negatives'
        )

    @staticmethod
    def plot_annotated_donut(
        labels: [str], data: [float], ax: None, palette="Blues", title: str = 'title'
    ):

        """ Donut chart with labels """

        if len(labels) != len(data):
            raise ValueError('Data and label lists most be of the same length.')

        color = sns.color_palette(palette, len(data))

        ax.set_title(title)

        wedges, texts = ax.pie(
            data, wedgeprops=dict(width=0.5), startangle=-40, colors=color
        )

        bbox_props = dict(
            boxstyle="square,pad=0.3", fc="w", ec="k", lw=0.72
        )

        kw = dict(
            arrowprops=dict(arrowstyle="-"),
            bbox=bbox_props, zorder=0, va="center"
        )

        for i, p in enumerate(wedges):

            ang = (p.theta2 - p.theta1) / 2. + p.theta1
            y = np.sin(np.deg2rad(ang))
            x = np.cos(np.deg2rad(ang))

            horizontalalignment = {-1: "right", 1: "left"}[int(np.sign(x))]
            connectionstyle = "angle,angleA=0,angleB={}".format(ang)
            kw["arrowprops"].update({"connectionstyle": connectionstyle})

            ax.annotate(
                labels[i], xy=(x, y), **kw,
                xytext=(1.35 * np.sign(x), 1.4 * y),
                horizontalalignment=horizontalalignment,
                fontsize=16
            )
