""" Modules to process program output in the pipelines """

import json
import pandas
import logging
import numpy as np
import seaborn as sns

from matplotlib import pyplot as plt
from pySankey import sankey

from pathlib import Path
from nanopath.utils import PoreLogger, create_fastx_index


class AssemblyProcessor(PoreLogger):

    """ Process results from metagenome assembly """

    def __init__(
          self,
          fasta: Path,
          assembler: str = 'flye',
          info_file: Path = None,
          verbose: bool = True
    ):

        PoreLogger.__init__(
            self, level=logging.INFO if verbose else logging.ERROR
        )

        self.fasta = fasta
        self.assembler = assembler

        self.fxi = None
        self.fxi_file = None

        self.info_file = info_file

        # Some basic filters to extract potential candidates from Flye

        self.min_bacteria = 500000  # minimum length to consider chromosome
        self.max_virus = 100000  # max length to consider linear virus
        self.min_virus_mult = 3  # minimum multiplicity to look for viruses

    def get_server_data(self) -> dict:

        return self.get_assembly_data()

    def get_assembly_data(self):

        """ What do we need to know about the assembly?

        Gets a bunch of data from the assembly output directory
        and processes it according to selected assembler - currently Flye.

        """

        if self.assembler == 'flye':
            self.logger.info('Read assembly information file from Flye')

            info = self.read_assembly_info()

            # Simple check for chromosomes / plasmids /viruses

            bacterial_chromosomes = info.loc[
                (info['length'] > self.min_bacteria) &
                (info['circular'] == '+'), :
            ]

            # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4379921/

            plasmids = info.loc[
                (info['length'] < self.min_bacteria) &
                (info['circular'] == '+'), :
            ]

            viruses = info.loc[
                (info['length'] < self.max_virus) &
                (info['circular'] == '-') &
                (info['multiplicity'] > self.min_virus_mult), :
            ]

            return {
                'contigs': self.data_to_report(data=info),
                'candidates': {
                    'bacteria': self.data_to_report(
                        data=bacterial_chromosomes
                    ),
                    'viruses': self.data_to_report(data=viruses),
                    'plasmids': self.data_to_report(data=plasmids)
                }
            }

    @staticmethod
    def data_to_report(data: pandas.DataFrame):

        """ Transform contig data to JSON for Vue application """

        return data.to_dict(orient='records')

    def read_fasta(self):

        """ Read the Fasta file into a FastxIndex """

        self.fxi, self.fxi_file = create_fastx_index(fastx=self.fasta)

    def read_assembly_info(self) -> pandas.DataFrame:

        """ Read the Flye assembly information file """

        return pandas.read_csv(
            self.info_file, header=0, sep='\t',
            names=[
                'name',
                'length',
                'coverage',
                'circular',
                'repeat',
                'multiplicity'
            ], usecols=[0, 1, 2, 3, 4, 5]
        )


class KrakenProcessor(PoreLogger):

    """ Process results from Kraken2 """

    def __init__(
        self,
        report: Path,
        reads: Path,
        level: str = 'S',
        verbose: bool = True
    ):

        PoreLogger.__init__(
            self, level=logging.INFO if verbose else logging.ERROR
        )

        self.report = report
        self.reads = reads

        self.level = level
        self.host = 'Homo sapiens'  # for now always human

        self.report_data = self.read_report()
        self.read_data = self.read_classifications()

        self.total_reads = self.report_data.loc[
            self.report_data['taxonomy'].isin(
                ('unclassified', 'root')
            ), 'reads'
        ].sum()

    def get_server_data(self):

        (bacteria_report, virus_report, contamination_report), \
            (bacteria_data, virus_data) = self.process_microbial()

        # Get server data format for donut visualization
        bacteria_server = self.report_to_json(df=bacteria_report)
        virus_server = self.report_to_json(df=virus_report)
        contamination_server = self.report_to_json(df=contamination_report)

        # Domain summary plot:
        host_percent, host_reads, _ = self.get_host()
        unclassified_percent, unclassified_reads, _ = self.get_unclassified()

        microbial_reads = bacteria_data[1] + virus_data[1]
        microbial_percent = bacteria_data[0] + virus_data[0]

        summary_server = [
            {
                'species': 'Host',
                'reads': host_reads,
                'percent': host_percent
            },
            {
                'species': 'Unclassified',
                'reads': unclassified_reads,
                'percent': unclassified_percent
            },
            {
                'species': 'Microbial',
                'reads': microbial_reads,
                'percent': microbial_percent,
            }
        ]

        return {
            'bacteria': bacteria_server,
            'viruses': virus_server,
            'contamination': contamination_server,
            'summary': summary_server
        }

    def get_unclassified(self) -> (float, int, pandas.DataFrame) or dict:

        """ Get unclassified reads and compute summaries from report

        :returns
            None, if no unclassified reads detected

            (percent reads, total reads, reads), if server_json = False
                reads are a subset of the read classifications from Kraken

            {percent: float, total: int, species: 'Unclassfied'} if server_json = True

        """

        unclassified = self.report_data[
            self.report_data['taxonomy'] == 'unclassified'
        ]

        nonsense = self.report_data[
            self.report_data['taxonomy'].isin(
                ('root', 'cellular organisms')
            )
        ]

        if nonsense.empty:
            total_nonsense = 0
            percent_nonsense = 0.
        else:
            total_nonsense = int(
                nonsense['direct'].sum()
            )
            percent_nonsense = float(
                nonsense['direct'].sum() / self.total_reads
            )

        if unclassified.empty:
            total_unclassified = 0
            percent_unclassified = 0.
        else:
            total_unclassified = int(
                unclassified['reads']
            )
            percent_unclassified = float(
                unclassified['reads'] / self.total_reads
            ) + float(percent_nonsense)

        final_unclassified = total_unclassified + total_nonsense

        if final_unclassified == 0:
            return 0., 0, pandas.DataFrame()

        reads = pandas.concat((
            self.read_data[self.read_data['classified'] == 'U'],
            self.read_data[self.read_data['taxon'].isin(
                ('root (taxid 1)', 'cellular organisms (taxid 131567)')  # GTDB [#6]
            )]
        ))

        return percent_unclassified, total_unclassified, reads

    def get_host(self) -> (float, int, pandas.DataFrame) or dict:

        """ Get host reads and compute summaries from report

        :returns percent reads, total reads, reads
            where reads are a subset of the read classifications from Kraken

        """

        host = self.report_data[
            self.report_data['taxonomy'] == self.host
        ]

        if host.empty:
            return 0., 0, []
        else:
            percent = float(host['reads'] / self.total_reads)
            total = int(host['reads'])

        read_idx = []
        for i, row in self.read_data.iterrows():
            name, taxon = row['read'], row['taxon']
            if self.host in taxon:
                read_idx.append(i)
        reads = self.read_data.iloc[read_idx]

        return percent, total, reads  # float, int, pandas.DataFrame

    def get_microbial(self) -> (tuple, tuple, tuple):

        """ Get microbial data for microbial domains from report

        Curently excludes eukaryotic pathogens; included:
        Bacteria, Viruses, Archaea

        :returns
            None, if no unclassified reads detected

            (percent reads, total reads, reads), if server_json = False
                reads are a subset of the read classifications from Kraken

            {percent: float, total: int, species: 'Host'} if server_json = True

        """

        bacteria, viruses, archaea = \
            self.get_subdomain(self.report_data, 'Bacteria'), \
            self.get_subdomain(self.report_data, 'Viruses'), \
            self.get_subdomain(self.report_data, 'Archaea')

        if bacteria is not None:
            bacteria_data = \
                float(bacteria.loc[0, 'reads'] / self.total_reads), \
                int(bacteria.loc[0, 'reads']), \
                bacteria
        else:
            self.logger.info(
                'Did not find any bacterial classifications in report.'
            )
            bacteria_data = (0., 0, pandas.DataFrame())

        if viruses is not None:
            viruses_data = \
                float(viruses.loc[0, 'reads'] / self.total_reads), \
                int(viruses.loc[0, 'reads']), \
                viruses
        else:
            self.logger.info(
                'Did not find any viral classifications in report.'
            )
            viruses_data = (0., 0, pandas.DataFrame())

        if archaea is not None:
            archaea_data = \
                float(archaea.loc[0, 'reads'] / self.total_reads), \
                int(archaea.loc[0, 'reads']), \
                archaea
        else:
            self.logger.info(
                'Did not find any archaeal classifications in report.'
            )
            archaea_data = (0., 0, pandas.DataFrame())

        return bacteria_data, viruses_data, archaea_data

    def process_microbial(self) -> (tuple, tuple) or None:

        """ Summarise Kraken report data for server output """

        # Get microbial
        bacteria_data, virus_data, archaea_data = self.get_microbial()

        # TODO: Identify contamination, needs work
        decontaminated, contamination = self.find_contamination(
            dict(
                Bacteria=bacteria_data[-1],
                Viruses=virus_data[-1],
                Archaea=archaea_data[-1]
            )
        )

        # Species classifications only!

        try:
            bacteria_report = decontaminated[0]
            bacteria_report = bacteria_report.loc[
              bacteria_report['level'] == self.level, :
            ]
        except IndexError:
            bacteria_report = pandas.DataFrame()

        try:
            virus_report = decontaminated[1]
            virus_report = virus_report.loc[
               virus_report['level'] == self.level, :
            ]
        except IndexError:
            virus_report = pandas.DataFrame()

        contamination_report = contamination.loc[
           contamination['level'] == self.level, :
        ]

        return (bacteria_report, virus_report, contamination_report), \
               (bacteria_data, virus_data)

    # Support methods

    def get_reads_from_report(self, report: pandas.DataFrame) -> pandas.DataFrame:

        read_idx = []
        allowed_taxa = report.taxonomy.tolist()
        for i, row in self.read_data.iterrows():
            name = row['taxon'].split('(')[0].strip()
            if name in allowed_taxa:
                read_idx.append(i)
        reads = self.read_data.iloc[read_idx].copy()

        reads.index.name = report.index.name

        return reads

    @staticmethod
    def find_contamination(data: dict) -> (tuple, pandas.DataFrame):

        """ Identify potential contaminants in Kraken report domains """

        decontaminated = []
        contaminated = []
        for name, df in data.items():
            if not df.empty:
                if name == 'Archaea':
                    # archaeal reads
                    contaminated.append(df)
                else:
                    # singleton reads
                    decon = df.loc[df['reads'] > 1, :]
                    conta = df.loc[df['reads'] <= 1, :]
                    decontaminated.append(decon)
                    contaminated.append(conta)

        return decontaminated, pandas.concat(contaminated)

    def report_to_json(self, df: pandas.DataFrame) -> list:

        if df.empty:
            return []
        else:
            return [
                {
                    'species': row['taxonomy'],
                    'percent': float(row['reads'] / self.total_reads),
                    'reads': row['reads']
                }
                for i, row in df.iterrows()
            ]

    def check_classifications(self, host_prefix='human'):

        # Get unclassified, check how many host, how many microbial

        upercent, utotal, ureads = self.get_unclassified()
        uhost, umicrobial = self._assess(ureads, host_prefix)

        self.logger.info(f'Unclassified reads: {uhost} host, {umicrobial} microbial')

        # Get reads classified as host, check how many host, microbial

        hpercent, htotal, hreads = self.get_host()
        hhost, hother = self._assess(hreads, host_prefix)

        htotal = hhost+hother

        self.logger.info(
            f'Reads classified as host: {hhost}/{htotal} true positives (host)'
        )
        self.logger.info(
            f'Reads classified as host: {hother}/{htotal} false positives (other)'
        )

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

            predicted_domain = self.get_superdomain(
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

        negative_data.to_csv(
            f'{prefix}.negatives.tsv', sep='\t', index=False
        )
        summary_proportions.to_csv(
            f'{prefix}.props.tsv', sep='\t', index=True
        )
        summary.to_csv(
            f'{prefix}.reads.tsv', sep='\t', index=True
        )

    def get_superdomain(self, taxon_name: str):

        """ Get a domain-specific subset of the report dataframe"""

        df = self.report_data.reset_index()

        try:
            taxon_index = [
                i for i, taxon in enumerate(
                    df['taxonomy']
                ) if taxon in taxon_name
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

    def get_subdomain(
        self, report: pandas.DataFrame, taxon: str
    ) -> pandas.DataFrame or None:

        """ Get a domain-specific subset of the report dataframe"""

        df = report.reset_index()

        try:
            taxon_index = int(df[df['taxonomy'] == taxon].index.values)
        except TypeError:
            self.logger.warning(f'Could not get subdomain of {taxon}')
            return None

        domain_index = 0
        for i, level in enumerate(
                df['level'][taxon_index:]
        ):
            if level == 'D' and i > 0:
                domain_index = taxon_index + i
                break

        if domain_index == 0:
            return report.iloc[taxon_index:, :].reset_index()
        else:
            return report.iloc[taxon_index:domain_index, :].reset_index()

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

        categories = data[true].unique().tolist() + \
                     data[predicted].unique().tolist()

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
            boxstyle="square,pad=0.3", fc="w", ec="k", lw=0
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
