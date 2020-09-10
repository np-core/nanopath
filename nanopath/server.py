from jinja2 import Template, Environment, FileSystemLoader
from weasyprint import HTML

from nanopath.pipelines import PathogenPipeline
from nanopath.utils import PoreLogger
import datetime
from pathlib import Path
from getpass import getuser
from uuid import uuid4

import base64
import logging
import yaml


class AppServer(PoreLogger):

    def __init__(self, pathogen_path: Path):

        PoreLogger.__init__(self, level=logging.INFO, name="AppServer")

        self.pipelines = {}

        self.pathogen_path = pathogen_path

        if self.pathogen_path:
            self.pipelines['pathogen'] = PathogenPipeline()

        pass

    # Individual pipeline collector functions

    def clean(self):

        for pipeline in self.pipelines.values():
            pipeline.clean()  # remove working dirs

    def collect_pathogen_results(
        self,
        groupby: str = r'barcode\d\d',
        host_species: str = "Homo sapiens",
        major_threshold: float or int = 0.1,
        minor_threshold: float = 0.5,
        level: str = "S"
    ):

        # Group files by regular expression

        db_reports, db_reads = self.pipelines['pathogen'].collect_results(
            path=self.pathogen_path, groupby=groupby, level=level
        )

        db_server = self.pipelines['pathogen'].get_server_data(
            host_label=host_species,
            major_threshold=major_threshold,
            minor_threshold=minor_threshold
        )

        # Total reads by database must all be the same (same input)
        # TODO: Add exception here if not all the same
        group_reads = {
            group: total for db, groups in db_reads.items()
            for group, total in groups.items()
        }

        return db_reports, group_reads, db_server

    @staticmethod
    def collect_pathogen_report(
        app_data: dict, db_reports: dict, group_reads: dict, outdir: Path
    ):

        report = PathogenReport(
            app_data=app_data,
            db_reports=db_reports,
            group_reads=group_reads
        )

        report.process_report(remove_host=False, remove_controls=False)
        report.write_report(fname=outdir / f"{report.sample_id}_{report.uid}")


class PathogenReport(PoreLogger):

    def __init__(
        self,
        app_data: dict,
        db_reports: dict,
        group_reads: dict
    ):

        PoreLogger.__init__(self, level=logging.INFO, name="PathogenReport")

        self.app_data = app_data  # report configuration sent by app
        self.db_reports = db_reports  # database report tables
        self.group_reads = group_reads  # total reads by group / sample

        self.config = dict()
        self.report_data = None

        self.template_path = Path(__file__).parent / 'templates' / 'reports'

        self.sample_id = self.app_data['header']['sample_id']

        self.uuid = str(uuid4())
        self.uid = self.uuid.split('-')[0]

    def _init_report_data(self):

        sample = self.app_data['config']['sample']
        negative_control = self.app_data['config']['negative_control']

        # Things that should not go wrong:

        try:
            threshold = float(
                self.app_data['config']['threshold']
            )
        except TypeError:
            self.logger.info(
                f'Could not parse threshold value from application'
            )
            raise

        try:
            total_reads = self.group_reads[sample]['total']
        except KeyError:
            self.logger.info(
                f'Could not find {sample} in {self.group_reads} for report'
            )
            raise

        # Threshold settings

        if threshold >= 1:
            _threshold_column = 'reads'
            threshold_reads = int(threshold)
            threshold_percent = round(threshold / total_reads, 6)
        else:
            _threshold_column = 'percent'
            threshold_reads = int(threshold*total_reads)
            threshold_percent = round(threshold, 6)

        # Report data object

        report_id = f"{self.uid} ({getuser()})"

        report_data = dict(
            data=list(),
            header=self.app_data['header'],
            footer={
                'pipeline':  self.app_data['config']['pipeline'],
                'version': self.app_data['config']['version'],
                'threshold_reads': threshold_reads,
                'threshold_percent': threshold_percent,
                'total_reads': f"{total_reads:,}",  # analysed / processed
                'sample': sample,
                'negative_control': negative_control,
                'report_id': report_id
            }
        )

        report_data['header']['run_complete'] = \
            self.app_data['config']['run_complete']
        report_data['header']['date'] = \
            datetime.datetime.today().strftime("%d/%m/%Y")

        return report_data, threshold, _threshold_column

    def process_report(
        self, remove_host: bool = True, remove_controls: bool = False
    ):

        sample = self.app_data['config']['sample']
        host_species = self.app_data['config']['host_species']
        negative_control = self.app_data['config']['negative_control']

        report_data, threshold, _threshold_column = self._init_report_data()

        # Process report data, extract relevant app settings first
        if isinstance(negative_control, str):
            negative_control = [negative_control]

        for db, reports in self.db_reports.items():
            self.logger.info(f"DB {db}")
            print(reports)
            db_sample_data = reports[sample]
            total_reads = self.group_reads[sample]['total']

            if remove_host:
                db_sample_data = db_sample_data[
                    db_sample_data['tax'] == host_species
                ]

            # Negative control filter block
            db_negative_control_taxa = []
            if negative_control is not None:
                for ng in negative_control:
                    db_negative_control_taxa += reports[ng].tax.tolist()
                if remove_controls:
                    db_sample_data = db_sample_data[
                        ~db_sample_data['tax'].isin(db_negative_control_taxa)
                    ]

            # Threshold filter block
            db_sample_data = db_sample_data[
                db_sample_data[_threshold_column] >= threshold
            ]

            # Data summary block
            db_sample_data = db_sample_data.sort_values(
                'reads', ascending=False
            )
            labels = db_sample_data.tax.tolist()
            reads = db_sample_data.reads.tolist()

            percent = [
                round((x/total_reads)*100, 2) for x in reads
            ]  # bracken classified over kraken total

            # YAML report configuration of Netflow

            try:
                database_name = self.config['database_labels'][db]
            except KeyError:
                self.logger.debug(
                    f'Could not detect database label entry for: {db}'
                )
                database_name = db

            try:
                species_flags = self.config['species_flags']
            except KeyError:
                self.logger.debug(
                    f'Could not detect species flags in report configuration'
                )
                species_flags = dict()

            # Flag block:
            flags = []
            for label in labels:
                try:
                    label_flag = species_flags[label]  # string
                except KeyError:
                    self.logger.debug(
                        f'Could not find {label} in species flag configuration'
                    )
                    label_flag = ""

                if label in db_negative_control_taxa:
                    label_flag = "IN NEGATIVE CONTROL"

                flags.append(label_flag)

            classified = sum(reads)  # bracken at species level

            species_data = [{
                'label': label,
                'flags':  flags[i],
                'reads': f"{reads[i]:,}",
                'percent': percent[i],
                } for i, label in enumerate(labels)
            ]
            report_data['data'].append({
                'db': database_name,
                'classified_reads': f"{classified:,}",
                'classified_percent': round((classified/total_reads)*100, 2),
                'species': species_data
            })

        self.report_data = report_data

        return report_data

    def read_config(self, file: Path):

        with file.open() as fin:
            config = yaml.load(fin, Loader=yaml.BaseLoader)

        try:
            self.config = config['reports']['pathogen']
        except KeyError:
            self.logger.debug(
                f'Could not detect entry `pathogen` in report configuration'
            )
            self.config = dict()

    def write_report(self, fname: Path):

        loader = FileSystemLoader(
            str(self.template_path).replace('\\', '/')
        )

        # Images

        with (self.template_path / 'qg.jpg').open('rb') as image_file:
            bstring = b"data:image/jpg;base64,"+base64.b64encode(
                image_file.read()
            )

        env = Environment(loader=loader)
        template = env.get_template('pathogen.html')
        rendered = template.render(
            header__logo=bstring.decode('utf-8'),
            header__date=self.report_data['header']['date'],
            header__patient_name=self.report_data['header']['patient_name'],
            header__dob=self.report_data['header']['dob'],
            header__location=self.report_data['header']['location'],
            header__sample_type=self.report_data['header']['sample_type'],
            header__completed=self.report_data['header']['run_complete'],
            header__sample_id=self.report_data['header']['sample_id'],
            header__sequenced_from=self.report_data['header']['sequenced_from'],
            header__requested_by=self.report_data['header']['requested_by'],
            header__contact=self.report_data['header']['contact'],
            data__databases=self.report_data['data'],
            footer__pipeline=self.report_data['footer']['pipeline'],
            footer__version=self.report_data['footer']['version'],
            footer__barcode=self.report_data['footer']['sample'],
            footer__report_id=self.report_data['footer']['report_id'],
            footer__total_reads=self.report_data['footer']['total_reads'],
            footer__threshold_reads=self.report_data['footer']['threshold_reads'],
            footer__threshold_percent=self.report_data['footer']['threshold_percent'],
            footer__negative_control=self.report_data['footer']['negative_control'],
        )

        with fname.with_suffix('.html').open('w') as fout:
            fout.write(rendered)

        HTML(string=rendered).write_pdf(
            fname.with_suffix('.pdf')
        )

        self.logger.info(f'Report written: {fname}')
