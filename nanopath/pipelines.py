import json
import logging
import pandas

from pathlib import Path
from nanopath.utils import PoreLogger
from nanopath.processors import KrakenProcessor, AssemblyProcessor


class NanoPathLive(PoreLogger):

    def __init__(self, path: Path, verbose: bool = True):

        PoreLogger.__init__(
            self, level=logging.INFO if verbose else logging.ERROR
        )

        self.path = path

        # Setup workflow output directory path files
        self.guppy_summary_files = list(
            self.path.glob('guppy/*.summary')
        )
        self.guppy_telemetry_files = list(
            self.path.glob('guppy/*.telemetry')
        )
        self.nanoq_stats_files = list(
            self.path.glob('server/*.stats')
        )
        self.kraken_read_files = list(
            self.path.glob('server/*.reads')
        )
        self.kraken_report_files = list(
            self.path.glob('server/*.report')
        )

    def update_run_view(self, **read_distribution_params) -> (list, list, list):

        read_length_dist, read_quality_dist = \
            self.get_read_distributions(**read_distribution_params)

        read_stats = self.get_read_statistics()
        general_telemetry, telemetry = self.get_guppy_telemetry()

        barcode_view_data = []
        for barcode, barcode_data in read_stats.groupby('barcode', sort=False):
            merged = barcode_data.merge(telemetry, on='fast5')
            data = merged[['fast5', 'reads', 'bp', 'time_stamp']] \
                .sort_values('time_stamp')

            barcode_yield_series = self.get_barcode_yield_series(data, telemetry)

            barcode_view_data.append(dict(
                label=barcode,
                reads=int(data.reads.sum()),
                yield_data=barcode_yield_series,
            ))

        return barcode_view_data, read_length_dist.tolist(), read_quality_dist.values.tolist()

    @staticmethod
    def get_barcode_yield_series(
        barcode_data: pandas.DataFrame, telemetry: pandas.DataFrame
    ):
        if len(barcode_data) == len(telemetry):
            # This barcode is present in all basecall files
            return [
                {'x': str(row['time_stamp']), 'y': int(row['bp'])}
                for _, row in barcode_data.iterrows()
            ]
        elif len(barcode_data) < len(telemetry):
            # This barcode is not present in all basecalled files, fill in
            # with previous file value or 0 if absent completely

            yield_series = []
            basepairs = 0
            barcode_data = barcode_data.reset_index()
            for i, time_stamp in enumerate(telemetry.time_stamp):

                if str(time_stamp) in barcode_data.time_stamp.tolist():
                    basepairs += int(
                        barcode_data[barcode_data['time_stamp'] == time_stamp]['bp']
                    )
                    yield_series.append(
                        {'x': str(time_stamp), 'y': basepairs}
                    )
                else:
                    yield_series.append(
                        {'x': str(time_stamp), 'y': basepairs}
                    )
        else:
            raise ValueError

        return yield_series

    def update_stats_view(
        self, read_stats: pandas.DataFrame = None
    ):

        """ Get a side bar view JSON for Dashboard """

        if read_stats is None:
            read_stats = self.get_read_statistics()

        online_view_data = dict(
            basecalled=0,
            reads=int(read_stats.reads.sum()),
            basepairs=int(read_stats.bp.sum()),
            longest=max(read_stats.longest),
            n50=0,
            mean_length=float(read_stats.mean_length.mean()),
            median_length=float(read_stats.mean_length.median()),
            mean_quality=float(read_stats.mean_quality.mean()),
            median_quality=float(read_stats.mean_quality.median()),
        )

        return online_view_data

    def update_telemetry_view(self):

        """ Get a side bar view JSON for Dashboard """

        return self.get_guppy_telemetry()

    def get_read_distributions(
        self,
        length_bin_min: int = 0,
        length_bin_max: int = 80000,
        length_bin_size: int = 1000,
        quality_bin_min: float = 0.0,
        quality_bin_max: float = 21.0,
        quality_bin_size: float = 1.0
    ) -> (pandas.Series, pandas.Series):

        try:
            read_stats_distributions = pandas.concat([
                self.read_summary_file(file) for file in self.guppy_summary_files
            ]).sort_values(['fast5']).set_index('fast5')
        except ValueError:
            raise ValueError(
                'Could not detect required sequencing'
                'summary files from process: Guppy'
            )

        length_bins = pandas.IntervalIndex(
            pandas.interval_range(
                start=length_bin_min,
                end=length_bin_max,
                freq=length_bin_size
            )
        )
        quality_bins = pandas.IntervalIndex(
            pandas.interval_range(
                start=quality_bin_min,
                end=quality_bin_max,
                freq=quality_bin_size
            )
        )

        read_length_bins = pandas.cut(
            read_stats_distributions.read_length, bins=length_bins
        ).value_counts().sort_index()
        read_quality_bins = pandas.cut(
            read_stats_distributions.read_quality, bins=quality_bins
        ).value_counts().sort_index()

        return read_length_bins, read_quality_bins

    def get_guppy_telemetry(self) -> (dict, pandas.DataFrame):
        try:
            guppy_telemetry = pandas.concat([
                self.read_telemetry_file(file) for file in self.guppy_telemetry_files
            ]).sort_values(['fast5']).set_index('fast5')

        except ValueError:
            raise ValueError('Could not detect required telemetry files from process: Guppy')

        return self.process_telemetry(guppy_telemetry)

    def get_read_statistics(self) -> pandas.DataFrame:

        """ Parse basecalling and read statistics for the dashboard run view """
        try:
            read_stats = pandas.concat([
                self.read_nanoq_stats(file) for file in self.nanoq_stats_files
            ]).sort_values(['fast5', 'barcode']).set_index('fast5')
        except ValueError:
            raise ValueError('Could not detect required files belonging to process: Nanoq')

        return read_stats

    @staticmethod
    def process_telemetry(telemetry: pandas.DataFrame) -> (dict, pandas.DataFrame):

        telemetry_processed = {}
        for name in telemetry.columns:
            _unique = telemetry[name].unique().tolist()
            telemetry_processed[name] = _unique

        return telemetry_processed, telemetry.sort_values('time_stamp').reset_index()

    def read_nanoq_stats(self, file: Path):

        df = pandas.read_csv(
            file, sep=' ', names=[
                'reads', 'bp', 'longest', 'shortest', 'mean_length',
                'median_length', 'mean_quality', 'median_quality'
            ]
        )

        return self.add_file_identifiers(df, file)

    def read_summary_file(self, file: Path):

        df = pandas.read_csv(
            file, sep='\t', header=0,
            usecols=['sequence_length_template', 'mean_qscore_template'],
        )
        df = df.rename(
            columns={
                'sequence_length_template': 'read_length',
                'mean_qscore_template': 'read_quality'
            }
        )

        return self.add_file_identifiers(df, file, barcode=False)

    def read_telemetry_file(self, file: Path):

        with file.open('r') as infile:
            telemetry = json.load(infile)

        for seg in telemetry:
            if seg['segment_number'] == 1:  # TODO CHECK THIS
                seg_data = dict(
                    sequencing_kit=seg['context_tags']['sequencing_kit'],
                    experiment_type=seg['context_tags']['experiment_type'],
                    basecall_config_filename=seg['context_tags']['basecall_config_filename'],
                    run_id=seg['run_id'],
                    analysis=seg['software']['analysis'],
                    version=seg['software']['version'],
                    device_type=seg['tracking_id']['device_type'],
                    device_id=seg['tracking_id']['device_id'],
                    exp_start_time=seg['tracking_id']['exp_start_time'],
                    flow_cell_id=seg['tracking_id']['flow_cell_id'],
                    flow_cell_product_code=seg['tracking_id']['flow_cell_product_code'],
                    hostname=seg['tracking_id']['hostname'],
                    operating_system=seg['tracking_id']['operating_system'],
                    protocol_group_id=seg['tracking_id']['protocol_group_id'],
                    sample_id=seg['tracking_id']['sample_id'],
                    time_stamp=seg['tracking_id']['time_stamp']
                )

                df = pandas.DataFrame.from_dict(seg_data, orient='index').T

                return self.add_file_identifiers(df, file, barcode=False)

    def add_file_identifiers(
        self,
        df: pandas.DataFrame,
        file: Path,
        barcode: bool = True
    ) -> pandas.DataFrame:

        if barcode:
            barcode = self.extract_barcode(file)
            df['barcode'] = [barcode for _ in range(len(df))]

        _, _, f5 = self.extract_identifier(file)
        df['fast5'] = [f5 for _ in range(len(df))]

        return df

    @staticmethod
    def extract_barcode(fp: Path):

        return fp.stem.split(".")[-1]

    @staticmethod
    def extract_identifier(fp: Path) -> (str, str, int):

        """ Get file identifier from Guppy: flowcell - run id - fast5 index """

        fc, rid, idx = fp.stem.split(".")[0].split("_")

        return fc, rid, int(idx)


class MetagenomePipeline(PoreLogger):

    """ Metagenome pipeline collection of analysis requests """

    def __init__(self, path: Path, outdir: Path, verbose: bool = True):

        PoreLogger.__init__(
            self, level=logging.INFO if verbose else logging.ERROR
        )

        self.tags = {
            'ont': True,
            'illumina': False,
            'dna': True,
            'rna': False,
            'controls': True
        }

        self.meta = {
            'request': "96ddcaa6-c1ff-41aa-81da-f99a2cb757b2",
            'pipeline': "np-dummy v0.0.0",
            'completed': "20-20-2020 20:20:20",
            'patient': "TSV-001",
            'sample': "Whole Blood",
            'protocol': "qg-dummy v0.0.0"
        }

        self.path = path
        self.outdir = outdir

        self.outdir.mkdir(parents=True, exist_ok=True)

        self.path_kraken = 'kraken'
        self.path_assembly = 'assembly'

        self.results = [
            p for p in self.path.glob('*') if p.is_dir()
        ]

        self.reports = {}

    def collect(self, request: str = None):

        """ Collect the results of the metagenome pipeline using processors """

        if request is None:
            results = self.results
        else:
            results = [self.path / request]

        for path in results:
            request_uuid = path.name

            kraken_server = self.collect_kraken(
                result_path=path, request=request_uuid
            )

            assembly_server = self.collect_assembly(
                result_path=path, request=request_uuid
            )

            self.reports[request_uuid] = self.generate_report(
                reads=kraken_server,
                assemblies=assembly_server,
                viruses=dict(),
                genotypes=dict(),
                controls=dict()
            )

        for req in self.reports.keys():
            self.write_report(request=req)

    def write_report(self, request: str):

            try:
                report = self.reports[request]
            except KeyError:
                raise KeyError(f'Could not find report for request: {request}')

            server_json = self.outdir / f'{request}.json'
            with server_json.open('w') as server_out:
                json.dump(report, server_out)

    def generate_report(
        self,
        reads: dict,
        assemblies: dict,
        viruses: dict,
        genotypes: dict,
        controls: dict
    ) -> dict:

        return {
            'tags': self.tags,
            'meta': self.meta,
            'data': {
                'reads': reads,
                'assemblies': assemblies,
                'viruses': viruses,
                'genotypes': genotypes,
                'controls': controls
            }
        }

    def collect_assembly(self, result_path: Path, request: str) -> dict:

        assembly_path = result_path / self.path_assembly
        if assembly_path.exists():
            ap = AssemblyProcessor(
                fasta=assembly_path / f'{request}.consensus.fasta',
                info_file=assembly_path / f'{request}.info.txt'
            )
            return ap.get_server_data()
        else:
            return dict()

    def collect_kraken(self, result_path: Path, request: str) -> dict:

        kraken_path = result_path / self.path_kraken
        if kraken_path.exists():
            kp = KrakenProcessor(
                report=kraken_path / f'{request}.report',
                reads=kraken_path / f'{request}.reads'
            )
            return kp.get_server_data()
        else:
            return dict()
