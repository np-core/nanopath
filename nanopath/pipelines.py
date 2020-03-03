import json
import logging


from pathlib import Path
from nanopath.utils import PoreLogger
from nanopath.processors import KrakenProcessor, AssemblyProcessor


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
