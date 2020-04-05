""" Server utilities on pathology network """

import os
import json
import shutil
import logging

from nanopath.utils import run_cmd, PoreLogger

from pathlib import Path
from flask_socketio import emit

from pysam import FastxFile, AlignmentFile
from uuid import uuid4


class ServerUtilities(PoreLogger):

    def __init__(self, host_index: Path = None):

        PoreLogger.__init__(self, level=logging.INFO)

        self.host_mmi = Path(
            "/home/esteinig/code/nanopath-app/np_server/data/human38.mmi"
        )

        self.remote_port = 8822
        self.remote_address = "jc225327@zodiac.hpc.jcu.edu.au"

        self.nanopore_runs = Path('/data/nanopath/port_out/sepsis')
        self.workflow_results = Path('/home/esteinig/np-core/test_ground')

        self.local_transfer = Path('/data/nanopath/port_in')
        self.remote_port_in = Path('/data/nanopath/port_in')

    def launch_request_pipeline(
        self,
        fastq: Path = Path("/data/nanopath/test/test_mix.fq"),
        pipeline: str = 'metagenome',
        map_qual: int = 30,
        build_container: bool = False,
        runtime_encryption: bool = False,
        remote_transfer: bool = False,
    ):

        """ Map reads to human reference genome and remove """

        build_dir = Path(f"{uuid4()}")
        build_dir.mkdir(exist_ok=True, parents=True)
        filtered = build_dir / f'{fastq.stem}.microbial.fq'

        # On server (remote) or host machine (local)
        run_cmd(
            f'minimap2 -ax map-ont {self.host_mmi} {fastq} | '
            f'samtools view -Sbq {map_qual} -  > {build_dir / "aln.bam"}', shell=True  # TODO Security shell=True
        )

        with AlignmentFile(f"{build_dir / 'aln.bam'}", "rb") as bamfile:
            to_remove = [a.query_name for a in bamfile]

        with FastxFile(fastq) as fin, \
                filtered.open("w") as fout:
            for read in fin:
                if read.name not in to_remove:
                    fout.write(str(read) + "\n")

        emit(
            'request_launch_request', {'step': 2, 'final_step': False}
        )

        if build_container:

            # TODO encrypt with pub key

            if runtime_encryption:
                pass
            else:
                shutil.copyfile(
                    Path(__file__).parent / 'containers' / 'Light',
                    build_dir / 'Singularity'
                )

            run_cmd(
                f"singularity build --fakeroot "
                f"{build_dir / 'test.sif'} "
                f"{build_dir / 'Singularity'}"
            )

            transfer_file = 'test.sif'
        else:
            transfer_file = filtered

        emit(
            'request_launch_request', {'step': 3, 'final_step': False}
        )

        if remote_transfer:
            run_cmd(
                f"scp -P {self.remote_port} {transfer_file} "
                f"{self.remote_address}:{self.remote_port_in / pipeline}"
            )
        else:
            run_cmd(
                f"cp {transfer_file} {self.local_transfer / pipeline}"
            )

        emit(
            'request_launch_request', {'step': 4, 'final_step': False}
        )

        shutil.rmtree(build_dir)

    def get_analysis_results(
        self, pipeline: str = 'sepsis'
    ):

        self.logger.info(
            f'Parsing results from {self.workflow_results / pipeline}'
        )

        analysis_results = []
        for f in (self.workflow_results / pipeline).glob('*.json'):
            with f.open('r') as server_json:
                result = json.load(server_json)
                analysis_results.append(result)

        # Sort by completed date:

        analysis_results = sorted(
            analysis_results, key=lambda x: x['meta']['completed']
        )

        return analysis_results



