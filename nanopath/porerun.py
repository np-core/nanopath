""" Module to implement a structure for sequencign runs includign live updates """

from datetime import datetime
from dateutil import parser

from pathlib import Path
from nanopath.utils import PoreLogger

import pandas
import logging


class PoreRun(PoreLogger):

    """ MinKNOW enabled sequencing run directory representation """

    def __init__(
        self, path: Path
    ):
        PoreLogger.__init__(self, level=logging.INFO)

        self.path = path
        self.sample = self.path.parent.name

        self.date, \
            self.rid, \
            self.device, \
            self.flowcell, \
            self.uid \
            = self.path.name.split('_')

        try:
            self.duty_time_file = list(
                self.path.glob('*_duty_time.csv')
            )[0]
        except KeyError:
            self.logger.info(
                f"Could not find 'duty time' file in run: {self.path}"
            )
            self.duty_time_file = None

        try:
            self.throughput_file = list(
                self.path.glob('*_throughput.csv')
            )[0]
        except KeyError:
            self.logger.info(
                f"Could not find 'throughput' file in run: {self.path}"
            )
            self.throughput_file = None

        try:
            self.sequencing_summary_file = list(
                self.path.glob('*_sequencing_summary.txt')
            )[0]
            self.sequencing_summary = self.read_sequencing_summary()
        except KeyError:
            self.logger.info(
                f"Could not find 'sequencing summary' file in run: {self.path}"
            )
            self.sequencing_summary, self.server, self.start_time =  None, None, None

        self.final_summary_file = self.path / 'final_summary.txt'
        if self.final_summary_file.exists():
            self.final_summary = self.read_final_summary()
        else:
            self.logger.info(
                f"Could not find 'final summary' file in run: {self.path}"
            )
            self.final_summary = None

        # Run complete?

        self.complete = True if self.final_summary else False

        # Basecalling paths

        self.fastq_pass = self.path / 'fastq_pass'
        if not self.fastq_pass.exists():
            self.fastq_pass = None

        self.fastq_fail = self.path / 'fastq_fail'
        if not self.fastq_fail.exists():
            self.fastq_fail = None

        self.fast5_pass = self.path / 'fast5_pass'
        if not self.fast5_pass.exists():
            self.fast5_pass = None

        self.fast5_fail = self.path / 'fast5_fail'
        if not self.fast5_fail.exists():
            self.fast5_fail = None

        # Barcodes and Fastq / Fast5 files

        self.barcodes = True

        self.fastq_files = []
        self.fast5_files = []

        self.barcode_fastq = {}
        self.barcode_fast5 = {}

    def read_sequencing_summary(self):

        summary_name = self.sequencing_summary_file.name.split("_")
        self.server, self.start_time = summary_name[0], summary_name[2]

        return pandas.read_csv(
            self.sequencing_summary_file, sep='\t'
        )

    def read_final_summary(self):

        return pandas.read_csv(
            self.sequencing_summary_file, sep='='
        )

    def get_server_data(self):

        pass

    def get_statistics(self):

        """ Parse the sequencing summary file to get current statistics """

        pass
