import subprocess
import shlex
import sys
import logging
import pyfastx
import random
from pathlib import Path


class PoreLogger:

    def __init__(self, level=logging.ERROR, file: Path = None):

        logging.basicConfig(
            level=level,
            format="[%(asctime)s] [%(name)-7s]     %(message)s",
            datefmt='%H:%M:%S',
            filename=file
        )

        self.logger = logging.getLogger('Sketchy')


def run_cmd(cmd, callback=None, watch=False, background=False, shell=False):

    """Runs the given command and gathers the output.

    If a callback is provided, then the output is sent to it, otherwise it
    is just returned.

    Optionally, the output of the command can be "watched" and whenever new
    output is detected, it will be sent to the given `callback`.

    Returns:
        A string containing the output of the command, or None if a `callback`
        was given.
    Raises:
        RuntimeError: When `watch` is True, but no callback is given.

    """

    if watch and not callback:
        raise RuntimeError(
            "You must provide a callback when watching a process."
        )

    output = None
    try:
        if shell:
            proc = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)
        else:
            proc = subprocess.Popen(shlex.split(cmd), stdout=subprocess.PIPE)

        if background:
            # Let task run in background and return pmid for monitoring:
            return proc.pid, proc

        if watch:
            while proc.poll() is None:
                line = proc.stdout.readline()
                if line != "":
                    callback(line)

            # Sometimes the process exits before we have all of the output, so
            # we need to gather the remainder of the output.
            remainder = proc.communicate()[0]
            if remainder:
                callback(remainder)
        else:
            output = proc.communicate()[0]
    except:
        err = str(sys.exc_info()[1]) + "\n"
        output = err

    if callback and output is not None:
        return callback(output)

    return output


class ArtificialMixture(PoreLogger):

    """ Composes artificial mixtures of host and pathogen sequence reads """

    def __init__(self, composition: dict, reads: int = 10000, verbose: bool = False):

        """ Composes artificial mixtures of host and pathogen reads

        :param composition: dictionary: {'host': {file='host.fq', proportion=0.90}}
        :param reads: total number of reads to compose
        """

        PoreLogger.__init__(self, level=logging.INFO if verbose else logging.ERROR)

        self.composition = composition
        self.reads = reads

        self.check_proportions()

        self.fastq: dict = self.prepare_fastq()

    def prepare_fastq(self) -> dict:

        """ Checks file paths of input files and creates indices """

        fastq = {}
        for organism, data in self.composition.items():
            file = data['file']
            file_path = Path(file)
            if not file_path.exists():
                raise ValueError(f'File {file_path} does not exist.')
            else:
                fastq[organism] = pyfastx.Fastq(file)

        self.logger.info('Prepared read files - proceeding')

        return fastq

    def check_proportions(self):

        """ Check that proportions in composition file sum to 1"""

        proportions = [
            v['proportion'] for k, v in self.composition.items()
        ]

        if sum(proportions) < 1.0:
            raise ValueError('Sum of proportions between host and pathogen must be 1.0.')
        elif sum(proportions) > 1.0:
            raise ValueError('Sum of proportions between host and pathogen allocations cannot exceed 1.0')
        else:
            self.logger.info('Sum of proportions equals 1.0 - proceeding')

    def compose(self, fout: Path, shuffle: bool = True):

        """ Compose an artifical mixture of reads

        Read names / decription headers are renamed according to sequentially numbered
        keys in the composition file, e.. saureus_0, suareus_1 ... to better distinguish
        between composition components later.

        :param fout: file path to output fastq
        :param shuffle: shuffle reads before writing

        """

        self.logger.info('Sample and mix read data')

        reads_out = []
        for organism, fastq in self.fastq.items():
            read_names = [read.name for read in fastq]  # need to solve iterator for sampling, names avoid memory
            sampled_read_names = self.sample(read_names, reads=int(
                self.composition[organism]['proportion']*self.reads)
            )  # check if integer conversion can reduce total reads

            read_strings = self.rename_headers(
                reads=[fastq[name] for name in sampled_read_names],
                organism=organism
            )

            reads_out += read_strings

        if shuffle:
            self.logger.info('Shuffle output reads')
            random.shuffle(reads_out)

        self.logger.info(f'Write reads to: {fout}')
        with fout.open('w') as out:
            for read_str in reads_out:
                out.write(read_str + '\n')

        self.clean()

    def clean(self):

        """ Clean up the Fastq index files from Pyfastx """

        for _, data in self.composition.items():
            index_file = Path(data['file'] + '.fxi')
            if index_file.exists():
                index_file.unlink()

    @staticmethod
    def rename_headers(reads: list, organism: str):

        """ Rename read headers from the Pyfastx reads (read-only) """

        i = 0
        read_strings = []
        for read in reads:
            read_str = read.raw.splitlines()
            read_str[0] = f'@{organism}_{i}'
            read_str = '\n'.join(read_str)
            read_strings.append(read_str)
            i += 1

        return read_strings

    @staticmethod
    def sample(fastq: list, reads: int = None, replacement: bool = False):

        """ Sample a list of Fastq reads / read names """

        if replacement:
            sampled_reads = random.choices(fastq, k=reads)
        else:
            sampled_reads = random.sample(fastq, k=reads)

        return sampled_reads




