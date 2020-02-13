import subprocess
import shlex
import sys
import logging
import pyfastx
import pandas
import random
from pathlib import Path
from pysam import AlignmentFile, FastxFile

import tempfile


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


class HostDecontaminator:

    def __init__(self, fastq: Path):

        self.fastq = fastq

    def decontaminate(
        self,
        method='minimap2',
        db: Path = None,
        host_mmi: Path = None,
        microbial_mmi: Path = None,
        mapqual: int = 60,
        threads: int = 3
    ):

        if method == 'minimap2':

            if host_mmi is None or microbial_mmi is None:
                raise ValueError('Need to specify host index file when using method: minimap2')

            # Reads mapping to human reference:
            with tempfile.TemporaryDirectory() as tmpdir:
                tmpdir = Path(tmpdir)

                run_cmd(
                    f'minimap2 -ax map-ont {host_mmi} {self.fastq} -o {tmpdir / "aln.sam"} -t {threads} &&'
                    f'samtools view -Sbq {mapqual} {tmpdir / "aln.sam"} -o {tmpdir / "aln.bam"}'
                )

                with AlignmentFile(f"{tmpdir / 'aln.bam'}", "rb") as bamfile:
                    human_reads_mapped = [a.query_name for a in bamfile]

            # Reads mapping to microbial reference:
            with tempfile.TemporaryDirectory() as tmpdir:
                tmpdir = Path(tmpdir)

                run_cmd(
                    f'minimap2 -ax map-ont {microbial_mmi} {self.fastq} -o {tmpdir / "aln.sam"} &&'
                    f'samtools view -Sbq {mapqual} {tmpdir / "aln.sam"} -o {tmpdir / "aln.bam"}'
                )

                with AlignmentFile(f"{tmpdir / 'aln.bam'}", "rb") as bamfile:
                    microbe_reads_mapped = [a.query_name for a in bamfile]

        elif method == 'minikraken2':

            if db is None:
                raise ValueError('Need to specify database file when using method: minikraken2')

            with tempfile.TemporaryDirectory() as tmpdir:
                tmpdir = Path(tmpdir)

                run_cmd(
                    f'kraken2 --db {db} --threads {threads} --output {tmpdir / "reads.out"}'
                    f'--report {tmpdir / "reads.report"} --use-names --memory-mapping {self.fastq}'
                )


"""
Produces simple Sankey Diagrams with matplotlib.
@author: Anneya Golob & marcomanz & pierre-sassoulas & jorwoods
                      .-.
                 .--.(   ).--.
      <-.  .-.-.(.->          )_  .--.
       `-`(     )-'             `)    )
         (o  o  )                `)`-'
        (      )                ,)
        ( ()  )                 )
         `---"\    ,    ,    ,/`
               `--' `--' `--'
                |  |   |   |
                |  |   |   |
                '  |   '   |
"""

from collections import defaultdict

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns


class PySankeyException(Exception):
    pass


class NullsInFrame(PySankeyException):
    pass


class LabelMismatch(PySankeyException):
    pass


def check_data_matches_labels(labels, data, side):
    if len(labels > 0):
        if isinstance(data, list):
            data = set(data)
        if isinstance(data, pd.Series):
            data = set(data.unique().tolist())
        if isinstance(labels, list):
            labels = set(labels)
        if labels != data:
            msg = "\n"
            if len(labels) <= 20:
                msg = "Labels: " + ",".join(labels) + "\n"
            if len(data) < 20:
                msg += "Data: " + ",".join(data)
            raise LabelMismatch('{0} labels and data do not match.{1}'.format(side, msg))


def sankey(left, right, leftWeight=None, rightWeight=None, colorDict=None,
           leftLabels=None, rightLabels=None, aspect=4, rightColor=False,
           fontsize=14, figureName=None, closePlot=False):
    '''
    Make Sankey Diagram showing flow from left-->right
    Inputs:
        left = NumPy array of object labels on the left of the diagram
        right = NumPy array of corresponding labels on the right of the diagram
            len(right) == len(left)
        leftWeight = NumPy array of weights for each strip starting from the
            left of the diagram, if not specified 1 is assigned
        rightWeight = NumPy array of weights for each strip starting from the
            right of the diagram, if not specified the corresponding leftWeight
            is assigned
        colorDict = Dictionary of colors to use for each label
            {'label':'color'}
        leftLabels = order of the left labels in the diagram
        rightLabels = order of the right labels in the diagram
        aspect = vertical extent of the diagram in units of horizontal extent
        rightColor = If true, each strip in the diagram will be be colored
                    according to its left label
    Ouput:
        None
    '''
    if leftWeight is None:
        leftWeight = []
    if rightWeight is None:
        rightWeight = []
    if leftLabels is None:
        leftLabels = []
    if rightLabels is None:
        rightLabels = []
    # Check weights
    if len(leftWeight) == 0:
        leftWeight = np.ones(len(left))

    if len(rightWeight) == 0:
        rightWeight = leftWeight

    plt.figure()
    plt.rc('text', usetex=False)
    plt.rc('font', family='serif')

    # Create Dataframe
    if isinstance(left, pd.Series):
        left = left.reset_index(drop=True)
    if isinstance(right, pd.Series):
        right = right.reset_index(drop=True)
    dataFrame = pd.DataFrame({'left': left, 'right': right, 'leftWeight': leftWeight,
                              'rightWeight': rightWeight}, index=range(len(left)))

    if len(dataFrame[(dataFrame.left.isnull()) | (dataFrame.right.isnull())]):
        raise NullsInFrame('Sankey graph does not support null values.')

    # Identify all labels that appear 'left' or 'right'
    allLabels = pd.Series(np.r_[dataFrame.left.unique(), dataFrame.right.unique()]).unique()

    # Identify left labels
    if len(leftLabels) == 0:
        leftLabels = pd.Series(dataFrame.left.unique()).unique()
    else:
        check_data_matches_labels(leftLabels, dataFrame['left'], 'left')

    # Identify right labels
    if len(rightLabels) == 0:
        rightLabels = pd.Series(dataFrame.right.unique()).unique()
    else:
        check_data_matches_labels(leftLabels, dataFrame['right'], 'right')
    # If no colorDict given, make one
    if colorDict is None:
        colorDict = {}
        palette = "hls"
        colorPalette = sns.color_palette(palette, len(allLabels))
        for i, label in enumerate(allLabels):
            colorDict[label] = colorPalette[i]
    else:
        missing = [label for label in allLabels if label not in colorDict.keys()]
        if missing:
            msg = "The colorDict parameter is missing values for the following labels : "
            msg += '{}'.format(', '.join(missing))
            raise ValueError(msg)

    # Determine widths of individual strips
    ns_l = defaultdict()
    ns_r = defaultdict()
    for leftLabel in leftLabels:
        leftDict = {}
        rightDict = {}
        for rightLabel in rightLabels:
            leftDict[rightLabel] = dataFrame[(dataFrame.left == leftLabel) & (dataFrame.right == rightLabel)].leftWeight.sum()
            rightDict[rightLabel] = dataFrame[(dataFrame.left == leftLabel) & (dataFrame.right == rightLabel)].rightWeight.sum()
        ns_l[leftLabel] = leftDict
        ns_r[leftLabel] = rightDict

    # Determine positions of left label patches and total widths
    leftWidths = defaultdict()
    for i, leftLabel in enumerate(leftLabels):
        myD = {}
        myD['left'] = dataFrame[dataFrame.left == leftLabel].leftWeight.sum()
        if i == 0:
            myD['bottom'] = 0
            myD['top'] = myD['left']
        else:
            myD['bottom'] = leftWidths[leftLabels[i - 1]]['top'] + 0.02 * dataFrame.leftWeight.sum()
            myD['top'] = myD['bottom'] + myD['left']
            topEdge = myD['top']
        leftWidths[leftLabel] = myD

    # Determine positions of right label patches and total widths
    rightWidths = defaultdict()
    for i, rightLabel in enumerate(rightLabels):
        myD = {}
        myD['right'] = dataFrame[dataFrame.right == rightLabel].rightWeight.sum()
        if i == 0:
            myD['bottom'] = 0
            myD['top'] = myD['right']
        else:
            myD['bottom'] = rightWidths[rightLabels[i - 1]]['top'] + 0.02 * dataFrame.rightWeight.sum()
            myD['top'] = myD['bottom'] + myD['right']
            topEdge = myD['top']
        rightWidths[rightLabel] = myD

    # Total vertical extent of diagram
    xMax = topEdge / aspect

    # Draw vertical bars on left and right of each  label's section & print label
    for leftLabel in leftLabels:
        plt.fill_between(
            [-0.02 * xMax, 0],
            2 * [leftWidths[leftLabel]['bottom']],
            2 * [leftWidths[leftLabel]['bottom'] + leftWidths[leftLabel]['left']],
            color=colorDict[leftLabel],
            alpha=0.99
        )
        plt.text(
            -0.05 * xMax,
            leftWidths[leftLabel]['bottom'] + 0.5 * leftWidths[leftLabel]['left'],
            leftLabel,
            {'ha': 'right', 'va': 'center'},
            fontsize=fontsize
        )
    for rightLabel in rightLabels:
        plt.fill_between(
            [xMax, 1.02 * xMax], 2 * [rightWidths[rightLabel]['bottom']],
            2 * [rightWidths[rightLabel]['bottom'] + rightWidths[rightLabel]['right']],
            color=colorDict[rightLabel],
            alpha=0.99
        )
        plt.text(
            1.05 * xMax,
            rightWidths[rightLabel]['bottom'] + 0.5 * rightWidths[rightLabel]['right'],
            rightLabel,
            {'ha': 'left', 'va': 'center'},
            fontsize=fontsize
        )

    # Plot strips
    for leftLabel in leftLabels:
        for rightLabel in rightLabels:
            labelColor = leftLabel
            if rightColor:
                labelColor = rightLabel
            if len(dataFrame[(dataFrame.left == leftLabel) & (dataFrame.right == rightLabel)]) > 0:
                # Create array of y values for each strip, half at left value,
                # half at right, convolve
                ys_d = np.array(50 * [leftWidths[leftLabel]['bottom']] + 50 * [rightWidths[rightLabel]['bottom']])
                ys_d = np.convolve(ys_d, 0.05 * np.ones(20), mode='valid')
                ys_d = np.convolve(ys_d, 0.05 * np.ones(20), mode='valid')
                ys_u = np.array(50 * [leftWidths[leftLabel]['bottom'] + ns_l[leftLabel][rightLabel]] + 50 * [rightWidths[rightLabel]['bottom'] + ns_r[leftLabel][rightLabel]])
                ys_u = np.convolve(ys_u, 0.05 * np.ones(20), mode='valid')
                ys_u = np.convolve(ys_u, 0.05 * np.ones(20), mode='valid')

                # Update bottom edges at each label so next strip starts at the right place
                leftWidths[leftLabel]['bottom'] += ns_l[leftLabel][rightLabel]
                rightWidths[rightLabel]['bottom'] += ns_r[leftLabel][rightLabel]
                plt.fill_between(
                    np.linspace(0, xMax, len(ys_d)), ys_d, ys_u, alpha=0.65,
                    color=colorDict[labelColor]
                )
    plt.gca().axis('off')
    plt.gcf().set_size_inches(6, 6)
    if figureName != None:
        plt.savefig("{}.png".format(figureName), bbox_inches='tight', dpi=150)
    if closePlot:
        plt.close()