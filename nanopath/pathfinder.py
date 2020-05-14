import subprocess
import dendropy
import pandas
import pysam
import shlex
import sys
import os

from random import shuffle
from pathlib import Path
from nanopath.plots import plot_date_randomisation
from nanopath.plots import plot_date_regression

import matplotlib.pyplot as plt


def run_cmd(cmd, callback=None, watch=False, background=False, shell=False):

    """ Runs the given command and gathers the output.

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

############################
# Survey Support Functions #
############################


def get_aspera_key() -> Path:

    """ Return path to aspera key for connections to ENA """

    return Path(__file__).parent / "resources" / "asperaweb_id_dsa.openssh"


def get_genome_sizes() -> pandas.DataFrame:

    """ Return a dataframe from the `genome.size` file in `pathfinder.resources`

    Genome sizes are computed as media genome size for given taxonomic
    identifier from the NCBI Prokaryot DB.

    :return Dataframe with one column size and row index taxid

    """

    genome_sizes = Path(__file__).parent / "resources" / "genome.sizes"
    genome_sizes = pandas.read_csv(genome_sizes, index_col=0)

    return genome_sizes

###############################
# Alignment support functions #
###############################


def remove_sample(alignment: Path, outfile: Path, remove: str or list) -> None:

    """ Remove any sequence from the alignment file by sequence names

    :param alignment: alignment file (.fasta)
    :param outfile: output file (.fasta)
    :param remove: sequence identifiers to remove

    :return:  None, outputs alignment file with sequences removed

    """

    if isinstance(remove, str):
        remove = [remove]

    with pysam.FastxFile(alignment) as fin, outfile.open('w') as fout:
            for entry in fin:
                if entry.name not in remove:
                    fout.write(str(entry) + '\n')

###################################
# Phylogenetics Support Functions #
###################################


def get_tree_dates(newick_file: Path) -> pandas.DataFrame:

    """ Get the leaf names and dates from the input tree

    :param newick_file: tree file in newick format

    :returns `pandas.DataFrame` with two columns: name, date

    """

    tree = dendropy.Tree.get(path=newick_file, schema="newick")

    return pandas.DataFrame(
        data=[taxon.label.split() for taxon in tree.taxon_namespace],
        columns=['name', 'date']
    )

##############################
# Phybeast Support Functions #
##############################


def phybeast_randomise_date_file(
    date_file: Path,
    output_file: Path = None
) -> pandas.DataFrame:

    """ Randomise order of dates in file

    DataFrame input options can be passed as **kwargs

    :param date_file: path to date file with columns: name and date
    :param output_file: path to tab-delimited output file for randomised dates

    :returns DataFrame with shuffled dates

    :raises ValueError if date and name not in column headers

    """

    df = pandas.read_csv(date_file, sep='\t')

    if 'date' not in df.columns or 'name' not in df.columns:
        raise ValueError('Could not find date and name in columns')

    # Suppress warning
    with pandas.option_context('mode.chained_assignment', None):
        dates = df.date
        shuffle(dates)
        df.date = dates

    if output_file is not None:
        df.to_csv(output_file, sep='\t', header=True, index=False)

    return df


def phybeast_prepare_metadata_file(
    meta_file: Path,
    prep: str = 'lsd2',
    output_file: Path = Path.cwd() / 'lsd2.meta',
    na_range: bool = True,
) -> None:

    """ Prepare the tab-delimited input file for pf-phybeast

    :param meta_file: tab-delimited meta data file with columns: name, date
    :param prep: output file type to prepare for: lsd2, treetime
    :param output_file: output file path
    :param na_range: replace missing dates with range between minimum and maximum date

    :returns None, writes to file :param out_file

    """

    df = pandas.read_csv(meta_file, sep='\t')

    if 'date' not in df.columns or 'name' not in df.columns:
        raise ValueError('Could not find date and name in columns')

    if prep == 'treetime':
        df.to_csv(output_file, header=True, sep=',', index=False)

    if prep == 'lsd2':
        with output_file.open('w') as outfile:
            outfile.write(
                f'{len(df)}\n'
            )

            if na_range:
                na_rep = f'b({min(df["date"])},{max(df["date"])})'
            else:
                na_rep = 'NA'

            # TODO: document no spaces in unique ids
            df[['name', 'date']].to_csv(
                outfile, header=False, sep=' ', index=False, na_rep=na_rep
            )


def phybeast_extract_rate(
    result_file: Path,
    prep: str = 'lsd2',
    output_file: Path = Path.cwd() / 'rate.txt'
) -> None:

    """ Prepare the tab-delimited input file for pf-phybeast

    :param result_file: path to summary output file from: lsd2 or treetime
    :param prep: output file type to prepare for: lsd2 or treetime
    :param output_file: Output file path

    :returns None, writes to file :param output_file

    """

    if prep == 'lsd2':

        with result_file.open('r') as infile, output_file.open('w') as outfile:
            for line in infile:
                if line.strip().startswith('rate'):
                    rate = float(
                        line.strip().split()[1].strip(',')
                    )
                    tmrca = float(
                        line.strip().split()[3].strip(',')
                    )
                    print(rate, tmrca)
                    outfile.write(f"{rate}\t{tmrca}\n")

    elif prep == 'treetime':

        with result_file.open('r') as infile, output_file.open('w') as outfile:
            for line in infile:
                if line.strip().startswith('--rate:'):
                    rate = float(
                        line.strip().split()[1]
                    )
                    print(rate)
                    outfile.write(f"{rate}\n")


def phybeast_plot_date_randomisation(
    replicate_file: Path,
    rate_file: Path,
    tree_file: Path,
    date_file: Path,
    tree_format: str = 'nexus',
    output_file: Path = Path("date_randomisation.png"),
) -> None:

    """ Plot distribution of date randomised substitution rates

    :param replicate_file: `rates.tab` with replicates from `DateRandomisationPlot`
    :param rate_file: `rate.txt` containing true rate from `MolecularClock`
    :param tree_file:
    :param date_file:
    :param tree_format:

    :param output_file: output plot file, format by extension

    :returns None, writes to file :param output_file

    """

    fig, axes = plt.subplots(ncols=2, figsize=(27.0, 9))
    ax1, ax2 = axes.flatten()  # observe order

    # Date Randomisation

    replicate_df = pandas.read_csv(
        replicate_file, sep='\t', names=['replicate', 'tmrca']
    )

    replicates = replicate_df.replicate.tolist()

    rate_df = pandas.read_csv(
         rate_file, sep='\t', names=['rate', 'tmrca']
    )

    rate = float(rate_df.iloc[0, 0])

    plot_date_randomisation(
        ax=ax1, replicates=replicates, rate=rate
    )

    # Regression

    tree = dendropy.Tree.get(path=tree_file, schema=tree_format)

    distances = []
    tips = []
    for tip in tree.leaf_node_iter():
        distances.append(
            tip.distance_from_root()
        )
        tips.append(
            str(tip.taxon).replace("'", "")
        )

    distance_df = pandas.DataFrame(
        {'name': tips, 'distance': distances}
    )

    date_df = pandas.read_csv(date_file, sep='\t')

    print(date_df)
    print(distance_df)

    regression_data = date_df.merge(distance_df).dropna(axis=0)

    print(regression_data)

    plot_date_regression(
        ax=ax2, x=regression_data.date.values, y=regression_data.distance.values
    )

    fig.savefig(output_file)

#####################################
# Pipeline Modules Helper Functions #
#####################################


def get_content(path):
    """Get content by directory / files if it contains files """
    content = {}
    for directory, subdirs, files in os.walk(path):
        if files:
            content[Path(directory)] = [
                Path(directory) / Path(file) for file in files
            ]

    return content


def get_id_from_fname(fname: str or Path, remove: str or list = None):
    """Helper function to deconstruct a filename
    since there is no scheme to IDs.
    :param fname: file basename.
    :param remove: substrings to be removed."""

    if not remove:
        return fname

    if isinstance(fname, Path):
        fname = str(fname.name)
    if isinstance(remove, list):
        for r in remove:
            fname = fname.replace(r, '')
        return fname
    elif isinstance(remove, str):
        return fname.replace(remove, '')
    else:
        raise ValueError


def get_subdict(key, dictionary):

    for k, v in dictionary.items():
        if k == key:
            yield v
        elif isinstance(v, dict):
            for result in get_subdict(key, v):
                yield result
        elif isinstance(v, list):
            for d in v:
                if isinstance(d, dict):
                    for result in get_subdict(key, d):
                        yield result


def retain_files(files: list, retain: str) -> list:
    """Helper function to retain files if sub str
    is contained in file name
    :param files: list of file paths
    :param retain: str in file name, to retain file
    """
    retained = []
    for file in files:
        if isinstance(file, Path):
            fname = str(file.name)
        else:
            fname = os.path.basename(file)
        if retain in fname:
            retained.append(file)

    return retained
