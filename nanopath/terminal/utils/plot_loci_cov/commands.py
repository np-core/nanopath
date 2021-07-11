import click
import pandas
from pathlib import Path
from matplotlib import pyplot as plt
import seaborn as sns


@click.command()
@click.option(
    "--cov_path", "-c", type=Path, help="Locus specific coverage files from BED file loop directory of .cov files"
)
@click.option(
    "--plot_file", "-p", type=Path, default="multi_locus_coverage.pdf",
    help="Plot output file [multi_locus_coverage.pdf]"
)
@click.option(
    "--tail_length", "-t", type=int, default=5000,
    help="Additional side sequence coverag added in samtools depth step at each side of target seq [5000 bp]"
)
def plot_loci_cov(
    cov_path, plot_file, tail_length
):

    """ Plot a multiple enriched loci from an adaptive smpling run """

    cov_files = list(cov_path.glob("*.cov"))

    nrow, ncol = len(cov_files)//2, 2

    fig, ax = plt.subplots(
        nrows=nrow, ncols=ncol, figsize=(
            nrow*7, ncol*4.5
        )
    )

    for locus in cov_files:
        coverage = pandas.read_csv(locus, sep="\t", header=None, names=["locus", 'position', 'coverage'])
        if tail_length > 0:
            coverage = coverage.iloc[tail_length:len(coverage)-tail_length]
        sns.barplot(data=coverage, x='position', y="coverage", legend=False, ax=ax)
        plt.title(locus.stem)

    plt.tight_layout()
    plt.ylabel("Coverage\n")
    plt.xlabel("\nPosition")

    fig.savefig(f'{plot_file}')

