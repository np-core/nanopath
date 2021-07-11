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
    "--tail_length", "-t", type=int, default=0,
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

    fidx = 0
    for i in range(nrow):
        for c in range(ncol):
            locus_cov = cov_files[fidx]
            coverage = pandas.read_csv(locus_cov, sep="\t", header=None, names=["locus", 'position', 'coverage'])
            if tail_length > 0:
                coverage = coverage.iloc[tail_length:len(coverage) - tail_length]
            print(coverage)
            sns.lineplot(x=coverage.position, y=coverage.coverage, ax=ax)
            plt.title(locus_cov.stem)
            fidx += 1
    

    plt.tight_layout()
    plt.ylabel("Coverage\n")
    plt.xlabel("\nPosition")

    fig.savefig(f'{plot_file}')

