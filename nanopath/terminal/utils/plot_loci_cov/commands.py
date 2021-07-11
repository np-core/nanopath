import click
import pandas
from pathlib import Path
from matplotlib import pyplot as plt
import seaborn as sns
from rich import Console

@click.command()
@click.option(
    "--cov_path", "-c", type=Path, help="Locus specific coverage files from BED file loop directory of .cov files"
)
@click.option(
    "--snps", "-s", type=Path, default=None,
    help="SNP coordinate file with columns: snp, bp and chr [snps.txt]"
)
@click.option(
    "--plot_file", "-p", type=Path, default="multi_locus_coverage.pdf",
    help="Plot output file [multi_locus_coverage.pdf]"
)
@click.option(
    "--tail_length", "-t", type=int, default=0,
    help="Additional side sequence coverag added in samtools depth step at each side of target seq [5000 bp]"
)
@click.option(
    "--bar", "-b", is_flag=True,
    help="Plot bar plot instead of line plot (slow)"
)
def plot_loci_cov(
    cov_path, plot_file, tail_length, snps, bar
):

    """ Plot a multiple enriched loci from an adaptive smpling run """

    cov_files = list(cov_path.glob("*.cov"))

    nrow, ncol = len(cov_files)//2, 2

    fig, ax = plt.subplots(
        nrows=nrow, ncols=ncol, figsize=(
            nrow*4.5, ncol*12
        )
    )

    if snps:
        snps = pandas.read_csv(snps, sep="\t", header=0)

    console = Console()

    fidx = 0
    for i in range(nrow):
        for c in range(ncol):
            locus_cov = cov_files[fidx]
            coverage = pandas.read_csv(locus_cov, sep="\t", header=None, names=["locus", 'position', 'coverage'])
            if tail_length > 0:
                coverage = coverage.iloc[tail_length:len(coverage) - tail_length]
            if bar:
                p = sns.barplot(x=coverage.position, y=coverage.coverage, color='gray', ax=ax[i][c])
            else:
                p = sns.lineplot(x=coverage.position, y=coverage.coverage, color='gray', ax=ax[i][c])
            p.set_xticks([])
            p.yaxis.get_major_locator().set_params(integer=True)
            if snps is not None:
                chrom_contig = '_'.join(locus_cov.stem.split("_")[0:2])  # must be chromosome contig name NC_
                snp = snps[snps['chr'] == chrom_contig].values[0]
                p.set_title(f"{snp[0]} @ {snp[2]}")
                if not bar:
                    p.axvline(x=int(snp[1]), color='r')
                console.print(
                    f"{snps[0]:<20}@{snps[2]:<8}@{snps[1]:<20}  A1: {snps[3]:<5} Odds: {snps[5]:<7}"
                )
            else:
                p.set_title(locus_cov.stem)
            fidx += 1


    plt.tight_layout()
    plt.ylabel("Coverage\n")
    plt.xlabel("\nPosition")

    fig.savefig(f'{plot_file}')

