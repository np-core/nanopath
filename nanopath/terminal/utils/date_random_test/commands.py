import click
import pandas
from pandas.errors import EmptyDataError
from pathlib import Path
from nanopath.utils import run_cmd
from nanopath.pathfinder import phybeast_randomise_date_file
from nanopath.pathfinder import phybeast_extract_rate
from matplotlib import pyplot as plt


@click.command()
@click.option(
    "--dates", "-d", default="dates.tsv",
    help="Input tab-delineated table with columns: name, date", type=Path,
)
@click.option(
    "--alignment", "-a", default="dates.core.tsv",
    help="Alignment in FASTA format", type=Path,
)
@click.option(
    "--tree", "-t", default="dates.newick",
    help="Tree in newick format", type=Path,
)
@click.option(
    "--replicates", "-r", default=100,
    help="Estimate distribution of randomly dated molecular clocks across"
         " this many replicates.", type=int,
)
@click.option(
    "--clock_args", default="--covariation --reroot least-squares",
    help="Treetime clock estimate arguments", type=Path,
)
@click.option(
    "--outdir", "-o", default="date_random_test", help="Output directory", type=Path,
)
@click.option(
    "--rate_file", default=None, help="Rate file, one per line", type=Path
)
@click.option(
    "--clock_rate_file", "-c", default=None,
    help="True clock rate, one line in a file", type=Path,
)
@click.option(
    "--limit", '-l', default=300, help="Limit of tries to compute clock estimates", type=int
)
def date_random_test(
    dates, tree, alignment, clock_args, clock_rate_file,
    replicates, outdir, rate_file, limit
):

    """ Date randomisation test using Treetime rate estimates """

    reps = 0  # Treetime can fail with randomised dates!

    outdir.mkdir(parents=True, exist_ok=True)

    if rate_file:
        rates = []
        for line in rate_file.open('r'):
            value = line.rstrip()
            rates.append(
                float(value)
            )
        if clock_rate_file:
            clock_rate = float(
                clock_rate_file.open("r").readline().strip()
            )
        else:
            clock_rate = 0

        fig, axis = plt.subplots(ncols=1, figsize=(15.0, 9.0))
        axis.hist(x=rates, color='gray', bins=50)
        axis.axvline(x=clock_rate, color='r', linewidth=3.0)
        plt.title(
            f'Distribution of date randomisations and true clock rate (red, core genome SNPs: {clock_rate})'
        )

        fig.savefig(outdir / "date_randomisation_test.png")
    else:
        # We use a while counter to list successful runs until the replicate limit
        date_replicates = []
        rate_replicates = []
        while reps < replicates:
            print(f'Randomising dates: {reps}')
            df = phybeast_randomise_date_file(
                date_file=dates, output_file=outdir / f"{reps}.dates.tsv"
            )

            try:
                run_cmd(
                    f"treetime clock --dates {outdir / f'{reps}.dates.tsv'} --tree {tree} "
                    f"--aln {alignment} {clock_args} > {outdir / f'{reps}.rate.txt'}", shell=True
                )
                print(outdir)
                date_replicates.append(df)
            except ValueError:
                print('Failed to run Treetime on randomised dates')
                raise

            phybeast_extract_rate(
                result_file=outdir / f"{reps}.rate.txt",
                prep='treetime',
                output_file=outdir / f"{reps}.rate.tsv"
            )
            try:
                rate = pandas.read_csv(outdir / f"{reps}.rate.tsv")
                rate_replicates.append(rate)
                reps += 1
            except EmptyDataError:
                print('Repeating randomised Treetime clock rate estimate')
                continue

            if reps == limit:
                break

        rate_distribution = pandas.concat(rate_replicates)

        with Path('rates.tsv').open("w") as outfile:
            for r in rate_distribution:
                outfile.write(f"{r}\n")
