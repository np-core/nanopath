import click
import pandas

from pathlib import Path

@click.command()
@click.option(
    '--data',
    '-d',
    type=Path,
    help='Data with columns: name,trait1,trait2 ...'
)
@click.option(
    '--alignment',
    '-a',
    type=Path,
    help='Path to alignment file to subset data with names'
)
@click.option(
    '--trait',
    '-r',
    type=str,
    help='Trait column in data to count values for'
)
def count_trait(data, alignment, trait):

    """ Quickly count unique values in a data file column """

    df = pandas.read_csv(data, sep="\t", header=0)

    if not alignment:
        print(df[trait].value_counts())
    else:
        with alignment.open('r') as aln:
            names = [line.strip().split().replace(">", "") for line in aln if line.startswith(">")]

        df = df[df['name'].isin(names)]

        print(df[trait].value_counts())
