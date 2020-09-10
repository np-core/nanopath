import click
import pandas

from pathlib import Path


@click.command()
@click.option(
    "--data", "-d", help="Data file with columns name and and to subset", type=Path,
)
@click.option(
    "--column", "-c", help="Column with subset identifiers", type=str, default='date'
)
def show_dates(data, column):

    """ Remove 'Reference' from Snippy alignment output file """

    df = pandas.read_csv(data, sep='\t', header=0, na_values=["-"])

    dates = df[column]

    print(min(dates), max(dates), max(dates) - min(dates), len(df))
