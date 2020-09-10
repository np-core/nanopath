import click
import pandas
from pathlib import Path


@click.command()
@click.option(
    "--data", "-d", type=Path, help="Data file with header"
)
@click.option(
    "--output", "-o", type=Path, help="Data file output"
)
@click.option(
    "--trait", "-t", type=str, default="name", help="Column name to sort by [name]"
)
@click.option(
    "--numeric", "-n", is_flag=True, help="Apply numeric string sorting [false]"
)
@click.option(
    "--sep", "-s", type=str, default="\t", help="Data frame delimiter [\t]"
)
def sort_traits(data, output, trait, numeric, sep):

    """ Join two dataframes to add selected trait columns """

    df = pandas.read_csv(data, sep=sep, header=0)

    print(df)

    df = df.sort_values(by=trait)

    print(df)

    df.to_csv(output, sep=sep, index=False)