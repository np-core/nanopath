import click
import pandas
from pathlib import Path


@click.command()
@click.option(
    "--df", "-d", type=Path, help="Data file with header"
)
@click.option(
    "--output", "-o", type=Path, help="Data file output"
)
@click.option(
    "--column", "-c", type=str, default="name", help="Column name to sort by, can be comma delimited string [name]"
)
@click.option(
    "--numeric", "-n", is_flag=True, help="Apply numeric string sorting [false]"
)
@click.option(
    "--sep_out", "-so", type=str, default="\t", help="Data frame delimiter output [\\t]"
)
@click.option(
    "--sep_in", "-si", type=str, default="\t", help="Data frame delimiter input [\\t]"
)
def sort_df(df, output, column, numeric, sep_in, sep_out):

    """ Join two dataframes to add selected trait columns """

    columns = column.split(',')

    df = pandas.read_csv(df, sep=sep_in, header=0)

    print(df)

    df = df.sort_values(by=columns)

    print(df)

    df.to_csv(output, sep=sep_out, index=False)