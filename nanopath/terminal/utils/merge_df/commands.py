
import click
import pandas
from pathlib import Path


@click.command()
@click.option(
    "--data1", "-d1", type=Path, help="Data file 1 with header"
)
@click.option(
    "--data2", "-d2", type=Path, help="Data file 2 with header"
)
@click.option(
    "--output", "-o", default='joined.tsv', help="Joined output dataframe", type=Path
)
@click.option(
    "--how", default='inner', help="Join mode: inner, outer ... [inner]", type=str,
)
@click.option(
    "--left_on", default='name', help="Merge on column in data2 [name]", type=str,
)
@click.option(
    "--right_on", default='name', help="Merge on column in data2 [name]", type=str,
)
@click.option(
    "--sep_in", "-si", default="\t", help="Input delimiter [\\t]", type=str,
)
@click.option(
    "--sep_out", "-so", default="\t", help="Output delimiter [\\t]", type=str,
)
def merge_df(data1, data2, output, how, left_on, right_on, sep_in, sep_out):

    """ Join two dataframes to add selected trait columns """

    df1 = pandas.read_csv(data1, sep=sep_in, header=0)
    df2 = pandas.read_csv(data2, sep=sep_out, header=0)

    print(df1)
    print(df2)

    df = df1.merge(df2, how=how, left_on=left_on, right_on=right_on)

    print(f"D1: {len(df1)} D2: {len(df2)} D: {len(df)}")

    if len(df1) != len(df2):
        diff = set(df1[left_on]).difference(df2[right_on])
        print(f'Difference between data frame indices: {diff}')

    df.to_csv(output, sep=sep_out, index=False)