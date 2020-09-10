import click
import pandas
import seaborn as sns
from pathlib import Path
from nanopath.utils import color_tree


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
    "--traits", default='name,date', help="List of trait columns to retain [name,date]", type=str,
)
@click.option(
    "--sep", default='\t', help="Delimiter for input dataframes", type=str
)
def join_traits(data1, data2, output, how, left_on, right_on, traits, sep):

    """ Join two dataframes to add selected trait columns """

    df1 = pandas.read_csv(data1, sep=sep, header=0)
    df2 = pandas.read_csv(data2, sep=sep, header=0)

    print(df1)
    print(df2)

    df = df1.merge(df2, how=how, left_on=left_on, right_on=right_on)

    traits = traits.split(',')

    for c in df.columns:
        if c not in traits:
            df.drop(columns=c, inplace=True)

    print(f"D1: {len(df1)} D2: {len(df2)} D: {len(df)}")

    if len(df1) != len(df2):
        diff = set(df1[left_on]).difference(df2[right_on])
        print(f'Difference between data frame indices: {diff}')

    df.to_csv(output, sep=sep, index=False)