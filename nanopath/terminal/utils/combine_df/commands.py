import click
import pandas
from pathlib import Path


@click.command()
@click.option(
    "--df", "-df", default=None, help="Input dataframes for combining", type=str,
)
@click.option(
    "--dir", "-d", default=None, help="Input directory containing dataframes for globbing", type=Path,
)
@click.option(
    "--glob", "-g", default=None, help="Input glob to collect data frames [*.csv]", type=str,
)
@click.option(
    "--output", "-o", default="combined.tsv", help="Output data frame [combined.tsv]", type=Path,
)
@click.option(
    "--axis", "-a", default=0, help="Axis to combine: 0 (join new rows), 1 (join new columns) [0]", type=int,
)
@click.option(
    "--sep_in", "-si", default="\t", help="Input delimiter [\\t]", type=str,
)
@click.option(
    "--sep_out", "-so", default="\t", help="Output delimiter [\\t]", type=str,
)
@click.option(
    "--extract", "-ex", default="",
    help="Extract categorical data to split from filename, replace this <str> then split file stem",
    type=str
)
@click.option(
    "--extract_split", "-es", default="", help="Extract categorical data split delimiter [_]", type=str
)
@click.option(
    "--extract_head", "-eh", default="", help="Extract categorical data header names [nf_id,model]", type=str
)
@click.option(
    "--clean", "-c", help="Subset model column if present to the first (single) category and remove model column", is_flag=True
)
@click.option(
    "--fill", "-f", help="Fill missing values with this value [raw]", type=str, default=None
)
def combine_df(df, dir, glob, output, sep_in, sep_out, axis, extract, extract_split, extract_head, clean, fill):

    """ Concatenate data frames with the same columns FASTA headers """

    if df is None:
        files = dir.glob(glob)
    else:
        files = df.split(',')

    dfs = []
    for f in files:
        df = pandas.read_csv(f, sep=sep_in)

        # Extract for Nextflow
        if extract:
            header = extract_head.split(",")
            cats = f.name.replace(extract, "").strip().split(extract_split)
            for i, cat in enumerate(cats):
                cat_name = header[i]  # must be same length as cats
                df[cat_name] = [cat for _ in range(len(df))]

        dfs.append(df)

    print(dfs)

    dfout = pandas.concat(dfs, axis=axis)

    if output.parent.exists():
        output.parent.mkdir(parents=True, exist_ok=True)

    print(dfout)

    if clean and 'model' in dfout.columns.tolist():
        extract_model = dfout['model'].unique().tolist()[0]
        dfout = dfout[dfout['model'] == extract_model]
        dfout = dfout.drop(columns="model")

    if fill is not None:
        dfout = dfout.fillna(fill)

    dfout.to_csv(f"{output}", sep=sep_out, index=False)
