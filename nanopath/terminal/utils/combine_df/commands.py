import click
import pandas
from pathlib import Path


@click.command()
@click.option(
    "--dir", "-d", default=".", help="Input directory containing dataframes for globbing", type=Path,
)
@click.option(
    "--glob", "-g", default="*.csv", help="Input glob to collect data frames [*.csv]", type=str,
)
@click.option(
    "--output", "-o", default="combined.tsv", help="Output data frame [combined.tsv]", type=Path,
)
@click.option(
    "--axis", "-a", default=0, help="Axis to comine: 0 columns (join new rows), 1 rows (join new columns) [0]", type=int,
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
def combine_df(dir, glob, output, sep_in, sep_out, axis, extract, extract_split, extract_head):

    """ Concatenate data frames with the same columns FASTA headers """

    dfs = []
    for f in dir.glob(glob):
        df = pandas.read_csv(f, sep=sep_in)

        # Extract for Nextflow
        if extract:
            header = extract_head.split(",")
            cats = f.stem.replace(extract, "").strip().split(extract_split)
            for i, cat in enumerate(cats):
                cat_name = header[i]  # must be same length as cats
                df[cat_name] = [cat for _ in range(len(df))]

        dfs.append(df)

    print(dfs)

    dfout = pandas.concat(dfs, axis=axis)

    if output.parent.exists():
        output.parent.mkdir(parents=True, exist_ok=True)

    print(dfout)

    dfout.to_csv(f"{output}", sep=sep_out)
