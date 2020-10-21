import click
import pandas
from pathlib import Path


@click.command()
@click.option(
    "--vcf", "-v", default=None, help="VCF of core genome SNPs from Snippy", type=Path,
)
@click.option(
    "--out", "-o", default=None, help="VCF for input to candidate calling in Clair", type=Path,
)
def snippy_core_to_clair(vcf, out):

    """ Reformat Snippy Core VCF for input to Clair candidate variant calling """

    # Use a simple DataFrame for this:

    names = ["CHROM", "POS", "ID",  "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"]

    df = pandas.read_csv(
        vcf, sep="\t", comment="#", header=None, names=names, usecols=names
    )

    df['CORE'] = ["1/1" for _ in df.POS]

    header = []
    with vcf.open("r") as header_read:
        for line in header_read:
            if line.startswith("#"):
                if not line.startswith("#CHROM"):
                    header.append(line)

    header = "".join(header)
    column_header = "#" + "\t".join(names) + "\tCORE\n"
    with out.open("w") as clair_out:
        clair_out.write(header+column_header)
        df.to_csv(clair_out, header=False, index=False, sep="\t")
