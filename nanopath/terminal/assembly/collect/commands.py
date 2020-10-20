import click
import pandas
from nanopath.pipelines import AssemblyPipeline
from pathlib import Path


@click.command()
@click.option(
    '--path',
    '-p',
    type=Path,
    help='Path to pipeline output directory containing results'
)
@click.option(
    '--exclude',
    '-e',
    default=None,
    type=Path,
    help='Path to file with one name per line to exclude (e.g. from manual curation)'
)
@click.option(
    '--outdir',
    '-o',
    default='collected',
    type=Path,
    help='Path to server collection output directory'
)
def collect(path, exclude, outdir):
    """ Collect output from assembly pipeline into a data table """

    ap = AssemblyPipeline(path=path, outdir=outdir)

    ref = ap.collect_genotypes(component='illumina')
    ont = ap.collect_genotypes(component='ont')
    hybrid = ap.collect_genotypes(component='hybrid')

    nanoq = ap.collect_statistics(mode="")

    print(nanoq)

    if exclude:
        excl = pandas.read_csv(
            exclude, header=None, names=['name']
        ).name.tolist()
    else:
        excl = None

    if ont is not None and hybrid is not None:
        ap.plot_genotype_heatmap(
            genotype=ont, genotype2=hybrid, reference=ref, exclude=excl
        )

    dnadiff = ap.collect_dnadiff(exclude=excl)

    if dnadiff is not None:
        ont_dnadiff = dnadiff.loc[dnadiff['branch'] == 'ont_medaka', :]\
            .merge(ont, on="name", how='inner').merge(nanoq, on="name")

        hybrid_dnadiff = dnadiff.loc[dnadiff['branch'] == 'hybrid_medaka', :]\
            .merge(hybrid, on="name").merge(nanoq, on="name")

        ont_dnadiff.to_csv(
            outdir / 'ont_vs_ref.tsv', sep='\t', index=False
        )

        hybrid_dnadiff.to_csv(
            outdir / 'hybrid_vs_ref.tsv', sep='\t', index=False
        )

    ref.to_csv(
        outdir / 'illumina_reference_genotypes.tsv', sep='\t', index=False
    )

