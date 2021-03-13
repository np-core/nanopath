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
    '--exclude_genotype',
    '-g',
    default=None,
    type=str,
    help='List of comma separated columns of genotypes to exclude'
)
@click.option(
    '--outdir',
    '-o',
    default='collected',
    type=Path,
    help='Path to server collection output directory'
)
def collect(path, exclude, exclude_genotype, outdir):
    """ Collect output from assembly pipeline into a data table """

    ap = AssemblyPipeline(path=path, outdir=outdir)

    if exclude_genotype is not None:
        exclude_genotype = exclude_genotype.split(',')

    ref = ap.collect_genotypes(component='illumina', exclude=exclude_genotype)
    ont = ap.collect_genotypes(component='ont', exclude=exclude_genotype)
    hybrid = ap.collect_genotypes(component='hybrid', exclude=exclude_genotype)
    unicycler = ap.collect_genotypes(component='unicycler', exclude=exclude_genotype)

    nanoq = ap.collect_statistics(mode="")

    if exclude:
        excl = pandas.read_csv(
            exclude, header=None, names=['name']
        ).name.tolist()
    else:
        excl = None

    print(unicycler)

    if ont is not None or hybrid is not None or unicycler is not None:
        ap.plot_genotype_heatmap(
            reference=ref, genotypes={'ont': ont, 'hybrid': hybrid, 'unicycler': unicycler}, exclude=excl
        )

    dnadiff = ap.collect_dnadiff(exclude=excl)

    if dnadiff is not None:
        ont_dnadiff = dnadiff.loc[dnadiff['branch'] == 'ont_medaka', :]\
            .merge(ont, on="name", how='inner').merge(nanoq, on="name")

        hybrid_dnadiff = dnadiff.loc[dnadiff['branch'] == 'hybrid_medaka', :]\
            .merge(hybrid, on="name").merge(nanoq, on="name")

        unicycler_dnadiff = dnadiff.loc[dnadiff['branch'] == 'hybrid_unicycler', :] \
            .merge(ont, on="name", how='inner').merge(nanoq, on="name")

        ont_dnadiff.to_csv(
            outdir / 'ont_vs_ref.tsv', sep='\t', index=False
        )

        hybrid_dnadiff.to_csv(
            outdir / 'hybrid_vs_ref.tsv', sep='\t', index=False
        )

        unicycler_dnadiff.to_csv(
            outdir / 'unicycler_vs_ref.tsv', sep='\t', index=False
        )

    ref.to_csv(
        outdir / 'illumina_reference_genotypes.tsv', sep='\t', index=False
    )

    nanoq.to_csv(
        outdir / 'read_qc_all.tsv', sep='\t', index=False
    )

