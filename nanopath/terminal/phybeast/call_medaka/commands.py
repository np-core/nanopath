import click

from pathlib import Path
from nanopath.processors import MedakaCore


@click.command()
@click.option(
    '--medaka_directory',
    type=Path,
    help='Path to directory containing output '
         'directories with VCFs from Megalodon'
)
@click.option(
    '--snippy_directory',
    type=Path,
    help='Path to directory containing output VCFs from Snippy'
)
@click.option(
    '--reference',
    type=Path,
    help='Reference genome used in Snippy and Megalodon [none]'
)
@click.option(
    '--prefix',
    type=str,
    default='medaka.hybrid',
    help='Prefix for main output files [megalodon]'
)
@click.option(
    '--outdir',
    type=Path,
    default='medaka_core',
    help='Path to output directory [medaka_core]'
)
@click.option(
    '--min_depth',
    type=int,
    default=10,
    help='Minimum read depth to include reference site in ONT alignment [10]'
)
@click.option(
    '--min_quality',
    type=int,
    default=30,
    help='Minimum Phred quality of genotype to include in Medaka [30]'
)
@click.option(
    '--include_reference',
    is_flag=True,
    help='Include the reference sequence in variant alignments [false]'
)
def call_medaka(
    medaka_directory,
    snippy_directory,
    reference,
    include_reference,
    min_quality,
    min_depth,
    outdir,
    prefix,
):

    """ Detect core SNPs from Medaka and Snippy """

    sc = MedakaCore(reference=reference, outdir=outdir, prefix=prefix)

    if medaka_directory:
        sc.parse_medaka(
            medaka_directory=medaka_directory,
            min_quality=min_quality,
            min_depth=min_depth
        )

    sc.parse_snippy(
        snippy_directory=snippy_directory
    )

    sc.call_hybrid_core(include_reference=include_reference)

    #
    # mc = MedakaCore()
    # mc.call_hybrid_core(
    #     snippy_directory=snippy_directory,
    #     reference=reference
    # )




    # sc.check_megalodon(max_missingness=max_missing_rate)
    #
    # sc.create_alignment(
    #     fasta=outdir / f'{prefix+"." if prefix else ""}core.fasta',
    #     merge=merge
    # )
    # if filter_invariant:
    #     sc.filter_invariant_sites(
    #         alignment=outdir / f'{prefix+"." if prefix else ""}core.fasta',
    #         fasta=outdir / f'{prefix+"." if prefix else ""}core.variants.fasta'
    #     )
    # if full and reference:
    #     sc.create_full_alignment(
    #         merge=merge,
    #         reference=reference,
    #         fasta=outdir / f'{prefix+"." if prefix else ""}full.fasta')
