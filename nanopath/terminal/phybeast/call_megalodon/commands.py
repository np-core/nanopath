import click

from pathlib import Path
from nanopath.processors import MegalodonCore


@click.command()
@click.option(
    '--megalodon_directory',
    type=Path,
    help='Path to directory containing output '
         'directories with VCFs from Megalodon'
)
@click.option(
    '--vcf_snippy',
    type=Path,
    help='Path to core variant output VCF file from Snippy'
)
@click.option(
    '--vcf_megalodon',
    type=Path,
    help='Path to single output VCF file from Megalodon'
)
@click.option(
    '--merge',
    is_flag=True,
    help='Merge the Megalodon calls with Snippy calls [false]'
)
@click.option(
    '--full',
    is_flag=True,
    help='Write the full alignment with substituted variants [false]'
)
@click.option(
    '--reference',
    type=Path,
    help='Reference genome for full alignment used in Snippy and Megalodon [none]'
)
@click.option(
    '--filter_invariant',
    is_flag=True,
    help='Filter invariant sites sites from core alignment [false]'
)
@click.option(
    '--outdir',
    type=Path,
    default='megalodon_core',
    help='Path to output directory [megalodon_core]'
)
@click.option(
    '--min_depth',
    type=int,
    default=5,
    help='Minimum read depth at genotype site to include [5]'
)
@click.option(
    '--min_likelihood',
    type=int,
    default=30,
    help='Minimum phred-scaled likelihood of genotype to include (0 - 999) [30]'
)
@click.option(
    '--max_missing_rate',
    type=float,
    default=0.05,
    help='Filter out samples with > missing rate [0.05]'
)
@click.option(
    '--prefix',
    type=str,
    default='megalodon',
    help='Prefix for main output files [megalodon]'
)
def call_megalodon(
    megalodon_directory,
    vcf_megalodon,
    vcf_snippy,
    merge,
    filter_invariant,
    full,
    reference,
    outdir,
    min_depth,
    min_likelihood,
    max_missing_rate,
    prefix,
):

    """ Merge variants detected by Megalodon with variants from Snippy """

    # TODO - check if there is discrepancy between Snippy VCF and Alignments
    # TODO - if so, use the VCF to generate the merged output alignments

    outdir.mkdir(parents=True, exist_ok=True)

    sc = MegalodonCore(core_vcf=vcf_snippy)

    if megalodon_directory:
        vcf_files = sorted([
            d / 'variants.sorted.vcf'
            for d in megalodon_directory.glob('*')
            if d.is_dir()
        ])
    else:
        vcf_files = [vcf_megalodon]

    for file in vcf_files:
        try:
            sc.parse_megalodon(
                vcf=file,
                min_depth=min_depth,
                min_quality=min_likelihood
            )
        except FileNotFoundError:
            sc.logger.info(
                f'Could not detect file: {file.absolute()}. Skipping.'
            )

    sc.check_megalodon(max_missingness=max_missing_rate)

    sc.create_alignment(
        fasta=outdir / f'{prefix+"." if prefix else ""}core.fasta',
        merge=merge
    )
    if filter_invariant:
        sc.filter_invariant_sites(
            alignment=outdir / f'{prefix+"." if prefix else ""}core.fasta',
            fasta=outdir / f'{prefix+"." if prefix else ""}core.variants.fasta'
        )
    if full and reference:
        sc.create_full_alignment(
            reference=reference,
            fasta=outdir / f'{prefix+"." if prefix else ""}full.fasta')
