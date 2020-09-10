import click
import pandas

from nanopath.processors import SnippySample, MedakaSample
from pathlib import Path


@click.command()
@click.option(
    '--dir_snippy',
    type=Path,
    help='Path to directory with VCFs from Snippy'
)
@click.option(
    '--dir_medaka',
    type=Path,
    help='Path to directory with VCFs from Medaka (same name as Snippy)'
)
@click.option(
    '--vcf_snippy',
    type=Path,
    help='Path to Snippy VCF'
)
@click.option(
    '--vcf_medaka',
    type=Path,
    help='Path to Medaka VCF'
)
@click.option(
    '--min_quality',
    type=int,
    default=30,
    help='Genotype quality to filter Medaka calls by [30]'
)
@click.option(
    '--prefix',
    type=str,
    default='diff',
    help='Prefix for table output [diff]'
)
@click.option(
    '--complex',
    is_flag=True,
    help='Break complex variants from Snippy (MNP, COMPLEX) into SNPs'
)
def snp_difference(
    dir_snippy, dir_medaka, vcf_snippy, vcf_medaka, min_quality, complex, prefix
):

    """ Compare SNPs called with Snippy vs. Medaka  """

    if vcf_snippy and vcf_medaka:
        comparisons = [(vcf_snippy, vcf_medaka)]
    elif dir_snippy and dir_medaka:

        if not dir_snippy.exists() or not dir_medaka.exists():
            raise ValueError('Could not find directories.')

        snippy_file_names = [
            f.name for f in sorted(dir_snippy.glob('*.vcf'))
        ]
        medaka_file_names = [
            f.name for f in sorted(dir_medaka.glob('*.vcf'))
        ]
        comparisons = [
            (dir_snippy / name, dir_medaka / name)
            for name in snippy_file_names if name in medaka_file_names
        ]
        print(
            f'Detected {len(comparisons)} matching VCFs '
            f'in {dir_snippy} and {dir_medaka}'
        )
    else:
        raise ValueError(
            'Specify either single VCFs or directories '
            'containing VCFs for comparison.'
        )

    results = []
    for comparison in comparisons:
        # Snippy is the reference to compare to

        snippy_vcf = comparison[0]
        medaka_vcf = comparison[1]

        snippy = SnippySample(vcf=snippy_vcf)

        # SNP only called in Medaka v.1.0.1: medaka consensus --haploid
        medaka = MedakaSample(vcf=medaka_vcf)
        medaka.filter(min_quality=min_quality)

        # Consider only SNPs not COMPLEX or MNP

        # Might result in slightly overestimating false positives as SNPs may
        # be hidden in the complex variant types from Snippy - these are
        # actually broken in the core genome computation (snp-sites) and
        # should therefore be included

        if complex:
            snippy.data = break_complex(snippy.data)
        else:
            snippy.data = snippy.data[snippy.data.snp == True]

        common = snippy.data.merge(
            medaka.data, on=['chromosome', 'position'], how='inner'
        )

        true_positives = len(common)
        false_positives = len(medaka.data) - len(common)
        false_negatives = len(snippy.data) - len(common)
        true_negatives = 0  # always 0

        accuracy = (true_positives+true_negatives) / \
            (true_positives+true_negatives+false_positives+false_negatives)

        precision = true_positives/(true_positives+false_positives)
        recall = true_positives/(true_positives+false_negatives)

        f1 = 2*(recall*precision)/(recall+precision)

        results.append({
            'name': snippy.name,
            'snippy': len(snippy.data),
            'medaka': len(medaka.data),
            'true_positives': len(common),
            'false_positives': false_positives,
            'false_negatives': false_negatives,
            'true_negatives': true_negatives,
            'accuracy': accuracy,
            'precision': precision,
            'recall': recall,
            'f1': f1
        })

    data = pandas.DataFrame(results).set_index('name')

    print(data)

    summary = pandas.DataFrame(
        [data.mean(axis=0), data.std(axis=0), data.sem(axis=0)],
        index=['mean', 'sd', 'sem']
    )

    data.to_csv(f'{prefix}.data.tsv', sep='\t', index=True)
    summary.to_csv(f'{prefix}.summary.tsv', sep='\t', index=True)

    print(summary)


def break_complex(snippy_data: pandas.DataFrame):

    """ Break complex variants called with Snippy

    Break complex variants (MNP, COMPLEX) called with Snippy
    into single nucleotide polymorphisms. Medaka does not have
    the capacity to call these types, but may call some SNPs that
    are hidden in them. In the core genome computations, these
    variants from Snippy are eventually extracted as SNPs into
    the final alignment, so it woudl be prudent to include them
    for comparison, otherwise the false positives will be
    slightly overestimated.

    """

    complex_variants = snippy_data[snippy_data.snp == False]

    print(complex_variants)
    broken_rows = []
    for _, row in complex_variants.iterrows():
        ref = list(row['ref'])
        call = list(row['call'])

        # Reference and called allele are
        # always the same length
        for i, base in enumerate(ref):
            if base != call[i]:
                # SNP detected
                new_row = row.copy()
                new_row['position'] = new_row['position']+i
                new_row['ref'] = base
                new_row['call'] = call[i]
                new_row['alt'] = call[i]
                new_row['snp'] = True
                broken_rows.append(new_row)

    complex_broken = pandas.DataFrame(broken_rows)
    noncomplex_variants = snippy_data[snippy_data.snp == True]

    return pandas.concat([noncomplex_variants, complex_broken])\
        .sort_values(['chromosome', 'position'])
