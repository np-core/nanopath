import click
import pandas
import logging
import pyfastx

from pathlib import Path
from pysam import VariantFile
from nanopath.utils import PoreLogger


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
    help='Minimum read depth at genotype site to include [0]'
)
@click.option(
    '--min_likelihood',
    type=int,
    default=30,
    help='Minimum phred-scaled likelihood of genotype to include (0 - 999) [0]'
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
@click.option(
    '--megalodon_prefix',
    type=str,
    default=None,
    help='Additional prefix for names of Megalodon samples [none]'
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
    megalodon_prefix
):

    """ Merge variants detected by Megalodon with variants from Snippy """

    # TODO - check if there is discrepancy between Snippy VCF and Alignments
    # TODO - if so, use the VCF to generate the merged output alignments

    sc = SnippyCore(vcf=vcf_snippy, outdir=outdir, prefix=prefix)

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
                missing='N',
                min_depth=min_depth,
                min_likelihood=min_likelihood,
                name_prefix=megalodon_prefix
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
            merge=merge,
            reference=reference,
            fasta=outdir / f'{prefix+"." if prefix else ""}full.fasta')


class SnippyCore(PoreLogger):

    def __init__(
        self,
        vcf: Path,
        outdir: Path,
        prefix: str
    ):

        PoreLogger.__init__(self, level=logging.INFO, name='SnippyCore')

        self.vcf = vcf
        self.outdir = outdir
        self.prefix = prefix

        self.megalodon = []  # filtered megalodon samples
        self.megalodon_pass = []  # passing megalodon samples

        self.genotypes = dict()  # snippy genotype list for fasta output
        self.candidates = pandas.DataFrame()  # to check megalodon completeness
        self.total = 0  # total number of core variants
        self.total_samples = 0  # total number of samples

        self.outdir.mkdir(parents=True, exist_ok=True)
        self.parse_snippy_vcf()

    def parse_snippy_vcf(self):

        # TODO: Snippy VCF must be ordered - should always be the case!

        self.logger.info('Parsing Snippy VCF file, extracting core variants')

        calls = []
        self.genotypes = {}
        self.total = 0
        for rec in VariantFile(self.vcf).fetch():
            chromosome = rec.chrom
            position = int(rec.pos)
            reference = rec.ref
            variants = rec.alts

            for sample, data in rec.samples.items():
                try:
                    genotype = data['GT'][0]  # tuple, but only one entry
                except IndexError:
                    raise IndexError(
                        f'Error parsing GT: {chromosome} - {position}'
                    )

                if genotype > 0:
                    # Variant called
                    try:
                        call = variants[genotype - 1]
                    except IndexError:
                        raise IndexError(
                            f'Error parsing ALT: {chromosome} - {position}'
                        )
                else:
                    call = reference

                if sample not in self.genotypes.keys():
                    self.genotypes[sample] = [call]
                else:
                    self.genotypes[sample].append(call)

            calls.append(dict(
                chromosome=chromosome,
                position=position,
                ref=reference
            ))
            self.total += 1
            self.total_samples = len(self.genotypes)

        self.logger.info(
            f'Parsed {self.total} variants in '
            f'{self.total_samples} samples from Snippy'
        )

        # Check that all genotypes have total length
        # matching number of core variants:
        for sample, genotype in self.genotypes.items():
            if len(genotype) != self.total:
                self.logger.error(
                    f"Genotype for sample {sample} "
                    f"should be length: {self.total} ("
                    f"currently: {len(genotype)}) "
                )

        # Should be sorted, but make sure it is
        self.candidates = pandas.DataFrame(calls)\
            .sort_values(['chromosome', 'position'])

    def parse_megalodon(
        self,
        vcf: Path,
        missing: str = 'N',
        min_depth: int = 5,
        min_likelihood: int = 30,
        name_prefix: str = None,
    ):
        mv = MegalodonSample(
            vcf=vcf,
            missing=missing,
            name_prefix=name_prefix,
            outdir=self.outdir / 'variant_calls'
        )
        mv.filter(min_depth=min_depth, min_likelihood=min_likelihood)

        self.megalodon.append(mv)

    def check_megalodon(self, max_missingness: float = 0.2):

        self.megalodon_pass = []
        for sample in self.megalodon:
            # First filter out > missingess as it can only increase
            if sample.missingness > max_missingness:
                self.logger.info(
                    f'Removed {sample.name}, missing rate '
                    f'{round(sample.missingness, 4)} > {max_missingness}'
                )
                continue
            # Compare to core SNPs and check if all candidates present
            # Should never be larger than candidate variants input
            if sample.total < self.total:
                self.logger.info(
                    f'Missing {self.total - sample.total} '
                    f'candidate variants in {sample.name} - filling and '
                    f'recomputing missing rate'
                )
                # Add missing candidates (N) and update missingness
                sample.fill_candidates(variants=self.candidates)

                # Filter again for missingness in case > threshold
                if sample.missingness > max_missingness:
                    self.logger.info(
                        f'Removed {sample.name}, missing rate after filling '
                        f'{sample.missingness} > {max_missingness}'
                    )
                    continue

            self.megalodon_pass.append(sample)

        self.logger.info(
            f'Retained {len(self.megalodon_pass)} samples from Megalodon'
        )

    def create_alignment(self, fasta: Path, merge: bool = True):

        """ Create alignment file from Megalodon calls

        :param fasta: output alignment in fasta format
        :param merge: merge output alignment with fasta
        :return:

        """

        self.logger.info(f'Creating alignment from Megalodon calls: {fasta}')

        with fasta.open('w') as fasta_out:
            for sample in self.megalodon_pass:
                fasta_out.write(
                    sample.get_fasta_entry()
                )

            if merge:
                self.logger.info(
                    f'Merging Snippy core variants with Megalodon calls: {fasta}'
                )
                for sample, genotype in self.genotypes.items():
                    record = f">{sample}\n{''.join(genotype)}\n"
                    fasta_out.write(record)

    def filter_invariant_sites(self, alignment: Path, fasta: Path):

        # Read alignment into a DataFrame and check columns for invariant sites
        # This is a bit weird, but ok for small (core) alignments

        self.logger.info(f'Removing invariant sites from: {alignment}')

        samples = []
        genotypes = []
        with alignment.open('r') as aln:
            for line in aln:
                if line.startswith('>'):
                    name = line.strip().replace('>', '')
                else:
                    genotypes.append(list(
                        line.strip()
                    ))
                    samples.append(name)

        df = pandas.DataFrame(data=genotypes, index=samples)

        invariant_sites = []
        # Full invariant sites without missing
        for site_index, col in df.iteritems():
            if col.nunique() == 1:
                self.logger.info(f'Found invariant site: {site_index}')
                invariant_sites.append(site_index)
            if col.nunique() == 2 and 'N' in col.unique():
                self.logger.info(f'Found invariant site: {site_index}')
                invariant_sites.append(site_index)

        self.logger.info(
            f'Removing {len(invariant_sites)} invariant sites from alignment'
        )

        df = df.drop(invariant_sites, axis=1)

        self.logger.info(f'Writing filtered alignment to: {fasta}')

        with fasta.open('w') as out:
            for sample, row in df.iterrows():
                record = f'>{sample}\n{"".join(row)}\n'
                out.write(record)

    def create_full_alignment(self, merge: bool, reference: Path, fasta: Path):

        """ Create a full alignment of the reference genome (substituted)

        :param merge:
        :param reference:
        :param fasta:
        :return:

        """

        self.logger.info(
            f'Creating full alignment from Megalodon calls: {fasta}'
        )
        ref = self.read_reference(fasta=reference)

        if merge:
            genotypes = self.genotypes.copy()
        else:
            genotypes = {}

        for sample in self.megalodon_pass:
            if sample.name in genotypes.keys():
                self.logger.error(
                    f'Megalodon sample {sample.name} is already in Snippy '
                    f'genotype list, this should not be the case.'
                )
                exit(1)
            genotypes[sample.name] = sample.get_genotype()

        # Insert calls into reference sequences and output to FASTA
        # Multiple chromosomes are concatenated into a single sequence
        with fasta.open('w') as alignment:
            for name, genotype in genotypes.items():
                self.logger.info(
                    f'Replacing variant sites in reference sequences for {name}'
                )
                variant_ref = dict()
                for chromosome in self.candidates.chromosome.unique():
                    try:
                        seq = list(ref[chromosome])
                    except KeyError:
                        self.logger.error(
                            f'Sequence {chromosome} not in reference.'
                        )
                        exit(1)

                    for i, row in self.candidates.iterrows():
                        position = row['position']
                        seq[position] = genotype[i]

                    variant_ref[chromosome] = "".join(seq)

                concat_seq = ''.join(
                    variant_ref[key] for key in sorted(variant_ref.keys())
                )
                record = f'>{name}\n{concat_seq}\n'
                alignment.write(record)

    @staticmethod
    def read_reference(fasta: Path) -> dict:

        return {
            name: seq for name, seq in
            pyfastx.Fasta(str(fasta), build_index=False)
        }


class MegalodonSample(PoreLogger):

    def __init__(
        self,
        vcf: Path,
        missing: str = 'N',
        outdir: Path = Path('calls'),
        name_prefix: str = None
    ):

        PoreLogger.__init__(self, level=logging.INFO, name='Megalodon')

        self.vcf = vcf
        self.missing = missing
        self.outdir = outdir

        self.data = pandas.DataFrame()
        self.data_filtered = pandas.DataFrame()

        # Stats for later use

        self.total = 0
        self.filtered = 0
        self.missingness = 0.

        # Use name of directory containing variant file
        self.name = vcf.absolute().parent.name

        if name_prefix:
            self.name = f'{name_prefix}_{self.name}'

        self.outdir.mkdir(parents=True, exist_ok=True)

        self.logger.info(f'Processing: {self.name}')
        self.parse()

    def parse(self):

        var = 0
        calls = []
        for rec in VariantFile(self.vcf).fetch():

            chromosome = rec.chrom
            position = int(rec.pos)
            reference = rec.ref
            variants = rec.alts  # tuple

            sample = rec.samples['SAMPLE']

            depth = int(sample['DP'])
            likelihood = int(sample['GQ'])  # haploid, not normalized

            try:
                genotype = sample['GT'][0]  # tuple, but only one entry
            except IndexError:
                raise IndexError(f'Error parsing GT: {chromosome} - {position}')

            if genotype > 0:
                # Variant called
                try:
                    variant = variants[genotype - 1]
                except IndexError:
                    raise IndexError(f'Error parsing ALT: {chromosome} - {position}')

                self.logger.debug(
                    f'Detected variant: {variant}/{reference} - {position}'
                    f' - {depth} - {likelihood}'
                )

                var += 1
                call = variant
                alt = variant
                is_variant = True
            else:
                call = reference
                alt = '-'
                is_variant = False

            calls.append(dict(
                position=position,
                call=call,
                ref=reference,
                alt=alt,
                chromosome=chromosome,
                read_depth=depth,
                likelihood=likelihood,
                is_variant=is_variant
            ))

        self.data = pandas.DataFrame(calls).sort_values('position')

        self.logger.info(
            f'Detected {var} SNPs called by Megalodon'
        )

        self.data.to_csv(
            f'{self.outdir / f"{self.name}.calls.raw.tsv"}', sep='\t', index=False
        )
        self.data[self.data['is_variant'] == True].to_csv(
            f'{self.outdir / f"{self.name}.variants.raw.tsv"}', sep='\t', index=False
        )

    def filter(self, min_depth: int = 10, min_likelihood: int = 30):

        self.data_filtered = self.data.copy()

        self.data_filtered.loc[
            (self.data_filtered['read_depth'] < min_depth) |
            (self.data_filtered['likelihood'] < min_likelihood),
            'call'
        ] = self.missing

        filtered = self.data_filtered.loc[
            self.data_filtered['call'] == self.missing
        ]

        filtered_variants = filtered[filtered['is_variant'] == True]
        filtered_references = filtered[filtered['is_variant'] == False]

        self.logger.info(
            f'Filtered {len(filtered)} / {len(self.data_filtered)} calls --> '
            f'{len(filtered_variants)} variants, {len(filtered_references)} '
            f'reference'
        )

        self.filtered = len(filtered)
        self.total = len(self.data_filtered)
        self.missingness = self.filtered / self.total

        self.logger.info(
            f'Missing rate after filtering: {round(self.missingness, 4)}'
        )

        filtered.to_csv(
            f'{self.outdir / f"{self.name}.all.filtered.tsv"}', sep='\t', index=False
        )
        self.data_filtered.to_csv(
            f'{self.outdir / f"{self.name}.calls.filtered.tsv"}', sep='\t', index=False
        )
        filtered_variants.to_csv(
            f'{self.outdir / f"{self.name}.variants.filtered.tsv"}',
            sep='\t', index=False
        )

        self.data_filtered = self.data_filtered.sort_values(
            ['chromosome', 'position']
        )  # make sure it's sorted

    def fill_candidates(self, variants: pandas.DataFrame):

        """ Compare variants from Snippy to this sample and fill missing """

        if self.data_filtered.empty:
            raise ValueError('Megalodon variants must first be filtered.')

        # Works only if the dataframes do not contain duplicates themselves
        # which is the case here:

        missing = pandas.concat([self.data_filtered, variants])\
            .drop_duplicates(subset=['chromosome', 'position'], keep=False)

        missing['call'] = [self.missing for _ in range(len(missing))]

        # Introduces NA into DataFrame!
        self.data_filtered = pandas.concat([self.data_filtered, missing])\
            .sort_values(['chromosome', 'position']).reset_index(drop=True)

        # Update stats
        self.total = len(self.data_filtered)
        self.filtered = len(self.data_filtered.loc[
            self.data_filtered['call'] == self.missing
        ])
        self.missingness = self.filtered / self.total

    def get_fasta_entry(self) -> str:

        header = f">{self.name}"
        seq = "".join(self.data_filtered.call)

        return f'{header}\n{seq}\n'

    def get_genotype(self) -> list:

        return self.data_filtered.call.to_list()
