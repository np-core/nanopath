""" Modules to process program output in the pipelines """

import json
import pandas
import logging
import numpy as np
import seaborn as sns

from matplotlib import pyplot as plt
from pySankey import sankey

from pathlib import Path
from nanopath.utils import create_fastx_index, run_cmd, PoreLogger

import dendropy

import pyfastx

try:
    from pysam import VariantFile, VariantRecord
except ModuleNotFoundError:
    pass # windows


def read_fasta(fasta: Path) -> dict:
    return {
        name: seq.upper() for name, seq in
        pyfastx.Fasta(str(fasta), build_index=False)  # capital bases
    }


class MegalodonCore:

    """ Class to parse Snippy core SNPs and merge with Megalodon calls """

    def __init__(self, core_vcf: Path = None):

        logging.basicConfig(
            level=logging.INFO,
            format=f"[%(asctime)s] [MegalodonCore]     %(message)s",
            datefmt='%H:%M:%S'
        )

        self.logger = logging.getLogger()

        self.core_vcf = core_vcf

        self.genotypes = dict()  # snippy genotype list for fasta output
        self.candidates = pandas.DataFrame()  # to check megalodon completeness
        self.total = 0  # total number of core variants
        self.total_samples = 0  # total number of samples

        self.megalodon = []  # filtered megalodon samples
        self.megalodon_pass = []  # passing megalodon samples

        if core_vcf is not None:
            self.parse_snippy_core_vcf()

    def parse_megalodon(
        self,
        vcf: Path,
        min_depth: int = 5,
        min_quality: int = 30
    ):
        mv = MegalodonSample(vcf=vcf)
        mv.filter(min_depth=min_depth, min_quality=min_quality)

        self.megalodon.append(mv)

    def check_megalodon(self, max_missingness: float = 0.05):

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

    def parse_snippy_core_vcf(self):

        # TODO: Snippy VCF must be ordered - should always be the case!

        self.logger.info('Processing: Snippy VCF')

        calls = []
        self.genotypes = {}
        self.total = 0
        for rec in VariantFile(self.core_vcf).fetch():
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
            f'Parsed {self.total} core variants across '
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

    def create_alignment(self, fasta: Path, merge: bool = True):

        """ Create alignment file from Megalodon calls

        :param fasta: output alignment in fasta format
        :param merge: merge output alignment with fasta
        :return:

        """

        self.logger.info(f'Creating alignment from calls: {fasta}')

        with fasta.open('w') as fasta_out:
            for sample in self.megalodon_pass:
                fasta_out.write(
                    sample.get_fasta_entry()
                )

            if merge:
                self.logger.info(
                    f'Merging Snippy variants to file: {fasta}'
                )
                for sample, genotype in self.genotypes.items():
                    record = f">{sample}\n{''.join(genotype)}\n"
                    fasta_out.write(record)

    def filter_invariant_sites(self, alignment: Path, fasta: Path):

        # Read alignment into a DataFrame and check columns for invariant sites

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

        if not genotypes:
            raise ValueError('Could not read alignment')

        df = pandas.DataFrame(data=genotypes, index=samples)

        self.logger.info(
            f'Input alignment length: {len(genotypes[0])}'
        )

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

        self.logger.info(
            f'Filtered alignment length: {len(df.columns)}'
        )
        self.logger.info(f'Writing filtered alignment to: {fasta}')

        with fasta.open('w') as out:
            for sample, row in df.iterrows():
                record = f'>{sample}\n{"".join(row)}\n'
                out.write(record)

    def create_full_alignment(self, reference: Path, fasta: Path):

        """ Create a full alignment of the reference genome (substituted)

        :param reference:
        :return:

        """

        self.logger.info(
            f'Creating full alignment: {fasta}'
        )

        ref = read_fasta(fasta=reference)
        genotypes = self.genotypes.copy()

        for sample in self.megalodon_pass:
            if sample.name in genotypes.keys():
                self.logger.error(
                    f'ONT sample {sample.name} is already in Snippy '
                    f'genotype list, this should not be the case.'
                )
                exit(1)
            genotypes[sample.name] = sample.get_genotypes()

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


class MedakaCore:

    def __init__(
        self,
        outdir: Path = None,
        reference: Path = None,
        prefix: str = 'core'
    ):

        logging.basicConfig(
            level=logging.INFO,
            format=f"[%(asctime)s] [MedakaCore]     %(message)s",
            datefmt='%H:%M:%S'
        )

        self.logger = logging.getLogger()

        self.reference = reference

        self.outdir = outdir
        self.prefix = prefix

        self.medaka = []  # medaka samples
        self.snippy = []  # snippy samples

    def parse_medaka(
        self,
        medaka_directory: Path,
        min_quality: int = 30,
        min_depth: int = 10
    ):
        """ Parse VCF file and corresponding aligned reads from Medaka

        :param medaka_directory: path to directory containing .vcf / .bam
        :param min_quality: minimum genotype quality to include from Medaka calls
        :param min_depth: minimum depth from read alignment against reference
            to consider a site present for core computation

        """

        # Process Medaka VCF and BAM file to define
        # excluded (non-core) positions
        for vcf in sorted([
            f for f in medaka_directory.glob('*.vcf')
        ]):
            bam_file = vcf.parent / f'{vcf.stem}.bam'
            ms = MedakaSample(
                vcf=vcf,
                bam=bam_file,
                min_depth=min_depth
            )
            ms.filter(min_quality=min_quality)

            self.logger.info(
                f'Processed {len(ms.data)} variants in sample: {ms.name}'
            )
            self.medaka.append(ms)

    def parse_snippy(self, snippy_directory: Path):

        """ Parse VCF files and corresponding aligned reference sequence from Snippy

        :param snippy_directory: path to directory containing .vcf / .aligned.fa

        """

        for vcf in sorted([
            f for f in snippy_directory.glob('*.vcf')
        ]):
            # Process Snippy VCF + aligned reference genome
            # to define excluded (non-core) positions
            alignment_fasta = vcf.parent / f'{vcf.stem}.aligned.fa'
            ss = SnippySample(vcf=vcf, aligned=alignment_fasta)

            self.logger.info(
                f'Processed {len(ss.data)} variants in sample: {ss.name}'
            )
            self.snippy.append(ss)

    def call_hybrid_core(
        self,
        include_reference: bool = True
    ):

        """ Determine core variants from Snippy and Medaka samples  """

        # Merge samples from Snippy and Medaka

        samples = self.snippy + self.medaka

        exclude = {}
        snp_positions = {}
        for sample in samples:
            # For each  chromosome add the excluded positions to a
            # dictionary of chromosome: excluded sites for all samples
            for chrom, excluded_positions in sample.excluded_positions.items():
                if chrom not in exclude.keys():
                    exclude[chrom] = excluded_positions
                else:
                    exclude[chrom] += excluded_positions

            # For each chromosome, add the called positions to a dictionary
            # of chromosome: called positions for all samples
            for chrom, positions in sample.get_snp_positions().items():
                if chrom in snp_positions.keys():
                    snp_positions[chrom] += positions
                else:
                    snp_positions[chrom] = positions

        # For each chromosomes, determine the sites to exclude across all samples
        # corresponding to called sites that fall into low coverage or gaps
        # in any sample (non core sites)
        snps_to_exclude = {}
        for chrom, excluded_sites in exclude.items():
            snps_to_exclude[chrom] = \
                set(snp_positions[chrom]).intersection(
                    pandas.Series(excluded_sites).unique().tolist()
                )

        # For each chromosome, determine the unique sites across all samples
        # and check how many of them fall into non-core sites, then
        # remove from unique sites to get the final core sites:
        core_sites = {}
        for chrom, to_exclude in snps_to_exclude.items():
            unique_snp_positions = pandas.Series(
                snp_positions[chrom]
            ).unique().tolist()

            self.logger.info(
                f'Detected a total of {len(unique_snp_positions)} '
                f'unique variants on {chrom}'
            )
            self.logger.info(
                f'Excluded {len(to_exclude)} variants due to '
                f'low coverage or gaps on {chrom}'
            )
            core_sites[chrom] = [
                p for p in unique_snp_positions if p not in to_exclude
            ]
            self.logger.info(
                f'Keeping {len(core_sites[chrom])} core variants on {chrom}'
            )

        # Create the filtered core genome sites for storage
        # in filtered data attribute, and output of sites to VCF
        all_core_data = []
        self.logger.info('Processing core variant calls')
        for sample in samples:
            filtered = []
            for chrom, data in sample.data.groupby('chromosome'):
                filtered.append(
                    data[data['position'].isin(core_sites[chrom])]
                )
            sample.data_core = pandas.concat(filtered)
            all_core_data.append(sample.data_core)

        # Output #

        # Get a table of unique core sites and their associated data:
        self.logger.info(
            f'Writing full core variant table to: {self.prefix}.tsv'
        )
        core_site_data = pandas.concat(all_core_data).drop_duplicates(
            ignore_index=True, subset=['chromosome', 'position']
        )

        core_site_table = core_site_data[
            ['chromosome', 'position', 'ref', 'alt']
        ].sort_values(['chromosome', 'position'])

        core_site_table.to_csv(f'{self.prefix}.tsv', sep='\t', index=False)

        # Read the reference sequence used for alignment
        ref = read_fasta(self.reference)

        # For each sample insert the called variant / reference allele
        # into the reference sequence, concatenate chromosomes and write
        # to file:

        self.logger.info(
            f'Writing full core variant alignment to: {self.prefix}.full.aln'
        )

        with open(f'{self.prefix}.full.aln', 'w') as full_alignment:
            for sample in samples:
                seq_concat = sample.replace_variants(reference=ref)
                full_alignment.write(f'>{sample.name}\n{seq_concat}\n')

            if include_reference:
                ref_seq = ''.join(ref[key]for key in sorted(ref.keys()))
                full_alignment.write(f'>Reference\n{ref_seq}\n')

        # Use snp-sites to extract the final core variant site alignment
        self.logger.info(
            f'Writing single nucleotide polymorphism alignment to: {self.prefix}.aln'
        )
        run_cmd(
            f'snp-sites -c -o {self.prefix}.aln {self.prefix}.full.aln'
        )

        self.logger.info(
            f'Writing single nucleotide polymorphism data to: {self.prefix}.vcf'
        )
        run_cmd(
            f'snp-sites -c -v -o {self.prefix}.vcf {self.prefix}.full.aln'
        )

        snp_count = sum(
            [1 for _ in VariantFile(f'{self.prefix}.vcf').fetch()]
        )
        self.logger.info(
            f'Final single nucleotide polymorphisms in alignment: {snp_count}'
        )


class Sample(PoreLogger):

    """ Base class for Snippy, Medaka and Megalodon per sample variant calls """

    def __init__(self, vcf: Path, stats: Path = None):

        self.vcf = vcf

        if stats:
            self.stats = self.read_pysamstats(file=stats)
        else:
            self.stats = stats

        self.features = pandas.DataFrame()  # snp random forest classifier feature storage [train, eval]
        self.filtered = pandas.DataFrame()  # filtered features where prediction SNP == True [eval]

        self.missing = 'N'

        self.data = pandas.DataFrame()  # all sites
        self.data_core = pandas.DataFrame()  # core sites only

        # Stats
        self.total = 0
        self.filtered = 0
        self.missingness = 0.

        self.name = vcf.stem

        logging.basicConfig(
            level=logging.INFO,
            format=f"[%(asctime)s] [{self.name}]     %(message)s",
            datefmt='%H:%M:%S'
        )

        PoreLogger.__init__(self)

    @staticmethod
    def read_pysamstats(file):

        return pandas.read_csv(file, sep="\t")

    def get_snp_positions(self) -> dict:

        """ Get SNP positions as dictionary of contigs - lists """

        return {
            chrom: data.position.tolist()
            for chrom, data in self.data.groupby('chromosome')
        }

    def replace_variants(self, reference: dict) -> str:

        """ Insert filtered (core site) variants into the reference sequence

        :param reference: str chromosome: str sequence
        :return: concatenated and variant-filled reference chromosome sequences
        """

        variant_reference_sequences = {}
        for chromosome, data in self.data_core.groupby('chromosome'):
            self.logger.debug(
                f'Inserting core variant calls into: {self.name} - {chromosome}'
            )

            try:
                seq = list(reference[chromosome])
            except KeyError:
                self.logger.error(
                    f'Sequence {chromosome} not in reference'
                )
                raise

            before = len(seq)

            for i, row in data.iterrows():
                position = row['position']
                call = list(row['call'])

                # Accounting for MNP and COMPLEX
                # 1-based position to sequence list index
                start_idx, end_idx = position-1, position-1+len(call)
                self.logger.debug(
                    f'Replacing: {position}, {seq[start_idx:end_idx]}, '
                    f'{call}, {start_idx}, {end_idx}'
                )
                seq[start_idx:end_idx] = call

            after = len(seq)

            # Make sure nothing dodgy happens during inserts
            try:
                assert before == after
            except AssertionError:
                self.logger.error(
                    'Something super dodgy happened during reference insertion'
                )
                raise

            variant_reference_sequences[chromosome] = "".join(seq)

        # Concatenate chromosomes and return as single sequence for snp-sites
        return ''.join(
            variant_reference_sequences[key]
            for key in sorted(variant_reference_sequences.keys())
        )

    def filter(self, min_quality: int = 30):

        # For Illumina, not ONT, I think

        self.data.loc[
            (self.data['quality'] < min_quality), 'call'
        ] = self.missing

        filtered = self.data.loc[
            self.data['call'] == self.missing
        ]

        self.filtered = len(filtered)
        self.total = len(self.data)
        self.missingness = self.filtered / self.total

        self.logger.info(
            f'Filtered {len(filtered)}/{len(self.data)} < Q{min_quality} '
            f'({round(self.missingness*100, 4)}%) - '
            f'keeping: {len(self.data) - len(filtered)}'
        )

        self.data = self.data[
            self.data['call'] != self.missing
        ]

        self.data = self.data.sort_values(
            ['chromosome', 'position']
        )  # make sure it's sorted

    def write_vcf(self, out: Path):

        vi = VariantFile(self.vcf)
        vo = VariantFile(str(out), 'w', header=vi.header)

        records_out = [
            (row['chromosome'], row['position']) for i, row in self.data.iterrows()
        ]

        written = 0
        for rec in vi.fetch():
            if (rec.chrom, int(rec.pos)) in records_out:  # TODO: does not consider Snippy broken COMPLEX
                written += 1
                vo.write(rec)

        self.logger.info(f'Wrote {written}/{len(records_out)} filtered records to file: {out}')


class ForestSample(Sample):

    """ Basic SNP samples only for parsing after Random Forest filter """

    def __init__(self, vcf: Path):

        Sample.__init__(self, vcf=vcf)

        self.parse()

    def parse(self):

        var = 0
        total = 0
        calls = []
        for rec in VariantFile(self.vcf).fetch():

            total += 1

            chromosome = rec.chrom
            position = int(rec.pos)
            reference = rec.ref
            variants = rec.alts  # tuple
            qual = float(rec.qual)

            # Snippy output VCFs only have variant sites: GT = 1/1
            try:
                call = variants[0]
            except IndexError:
                self.logger.debug(
                    f'Could not detect alts for variant in '
                    f'file {self.vcf} - position {position}'
                )
                raise

            # Only applies to SNP
            if len(call) == 1 and len(reference) == 1:
                calls.append(dict(
                    chromosome=chromosome,
                    position=position,
                    call=call,
                    ref=reference,
                    alt=call,
                    quality=qual,
                    snp=True
                ))

                var += 1

        self.data = pandas.DataFrame(calls).sort_values(
            ['chromosome', 'position']
        )


class SnippySample(Sample):

    """ Snippy samples are called de novo from Illumina """

    def __init__(
        self,
        vcf: Path,
        aligned: Path = None,
        break_complex: bool = True
    ):

        Sample.__init__(self, vcf=vcf)

        self.aligned = aligned
        if self.aligned:
            self.excluded_positions: dict = self.get_excluded_positions()

        self.logger.debug(f'Processing Snippy variant calls: {self.vcf}')
        self.parse()

        if break_complex:
            self.logger.debug(f'Breaking COMPLEX variants (excludes INDEL): {self.vcf}')
            self.data = self.break_complex()
        else:
            self.logger.debug(
                f'Not breaking COMPLEX variants (excludes COMPLEX, MNP, INDEL) '
                f'retaining SNPs only: {self.vcf}'
            )
            self.data = self.data[self.data['snp'] == True]

    def break_complex(self):

        """ Break complex variants called with Snippy

        Break complex variants (MNP, COMPLEX) called with Snippy
        into single nucleotide polymorphisms. Medaka does not have
        the capacity to call these types, but may call some SNPs that
        are hidden in them. In the core genome computations, these
        variants from Snippy are eventually extracted as SNPs into
        the final alignment, so it would be prudent to include them
        for comparison, otherwise false positives will be
        slightly overestimated.

        """

        # Consider only SNPs not COMPLEX or MNP

        # Might result in slightly overestimating false positives as SNPs may
        # be hidden in the complex variant types from Snippy - these are
        # actually broken in the core genome computation (snp-sites) and
        # should therefore be included

        complex_variants = self.data[self.data.snp == False]

        broken_rows = []
        for _, row in complex_variants.iterrows():
            ref = list(row['ref'])
            call = list(row['call'])

            # Reference and called allele are
            # always the same length from parsing
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
        noncomplex_variants = self.data[self.data.snp == True]

        return pandas.concat([noncomplex_variants, complex_broken])\
            .sort_values(['chromosome', 'position'])

    def parse(self):

        var = 0
        total = 0
        calls = []
        for rec in VariantFile(self.vcf).fetch():

            total += 1

            chromosome = rec.chrom
            position = int(rec.pos)
            reference = rec.ref
            variants = rec.alts  # tuple
            qual = float(rec.qual)

            info = rec.info
            print(info)
            if 'snp' in info['TYPE'] \
                or 'mnp' in info['TYPE'] \
                    or 'complex' in info['TYPE']:

                # Snippy output VCFs only have variant sites: GT = 1/1
                try:
                    call = variants[0]
                except IndexError:
                    self.logger.debug(
                        f'Could not detect alts for variant in '
                        f'file {self.vcf} - position {position}'
                    )
                    raise

                # Check if COMPLEX type reference allele length is
                # the same length as called variant allele, if not
                # it should not be included (INDEL)

                # Only applies to COMPLEX, not SNP or MNP
                if len(call) != len(reference):
                    continue

                calls.append(dict(
                    chromosome=chromosome,
                    position=position,
                    call=call,
                    ref=reference,
                    alt=call,
                    quality=qual,
                    snp=True if 'snp' in info['TYPE'] else False
                ))

                var += 1

        self.data = pandas.DataFrame(calls).sort_values(
            ['chromosome', 'position']
        )

    def get_excluded_positions(self):

        """ Get excluded positions based on aligned reference sequence """

        aligned_reference = read_fasta(self.aligned)

        excluded = {}
        for chrom, seq in aligned_reference.items():
            # List of positions containing gaps (-) or
            # regions of low coverage (N) determined by Snippy
            excluded[chrom] = [
                i for i, s in enumerate(seq, 1) if s == 'N' or s == '-'  # TODO: fix here not in ATCG
            ]  # 1-based indexing of sites to match VCF

        return excluded


class ClairSample(Sample):

    """ Clair samples are called de novo from ONT"""

    def __init__(
        self, vcf: Path, stats: Path = None
    ):
        Sample.__init__(self, vcf=vcf, stats=stats)

        self.parse()

    def parse(self):

        var = 0
        total = 0
        calls = []
        for rec in VariantFile(self.vcf).fetch():

            total += 1
            chromosome = rec.chrom
            position = int(rec.pos)
            reference = rec.ref
            variants = rec.alts  # tuple
            qual = float(rec.qual)

            # Clair haploid mode variants
            try:
                call = variants[0]
            except IndexError:
                self.logger.debug(
                    f'Could not detect alts for variant in '
                    f'file {self.vcf} - position {position}'
                )
                raise

            if len(reference) == 1 and len(call) == 1:
                # SNPs only
                calls.append(dict(
                    position=position,
                    call=call,
                    ref=reference,
                    alt=call,
                    chromosome=chromosome,
                    quality=qual,
                    snp=True
                ))

                var += 1
            else:
                # NO INDELS
                continue

        self.data = pandas.DataFrame(calls).sort_values(
            ['chromosome', 'position']
        )


class MedakaSample(Sample):

    """ Medaka samples are called de novo from ONT """

    def __init__(
        self,
        vcf: Path,
        stats: Path = None,
        bam: Path = None,
        min_depth: int = 10,
    ):

        Sample.__init__(self, vcf=vcf, stats=stats)

        self.bam = bam
        self.min_depth = min_depth

        # Filter reference sequence by depth
        if bam:
            self.excluded_positions: dict = self.get_excluded_positions()

        self.logger.debug(f'Processing Medaka variant calls: {self.name}')
        self.parse()

    def parse(self):

        var = 0
        total = 0
        calls = []
        for rec in VariantFile(self.vcf).fetch():

            total += 1
            chromosome = rec.chrom
            position = int(rec.pos)
            reference = rec.ref
            variants = rec.alts  # tuple
            qual = float(rec.qual)

            # Medaka output VCFs only have variant sites: GT = 1/1
            try:
                call = variants[0]
            except IndexError:
                self.logger.debug(
                    f'Could not detect alts for variant in '
                    f'file {self.vcf} - position {position}'
                )
                raise

            calls.append(dict(
                position=position,
                call=call,
                ref=reference,
                alt=call,
                chromosome=chromosome,
                quality=qual
            ))

            var += 1

        self.data = pandas.DataFrame(calls).sort_values(
            ['chromosome', 'position']
        )

    def get_excluded_positions(self) -> dict:

        """ Assess missing and excluded sites from Medaka read alignment

        Analogous to the {id}.aligned.fa created by Snippy
        where missing sites (depth = 0) are gapped (-) and
        sites with low coverage are replaced with missing (N)

        Dependency: samtools depth

        returns: dictionary with str chrom: list of excluded sites
        """

        self.logger.info(
            f'Core site assessment for sample: {self.name}'
        )

        run_cmd(
            f'samtools depth -a {self.bam} > {self.name}.depth.txt', shell=True
        )  # Include all sites with zero read depth

        depth = pandas.read_csv(
            f'{self.name}.depth.txt', sep='\t', header=None,
            names=['chromosome', 'position', 'depth']
        )

        self.logger.debug(
            f'Determining excluded positions for sample: {self.name}'
        )

        excluded = {}
        for chrom, data in depth.groupby('chromosome'):

            # Get positions where depth is 0 or < min_depth
            unsupported = data[
                (data['depth'] == 0) | (data['depth'] < self.min_depth)
            ]
            excluded[chrom] = unsupported.position.tolist()
            self.logger.debug(
                f'Found {len(unsupported)} unsupported positions on {chrom}'
            )

        Path(f'{self.name}.depth.txt').unlink()

        return excluded


class MegalodonSample(Sample):

    """ Megalodon samples are called from candidate SNPs detected with Snippy """

    def __init__(self, vcf: Path):
        Sample.__init__(self, vcf=vcf)

        self.logger.info(f'Processing Megalodon variant calls:: {self.name}')
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
            quality = int(sample['GQ'])  # haploid, not normalized

            # Megalodon output VCFs can have reference or variant sites
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
                    f' - {depth} - {quality}'
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
                quality=quality,
                is_variant=is_variant
            ))

        self.data = pandas.DataFrame(calls).sort_values('position')

        self.logger.info(
            f'Detected {var} SNPs called by Megalodon'
        )

    def filter(self, min_depth: int = 10, min_quality: int = 30):

        self.data_core = self.data.copy()

        self.data_core.loc[
            (self.data_core['read_depth'] < min_depth) |
            (self.data_core['quality'] < min_quality),
            'call'
        ] = self.missing

        filtered = self.data_core.loc[
            self.data_core['call'] == self.missing
            ]

        filtered_variants = filtered[filtered['is_variant'] == True]
        filtered_references = filtered[filtered['is_variant'] == False]

        self.logger.info(
            f'Filtered {len(filtered)} / {len(self.data_core)} calls --> '
            f'{len(filtered_variants)} variants, {len(filtered_references)} '
            f'reference'
        )

        self.filtered = len(filtered)
        self.total = len(self.data_core)
        self.missingness = self.filtered / self.total

        self.logger.info(
            f'Missing rate after filtering: {round(self.missingness, 4)}'
        )

        self.data_core = self.data_core.sort_values(
            ['chromosome', 'position']
        )  # make sure it's sorted

    def fill_candidates(self, variants: pandas.DataFrame):

        """ Compare core variants from Snippy to this sample and fill missing """

        if self.data_core.empty:
            raise ValueError('Megalodon variants must first be filtered')

        # Works only if the dataframes do not contain duplicates themselves
        # which is always the case here:

        missing = pandas.concat([self.data_core, variants])\
            .drop_duplicates(subset=['chromosome', 'position'], keep=False)

        missing['call'] = [self.missing for _ in range(len(missing))]

        # Introduces NA into DataFrame!
        self.data_core = pandas.concat([self.data_core, missing])\
            .sort_values(['chromosome', 'position']).reset_index(drop=True)

        # Update stats
        self.total = len(self.data_core)
        self.filtered = len(
            self.data_core.loc[self.data_core['call'] == self.missing]
        )
        self.missingness = self.filtered / self.total

    def get_fasta_entry(self) -> str:

        header = f">{self.name}"
        seq = "".join(self.data_core.call)

        return f'{header}\n{seq}\n'

    def get_genotype(self) -> list:

        return self.data_core.call.to_list()


class SnippyAnnotator:

    """ Annotate SNPs using a reference genome and VCF from Snippy """

    def __init__(self, reference: Path, vcf: Path = None, gubbins: bool = False):

        """

        :param reference: reference genome with coding region annotations (GBK)
        :param vcf: VCF file from Snippy-Core or Gubbins
        :param gubbins: indicate whether VCF file is from Gubbins
            If true, removes all variants that contain any other site than

        """

        self.reference = reference
        self.vcf = vcf

    def find_branching(self, nexus: Path, nodes: str or list):

        """ Annotate branching SNPs from Treetime

        :param nexus: path to tree file with mutation annotations from Treetime
        :param nodes: node or list of nodes to extract mutations to annotate

        """

        if isinstance(nodes, str):
            nodes = [nodes]

        tree = dendropy.Tree.get(
            file=nexus.open(), schema="nexus", preserve_underscores=True
        )

        for n in tree.preorder_node_iter():
            if n.label in nodes:
                s = ""
                for a in n.annotations:
                    # Timetree has some weirdannotation parsing,
                    # so merge into string and extract mutations
                    s += a.name + a.value + ","

                mutations = s[s.find('"')+1:s.rfind('"')].split(',')
                print(mutations)


class AssemblyProcessor(PoreLogger):

    """ Process results from metagenome assembly """

    def __init__(
          self,
          fasta: Path,
          assembler: str = 'flye',
          info_file: Path = None,
          verbose: bool = True
    ):

        PoreLogger.__init__(
            self, level=logging.INFO if verbose else logging.ERROR
        )

        self.fasta = fasta
        self.assembler = assembler

        self.fxi = None
        self.fxi_file = None

        self.info_file = info_file

        # Some basic filters to extract potential candidates from Flye

        self.min_bacteria = 500000  # minimum length to consider chromosome
        self.max_virus = 100000  # max length to consider linear virus
        self.min_virus_mult = 3  # minimum multiplicity to look for viruses

    def get_server_data(self) -> dict:

        return self.get_assembly_data()

    def get_assembly_data(self) -> dict:

        """ What do we need to know about the assembly?

        Gets a bunch of data from the assembly output directory
        and processes it according to selected assembler - currently Flye.

        """

        if self.assembler == 'flye':
            self.logger.info('Read assembly information file from Flye')

            try:
                info = self.read_assembly_info()
            except FileNotFoundError:
                return {}

            # Simple check for chromosomes / plasmids /viruses

            bacterial_chromosomes = info.loc[
                (info['length'] > self.min_bacteria) &
                (info['circular'] == '+'), :
            ]

            # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4379921/

            plasmids = info.loc[
                (info['length'] < self.min_bacteria) &
                (info['circular'] == '+'), :
            ]

            viruses = info.loc[
                (info['length'] < self.max_virus) &
                (info['circular'] == '-') &
                (info['multiplicity'] > self.min_virus_mult), :
            ]

            return {
                'contigs': self.data_to_report(data=info),
                'candidates': {
                    'bacteria': self.data_to_report(
                        data=bacterial_chromosomes
                    ),
                    'viruses': self.data_to_report(data=viruses),
                    'plasmids': self.data_to_report(data=plasmids)
                }
            }

    @staticmethod
    def data_to_report(data: pandas.DataFrame):

        """ Transform contig data to JSON for Vue application """

        return data.to_dict(orient='records')

    def read_fasta(self):

        """ Read the Fasta file into a FastxIndex """

        self.fxi, self.fxi_file = create_fastx_index(fastx=self.fasta)

    def read_assembly_info(self) -> pandas.DataFrame:

        """ Read the Flye assembly information file """

        return pandas.read_csv(
            self.info_file, header=0, sep='\t',
            names=[
                'name',
                'length',
                'coverage',
                'circular',
                'repeat',
                'multiplicity'
            ], usecols=[0, 1, 2, 3, 4, 5]
        )


class KrakenProcessor(PoreLogger):

    """ Process results from Kraken2 """

    def __init__(
        self,
        report: Path,
        reads: Path,
        level: str = 'S',
        verbose: bool = True
    ):

        PoreLogger.__init__(
            self, level=logging.INFO if verbose else logging.ERROR
        )

        self.report = report
        self.reads = reads

        self.level = level
        self.host = 'Homo sapiens'  # for now always human

        # Reads below which assignments are considered contamination
        # Since throughput may vary a proportional threshold might not
        # be applicable. In addition, septic agent in Lachlan's paper
        # was identified with 18 reads. Default consider a minimum of
        # 5 reads evidence for real taxonomic assignments.

        # TODO: Improve contamination assignments!
        self.contamination_threshold = 5

        self.report_data = self.read_report()
        self.read_data = self.read_classifications()

        self.total_reads = self.report_data.loc[
            self.report_data['taxonomy'].isin(
                ('unclassified', 'root')
            ), 'reads'
        ].sum()

    def get_server_data(self):

        (bacteria_report, virus_report, contamination_report), \
            (bacteria_data, virus_data) = self.process_microbial()

        # Get server data format for donut visualization
        bacteria_server = self.report_to_json(df=bacteria_report)
        virus_server = self.report_to_json(df=virus_report)
        contamination_server = self.report_to_json(df=contamination_report)

        # Domain summary plot:
        host_percent, host_reads, _ = self.get_host()
        unclassified_percent, unclassified_reads, _ = self.get_unclassified()

        microbial_reads = bacteria_data[1] + virus_data[1]
        microbial_percent = bacteria_data[0] + virus_data[0]

        summary_server = [
            {
                'species': 'Host',
                'reads': host_reads,
                'percent': host_percent
            },
            {
                'species': 'Unclassified',
                'reads': unclassified_reads,
                'percent': unclassified_percent
            },
            {
                'species': 'Microbial',
                'reads': microbial_reads,
                'percent': microbial_percent,
            }
        ]

        return {
            'bacteria': bacteria_server,
            'viruses': virus_server,
            'contamination': contamination_server,
            'summary': summary_server
        }

    def get_unclassified(self) -> (float, int, pandas.DataFrame) or dict:

        """ Get unclassified reads and compute summaries from report

        :returns
            None, if no unclassified reads detected

            (percent reads, total reads, reads), if server_json = False
                reads are a subset of the read classifications from Kraken

            {percent: float, total: int, species: 'Unclassfied'} if server_json = True

        """

        unclassified = self.report_data[
            self.report_data['taxonomy'] == 'unclassified'
        ]

        nonsense = self.report_data[
            self.report_data['taxonomy'].isin(
                ('root', 'cellular organisms', 'other sequences')
            )
        ]

        if nonsense.empty:
            total_nonsense = 0
            percent_nonsense = 0.
        else:
            total_nonsense = int(
                nonsense['direct'].sum()
            )
            percent_nonsense = float(
                nonsense['direct'].sum() / self.total_reads
            )

        if unclassified.empty:
            total_unclassified = 0
            percent_unclassified = 0.
        else:
            total_unclassified = int(
                unclassified['reads']
            )
            percent_unclassified = float(
                unclassified['reads'] / self.total_reads
            ) + float(percent_nonsense)

        final_unclassified = total_unclassified + total_nonsense

        if final_unclassified == 0:
            return 0., 0, pandas.DataFrame()

        reads = pandas.concat((
            self.read_data[self.read_data['classified'] == 'U'],
            self.read_data[self.read_data['taxon'].isin(
                (
                 'root (taxid 1)',
                 'cellular organisms (taxid 131567)',
                 'other sequences (taxid 28384)'
                 )  # TODO GTDB compatibility [#6]
            )]
        ))

        return percent_unclassified, total_unclassified, reads

    def get_host(self) -> (float, int, pandas.DataFrame) or dict:

        """ Get host reads and compute summaries from report

        :returns percent reads, total reads, reads
            where reads are a subset of the read classifications from Kraken

        """

        host = self.report_data[
            self.report_data['taxonomy'] == self.host
        ]

        if host.empty:
            return 0., 0, []
        else:
            percent = float(host['reads'] / self.total_reads)
            total = int(host['reads'])

        read_idx = []
        for i, row in self.read_data.iterrows():
            name, taxon = row['read'], row['taxon']
            if self.host in taxon:
                read_idx.append(i)
        reads = self.read_data.iloc[read_idx]

        return percent, total, reads  # float, int, pandas.DataFrame

    def get_microbial(self) -> (tuple, tuple, tuple):

        """ Get microbial data for microbial domains from report

        Curently excludes eukaryotic pathogens; included:
        Bacteria, Viruses, Archaea

        :returns
            None, if no unclassified reads detected

            (percent reads, total reads, reads), if server_json = False
                reads are a subset of the read classifications from Kraken

            {percent: float, total: int, species: 'Host'} if server_json = True

        """

        bacteria, viruses, archaea = \
            self.get_subdomain(self.report_data, 'Bacteria'), \
            self.get_subdomain(self.report_data, 'Viruses'), \
            self.get_subdomain(self.report_data, 'Archaea')

        if bacteria is not None:
            bacteria_data = \
                float(bacteria.loc[0, 'reads'] / self.total_reads), \
                int(bacteria.loc[0, 'reads']), \
                bacteria
        else:
            self.logger.info(
                'Did not find any bacterial classifications in report.'
            )
            bacteria_data = (0., 0, pandas.DataFrame())

        if viruses is not None:
            viruses_data = \
                float(viruses.loc[0, 'reads'] / self.total_reads), \
                int(viruses.loc[0, 'reads']), \
                viruses
        else:
            self.logger.info(
                'Did not find any viral classifications in report.'
            )
            viruses_data = (0., 0, pandas.DataFrame())

        if archaea is not None:
            archaea_data = \
                float(archaea.loc[0, 'reads'] / self.total_reads), \
                int(archaea.loc[0, 'reads']), \
                archaea
        else:
            self.logger.info(
                'Did not find any archaeal classifications in report.'
            )
            archaea_data = (0., 0, pandas.DataFrame())

        return bacteria_data, viruses_data, archaea_data

    def process_microbial(self) -> (tuple, tuple) or None:

        """ Summarise Kraken report data for server output """

        # Get microbial
        bacteria_data, virus_data, archaea_data = self.get_microbial()

        # TODO: Identify contamination, needs work
        decontaminated, contamination = self.find_contamination(
            dict(
                Bacteria=bacteria_data[-1],
                Viruses=virus_data[-1],
                Archaea=archaea_data[-1]
            )
        )

        # Species classifications only!

        try:
            bacteria_report = decontaminated[0]
            bacteria_report = bacteria_report.loc[
              bacteria_report['level'] == self.level, :
            ]
        except IndexError:
            bacteria_report = pandas.DataFrame()

        try:
            virus_report = decontaminated[1]
            virus_report = virus_report.loc[
               virus_report['level'] == self.level, :
            ]
        except IndexError:
            virus_report = pandas.DataFrame()

        contamination_report = contamination.loc[
           contamination['level'] == self.level, :
        ]

        return (bacteria_report, virus_report, contamination_report), \
               (bacteria_data, virus_data)

    # Support methods

    def get_reads_from_report(self, report: pandas.DataFrame) -> pandas.DataFrame:

        read_idx = []
        allowed_taxa = report.taxonomy.tolist()
        for i, row in self.read_data.iterrows():
            name = row['taxon'].split('(')[0].strip()
            if name in allowed_taxa:
                read_idx.append(i)
        reads = self.read_data.iloc[read_idx].copy()

        reads.index.name = report.index.name

        return reads

    def find_contamination(self, data: dict) -> (tuple, pandas.DataFrame):

        """ Identify potential contaminants in Kraken report domains """

        decontaminated = []
        contaminated = []
        for name, df in data.items():
            if not df.empty:
                if name == 'Archaea':
                    # archaeal reads
                    contaminated.append(df)
                else:
                    # singleton reads
                    decon = df.loc[df['reads'] > self.contamination_threshold, :]
                    conta = df.loc[df['reads'] <= self.contamination_threshold, :]
                    decontaminated.append(decon)
                    contaminated.append(conta)

        return decontaminated, pandas.concat(contaminated)

    def report_to_json(self, df: pandas.DataFrame) -> list:

        if df.empty:
            return []
        else:
            return [
                {
                    'species': row['taxonomy'],
                    'percent': float(row['reads'] / self.total_reads),
                    'reads': row['reads']
                }
                for i, row in df.iterrows()
            ]

    def check_classifications(self, host_prefix='human'):

        # Get unclassified, check how many host, how many microbial

        upercent, utotal, ureads = self.get_unclassified()
        uhost, umicrobial = self._assess(ureads, host_prefix)

        self.logger.info(f'Unclassified reads: {uhost} host, {umicrobial} microbial')

        # Get reads classified as host, check how many host, microbial

        hpercent, htotal, hreads = self.get_host()
        hhost, hother = self._assess(hreads, host_prefix)

        htotal = hhost+hother

        self.logger.info(
            f'Reads classified as host: {hhost}/{htotal} true positives (host)'
        )
        self.logger.info(
            f'Reads classified as host: {hother}/{htotal} false positives (other)'
        )

    def assess_composition(self, compose_file: Path, host_prefix='human'):

        """ Assess a composed mixture for host / microbial classification

         Given a composition file matching the reads from the input
         classification, investigate how each species in the mixture
         is classified.

        :param compose_file: composition file JSON format from `np utils compose`
        :param host_prefix: the prefix usedi n the composition file as host reads

        """

        # Read the compose file for the report and

        with compose_file.open() as jin:
            composition = json.load(jin)

        if host_prefix not in composition.keys():
            raise ValueError(f'Host prefix {host_prefix} not in compose file.')

        return self._assess_composition(
            composition=composition, host_prefix=host_prefix
        )

    def _assess_composition(self, composition: dict, host_prefix: str):

        true_positives, negatives, unclassified, nonsense = {}, {}, {}, {}
        negative_data = []

        for i, row in self.read_data.iterrows():
            name, taxon = row['read'], row['taxon']

            try:
                prefix = name.split('_')[0]
            except IndexError:
                self.logger.error(
                    'Something is wrong with the compose file prefixes.'
                )
                raise

            try:
                prefix_data = composition[prefix]
            except KeyError:
                self.logger.error(
                    f'Could not find {prefix} in composition.'
                )
                raise

            taxon_name = prefix_data['taxon_name']

            if taxon_name in taxon:
                if prefix not in true_positives.keys():
                    true_positives[prefix] = 1
                else:
                    true_positives[prefix] += 1
            elif 'unclassified' in taxon:
                if prefix not in unclassified.keys():
                    unclassified[prefix] = 1
                else:
                    unclassified[prefix] += 1
            elif 'root' in taxon or 'cellular organism' in taxon:
                if prefix not in nonsense.keys():
                    nonsense[prefix] = 1
                else:
                    nonsense[prefix] += 1
            else:
                if prefix not in negatives.keys():
                    negatives[prefix] = 1
                else:
                    negatives[prefix] += 1

                negative_data.append((name, taxon))

        totals = {}  # Fill in null counts and get total count
        for prefix in composition.keys():
            for count in (true_positives, negatives, unclassified, nonsense):
                if prefix not in count.keys():
                    count[prefix] = 0

                if prefix not in totals.keys():
                    totals[prefix] = count[prefix]
                else:
                    totals[prefix] += count[prefix]

        negative_data = pandas.DataFrame(negative_data)

        summary_data = (true_positives, negatives, unclassified, nonsense, totals)

        summary = pandas.DataFrame(
            summary_data,
            index=['positives', 'negatives', 'unclassified', 'nonsense', 'totals']
        )

        summary_proportions = (summary / summary.loc['totals', :])
        summary_proportions = summary_proportions.drop('totals', axis=0)

        # Assess negative data

        microbial, conflict, domain = [], [], []

        negative_data.columns = ['prefix', 'predicted']
        for i, row in negative_data.iterrows():
            prefix = row['prefix'].split('_')[0]

            if prefix == host_prefix:
                true_domain = 'Host'
            else:
                true_domain = 'Microbial'  # Bacteria, Archaea, Viruses !

            predicted_domain = self.get_superdomain(
                taxon_name=row['predicted']
            )

            microbial.append(true_domain)
            domain.append(predicted_domain)

            if predicted_domain in ('Bacteria', 'Archaea', 'Viruses') \
                    and true_domain == 'Microbial':
                conflict.append('no')
            elif predicted_domain == 'Viruses' and 'Human gammaherpesvirus 4' in row['predicted']:
                conflict.append('retro')
            elif predicted_domain == 'Viruses' and true_domain == 'Host':
                conflict.append('yes')
            elif predicted_domain in ('Bacteria', 'Archaea', 'Viruses') \
                    and true_domain == 'Host':
                conflict.append('yes')
            elif predicted_domain == 'Host' and true_domain == 'Microbial':
                conflict.append('yes')
            else:
                conflict.append('unknown')

        negative_data['predicted_domain'] = domain
        negative_data['prefix_class'] = microbial
        negative_data['conflict'] = conflict

        negative_data.index.name = 'negative'
        negative_data.sort_values(
            ['conflict'], inplace=True, ascending=False
        )

        return summary, summary_proportions, negative_data

    def create_summary_plot(
        self,
        summary,
        summary_proportions,
        negative_data,
        prefix: str = 'summary',
        palette: str = 'Blues'
    ):

        fig, axes = plt.subplots(
            nrows=3, ncols=1, figsize=(1 * 12, 3 * 9)
        )

        fig.subplots_adjust(hspace=0.8)

        if len(axes) != len(summary_proportions.columns):
            raise ValueError('Axes must be of same length as columns.')

        summary_proportions.loc['negatives', :] = \
            summary_proportions.loc['nonsense', :] + \
            summary_proportions.loc['negatives', :]

        summary_proportions.drop('nonsense', axis=0, inplace=True)

        for i, column in enumerate(summary_proportions):
            c = summary_proportions[column]
            labels = [
                name + f" [{round(c[i]*100, 2)}%]"
                for i, name in enumerate(c.index)
            ]
            self.plot_annotated_donut(
                labels=labels, data=c.tolist(),
                palette=palette, ax=axes[i], title=column
            )

        plt.tight_layout()
        fig.savefig(f'{prefix}.png')

        self.plot_sankey(data=negative_data, prefix=prefix)

        negative_data.to_csv(
            f'{prefix}.negatives.tsv', sep='\t', index=False
        )
        summary_proportions.to_csv(
            f'{prefix}.props.tsv', sep='\t', index=True
        )
        summary.to_csv(
            f'{prefix}.reads.tsv', sep='\t', index=True
        )

    def get_superdomain(self, taxon_name: str):

        """ Get a domain-specific subset of the report dataframe"""

        df = self.report_data.reset_index()

        try:
            taxon_index = [
                i for i, taxon in enumerate(
                    df['taxonomy']
                ) if taxon in taxon_name
            ][0]
        except IndexError:
            self.logger.error('Something went wrong in _get_super_domain')
            raise

        super_levels = df['level'][0:taxon_index].tolist()
        super_levels.reverse()

        domain_index = 0
        for i, level in enumerate(super_levels):
            if level == 'D' and i > 0:
                domain_index = taxon_index - i - 1
                break

        domain_row = df.iloc[domain_index, :]

        return str(domain_row['taxonomy'])

    def get_subdomain(
        self, report: pandas.DataFrame, taxon: str
    ) -> pandas.DataFrame or None:

        """ Get a domain-specific subset of the report dataframe"""

        df = report.reset_index()

        try:
            taxon_index = int(df[df['taxonomy'] == taxon].index.values)
        except TypeError:
            self.logger.warning(f'Could not get subdomain of {taxon}')
            return None

        domain_index = 0
        for i, level in enumerate(
                df['level'][taxon_index:]
        ):
            if level == 'D' and i > 0:
                domain_index = taxon_index + i
                break

        if domain_index == 0:
            return report.iloc[taxon_index:, :].reset_index()
        else:
            return report.iloc[taxon_index:domain_index, :].reset_index()

    @staticmethod
    def _assess(reads: pandas.DataFrame, host_prefix: str):

        host = 0
        microbial = 0
        for name in reads['read']:
            if name.startswith(host_prefix):
                host += 1
            else:
                microbial += 1

        return host, microbial

    def read_report(self):

        df = pandas.read_csv(
            self.report, header=None, sep="\t", names=[
                "percent", "reads", "direct", "level", "taxid", "taxonomy",
            ]
        )
        df.taxonomy = df.taxonomy.str.strip()  # strip indentation

        return df

    def read_classifications(self):
        return pandas.read_csv(
            self.reads, header=None, sep="\t", names=[
                'classified', 'read', 'taxon', 'node', 'path'
            ]
        )

    @staticmethod
    def plot_sankey(
        data: pandas.DataFrame,
        true: str = 'prefix_class',
        predicted: str = 'predicted_domain',
        palette: str = 'Blues',
        prefix: str = None,
    ):

        categories = data[true].unique().tolist() + \
                     data[predicted].unique().tolist()

        colors = {
            categories[i]: c for i, c in enumerate(
                sns.color_palette(palette, len(categories))
            )
        }

        sankey.sankey(
            data[true],
            data[predicted],
            aspect=40,
            colorDict=colors,
            fontsize=12,
            figure_name=prefix if prefix is not None else 'negatives'
        )

    @staticmethod
    def plot_annotated_donut(
        labels: [str], data: [float], ax: None, palette="Blues", title: str = 'title'
    ):

        """ Donut chart with labels """

        if len(labels) != len(data):
            raise ValueError('Data and label lists most be of the same length.')

        color = sns.color_palette(palette, len(data))

        ax.set_title(title)

        wedges, texts = ax.pie(
            data, wedgeprops=dict(width=0.5), startangle=-40, colors=color
        )

        bbox_props = dict(
            boxstyle="square,pad=0.3", fc="w", ec="k", lw=0
        )

        kw = dict(
            arrowprops=dict(arrowstyle="-"),
            bbox=bbox_props, zorder=0, va="center"
        )

        for i, p in enumerate(wedges):

            ang = (p.theta2 - p.theta1) / 2. + p.theta1
            y = np.sin(np.deg2rad(ang))
            x = np.cos(np.deg2rad(ang))

            horizontalalignment = {-1: "right", 1: "left"}[int(np.sign(x))]
            connectionstyle = "angle,angleA=0,angleB={}".format(ang)
            kw["arrowprops"].update({"connectionstyle": connectionstyle})

            ax.annotate(
                labels[i], xy=(x, y), **kw,
                xytext=(1.35 * np.sign(x), 1.4 * y),
                horizontalalignment=horizontalalignment,
                fontsize=16
            )
