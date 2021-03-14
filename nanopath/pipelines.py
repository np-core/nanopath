import logging
import copy
import shutil
import re

from tqdm import tqdm

import dataclasses
import pandas
import json
import math

from dataclasses import dataclass
from typing import Callable
from pathlib import Path
from nanopath.utils import PoreLogger, run_cmd
from nanopath.tools.kraken import merge_reports
from nanopath.pathfinder import retain_files, get_id_from_fname, get_subdict
from pandas.core.groupby import DataFrameGroupBy
from matplotlib import pyplot as plt

import multiprocessing as mp
import seaborn as sns


class PipelineError(Exception):
    pass


#################################
# Data containers for Nextflows #
#################################


class ResultData:
    """ Methods for data containers. """

    def __init__(self):

        self.iid: pandas.DataFrame = pandas.DataFrame(
            columns=['iid']
        )
        self.fasta: pandas.DataFrame = pandas.DataFrame(
            columns=['fasta']
        )
        self.fastq: pandas.DataFrame = pandas.DataFrame(
            columns=['forward', 'reverse']
        )

    def __copy__(self):
        """ Overwritten in subclasses """

    def __iter__(self):

        for attr, df in self.__dict__.items():
            yield attr, df

    def __add__(self, other):

        for attr, data in self:
            data_other = getattr(other, attr)

            merged = pandas.concat(
                [data, data_other], sort=False
            )

            setattr(self, attr, merged)

        self.update_iid()

        return self

    def __getitem__(self, iid):

        return self.iid[iid]

    def update_iid(self, iids=None):

        """ Update list of unique IIDs in current data container """

        if iids is None:
            iids = pandas.Series(
                [iid for attr, df in self for iid in df.index.unique().tolist()]
            ).unique()  # Union

        self.iid = pandas.DataFrame(data={'iid': iids}, index=iids)

    def remove(self, remove_dict: dict, retain: bool = False):

        for process, df in self:

            try:
                remove_dict[process]
            except KeyError:
                # Skip self.iid, self.fasta, self.fastq
                continue

            if retain:
                removed_df = df[
                    df.index.isin(
                        remove_dict[process]
                    )
                ]
            else:
                removed_df = df.drop(
                    remove_dict[process], errors='ignore'
                )

            setattr(self, process, removed_df)

        self.update_iid()

    def complete(
        self,
        at: list or None = ('kraken', 'mlst')
    ):

        """ Return data object with the set intersection between all IIDs
        in the analyses results; usually after filtering, can be narrowed
        to only some results by default to Kraken and MLST """

        if at is not None:
            allowed = at
        else:
            allowed = [process for process, _ in self]

        iids = [
            set(df.index.tolist()) for process, df in self
            if len(df) > 0 and process in allowed
        ]

        intersect = list(set.intersection(*iids))

        self.remove(
            {process: intersect for process, _ in self},
            retain=True
        )

    def write(self, outdir: str or Path) -> None:

        outdir = Path(outdir) \
            if isinstance(outdir, str) else outdir

        outdir.mkdir(parents=True, exist_ok=True)

        for process, data in self:
            data.to_csv(
                (outdir / f'{process}.csv'),
                index=True, header=True, sep=','
            )

    def read(self, results: str or Path):

        results = Path(results) \
            if isinstance(results, str) else results

        for file in results.glob('*.csv'):
            setattr(self, file.stem, pandas.read_csv(
                file, index_col=0, header=0, sep=','
            ))

        self.update_iid()

    def subset(self, iids, inplace=False):

        """ Subset the Data by a list of isolate IDs.

        If IDs are not in the data frame index, a row of null values
        for this ID is introduced into the data frame.
        """

        if inplace:
            sub = self
        else:
            sub = copy.copy(self)

        for attr, df in self:
            subset = df[df.index.isin(iids)]
            missing = list(set(iids).difference(set(subset.index)))
            if missing:
                subset = pandas.concat(
                    (
                        subset,
                        pandas.DataFrame(index=missing, columns=df.columns)
                     )
                )

            setattr(sub, attr, subset)

        sub.update_iid()

        return sub

    def common(self) -> list:

        """ Get the intersection between all IIDs in the Data """

        return list(
            set(self.iid).intersection(*[df.index.tolist() for _, df in self])
        )

    def groupby(self, attr, column, set_index=False, **kwargs) -> [DataFrameGroupBy]:

        """ Map the column values to each corresponding index
        in all other data attributes and group the data frames
        by the column values. This allows e.g. to group all data
        sequence type etc.
        """

        if not hasattr(self, attr):
            raise AttributeError(
                f'ResultData has no attribute: {attr}'
            )

        # Map sequence type to index in each DataFrame as new column
        # and then group every attribute data frame by this column:

        data = getattr(self, attr)
        data_series = data[column]

        for attr, df in self:
            df[column] = df.index.map(data_series)

            if set_index:
                df = df.set_index([df.index, column])
                setattr(self, attr, df)

            groups = df.groupby(by=column, **kwargs)
            yield groups, attr

    def select(self, attr, column, min_count=None,
               sample=None, values=None):

        data = getattr(self, attr)

        # TODO: Both

        if values:
            iids = data[
                data[column].isin(values)
            ].index.unique().tolist()

        elif min_count:
            counts = data[column].value_counts()
            largest = counts[counts > min_count]

            subset = largest.index.to_series()  # Index: unique values counted

            sub = data[
                data[column].isin(subset)
            ]

            iids = sub.index.unique().tolist()

            if sample:

                iids = pandas.concat([
                    group.sample(n=sample) for i, group in sub.groupby(column)
                ]).index

        else:
            iids = data[column].index.unique().tolist()

        return self.subset(iids)

    def link_fasta(self, fdir: str or Path, symlink: bool = True,
                   progbar: bool = True, index: dict = None):

        fdir = self._check_dir(fdir)

        # Can we predict read length distribution from a smear
        # of gel run of HMW DNA

        for fasta in tqdm(
                self.fasta.fasta,
                disable=not progbar,
                total=len(self.fasta)
        ):
            if index:
                try:
                    name = index[Path(fasta).stem] + '.fasta'
                except TypeError:
                    print(index[Path(fasta).stem])
                    continue
                except KeyError:
                    print('Could not find entry in FASTA Index.')
                    continue
            else:
                name = Path(fasta).name

            if symlink:
                (fdir / name).symlink_to(fasta)
            else:
                shutil.copy(
                    fasta, str(fdir / name)
                )

    def link_fastq(self, fdir: str or Path, symlink: bool = True,
                   progbar: bool = True):

        fdir = self._check_dir(fdir)

        for fwd, rv in tqdm(
            self.fastq.itertuples(index=False),
            disable=not progbar,
            total=len(self.fastq)
        ):
            for fastq in (fwd, rv):
                if symlink:
                    (fdir / Path(fastq).name).symlink_to(fastq)
                else:
                    shutil.copy(
                        fastq, str(fdir)
                    )

    @staticmethod
    def _check_dir(fdir: str or Path) -> Path:

        fdir = Path(fdir) if isinstance(fdir, str) else fdir

        if not fdir.exists():
            fdir.mkdir(parents=True)

        return fdir

    # File operations

    def _add_fasta(self, fasta_dir: Path, extension: str = '.fasta'):

        """ Add file paths to data frame in attribute: `fasta` """

        # Clear FASTA, weirdly in loops files accumulate
        # even if new instance of this is created

        self.fasta = pandas.DataFrame()

        for fpath in fasta_dir.glob(f'*{extension}'):
            iid = fpath.name.replace(extension, '')
            self.fasta.at[iid, 'fasta'] = fpath

    def _add_fastq(self, fastq_dir: Path, forward_tail: str = '_1',
                   reverse_tail: str = '_2', extension: str = '.fastq.gz'):

        self.fastq = pandas.DataFrame()

        for fpath in fastq_dir.glob(f'*{forward_tail}{extension}'):
            iid = fpath.name.replace(f'{forward_tail}{extension}', '')
            if iid in self.fastq.index:
                self.fastq.at[iid, 'forward'] = fpath
            else:
                pass

        for fpath in fastq_dir.glob(f'*{reverse_tail}{extension}'):
            iid = fpath.name.replace(f'{reverse_tail}{extension}', '')
            if iid in self.fastq.index:
                self.fastq.at[iid, 'reverse'] = fpath
            else:
                pass

    def get_file_paths(
        self,
        result_path: str or Path,
        fasta_dir: str = None,
        fastq_dir: str = None
    ) -> None:

        if fasta_dir is None and fastq_dir is None:
            raise ValueError(
                'Please specify a directory name for FASTA or FASTQ.'
            )

        print('Getting file paths:', result_path)

        Path(result_path) if isinstance(result_path, str) else result_path

        if fasta_dir:
            self._add_fasta(result_path / fasta_dir)
        if fastq_dir:
            self._add_fastq(result_path / fastq_dir)

        self.update_iid()

    def by_iid(self) -> pandas.DataFrame:

        """ Returns a bool summary of IIDs in each process """

        exclude_attrs = ('fasta', 'fastq', 'iid')

        df = pandas.DataFrame(
            index=self.iid.iid,
            columns=[attr for attr, _ in self if attr not in exclude_attrs]
        )

        for attr, data in self:
            if attr in exclude_attrs:
                continue
            for iid in self.iid.iid:
                if iid in data.index:
                    df.at[iid, attr] = True
                else:
                    df.at[iid, attr] = False

        return df


@dataclass
class SurveyProcessSetting:
    """ Template for processe pipeline data """

    files: list = None
    parse_func: Callable = None
    process_func: Callable = None
    remove: str = '.extension'


class SurveyData(ResultData):

    def __init__(self):

        ResultData.__init__(self)

        self.mlst = pandas.DataFrame()
        self.kraken = pandas.DataFrame()
        self.mash = pandas.DataFrame()
        self.kleborate = pandas.DataFrame()
        self.sccion = pandas.DataFrame()

        self.abricate_resistance = pandas.DataFrame()
        self.abricate_virulence = pandas.DataFrame()
        self.abricate_plasmid = pandas.DataFrame()

        self.mykrobe_phenotype = pandas.DataFrame()
        self.mykrobe_genotype = pandas.DataFrame()
        self.mykrobe_lineage = pandas.DataFrame()

    def __copy__(self):

        sd = SurveyData()
        for attr, data in self:
            setattr(sd, attr, data)
        return sd

    @property
    def empty(self):

        check = [data.empty for attr, data in self]

        if all(check):
            return True
        else:
            return False


class SurveyResult:

    """ Parse results from: pf-core/pf-survey """

    def __init__(self, path: Path = None):

        self.path = path
        self.data: SurveyData = SurveyData()

    def __add__(self, other):

        # Enable sum() with __radd__
        # by skipping starting condition:
        if other == 0:
            return self
        else:
            if not isinstance(other, SurveyResult):
                raise TypeError('Only other SurveyResult objects can be added.')
            self.path = None
            self.data += other.data

            return self

    def __radd__(self, other):

        return self.__add__(other)

    # Main access methods

    def parse(
        self
    ):
        sp = SurveyParser(path=self.path)

        for process, df in sp.parse():
            setattr(self.data, process, df)

        # self.data.update_iid()

        # print(self.data)

    def filter(
        self,
        taxon: int or str = None,
        level: str = "S",
        species: str = None,
        purity: float = 0.8,
        contamination: float = 0.02,
        gene_identity: float = 0.9,
        gene_coverage: float = 0.9,
        complete_lineage: bool = True
    ):

        sf = SurveyFilter(
            taxon=taxon,
            level=level,
            species=species,
            purity=purity,
            contamination=contamination,
            gene_identity=gene_identity,
            gene_coverage=gene_coverage,
            complete_lineage=complete_lineage
        )

        self.data.remove(
            dict(
                kraken=sf.clean_kraken(self.data.kraken),
                mlst=sf.clean_mlst(self.data.mlst)
            )
        )


#################
# Survey Filter #
#################

@dataclasses.dataclass
class SurveyFilter:

    purity: float = 0.8
    taxon: int or str = None
    level: str = "S"
    species: str = None
    contamination: float = 0.02
    gene_identity: float = 0.9
    gene_coverage: float = 0.9
    complete_lineage: bool = True

    def clean_kraken(
        self, data: pandas.DataFrame
    ) -> list:

        """ Clean and condense the taxonomic read assignments by Kraken2

        Inspects each isolate's results from the agglomerated data frame
        and identify the following:


            .. py:attr:`pathfinder.pipelines.SurveyFilter.contamination`

                Threshold for removing contamination with this percentage of
                reads classified in top species besides the most represented
                species.

            .. py:attr:`pathfinder.pipelines.SurveyFilter.purity`

                Threshold for retaining samples with at least this percentage
                of reads assigned to the top species.

        :returns A list of isolate identifiers to be removed at the
            end of the cleaning process.

        """

        isolates = data.groupby(by=data.index)

        uncertain, contaminated, misidentified = 0, 0, 0,

        to_remove = list()
        for _, g in isolates:

            if not g.empty:
                iid = g.index[0]

                g = g.loc[g['level'] == self.level]

                g = g.nlargest(3, 'percent')

                if g['percent'].iloc[0] < self.purity * 100:
                    to_remove.append(iid)
                    uncertain += 1

                if (g['percent'][1:] > self.contamination * 100).any():
                    contaminated += 1
                    to_remove.append(iid)

                if self.taxon:
                    if isinstance(self.taxon, str):
                        column = 'taxonomy'
                    else:
                        column = 'taxid'

                    if g[column].iloc[0] != self.taxon:
                        to_remove.append(iid)
                        misidentified += 1

        return list(
            pandas.Series(to_remove).unique()
        )

    def clean_mlst(self, data: pandas.DataFrame) -> list:

        if self.complete_lineage:
            to_remove = data[data['sequence_type'] == '-'].index.tolist()

            if self.species:
                to_remove += data[data['species'] != self.species].index.tolist()

            return to_remove
        else:
            return list()


#################
# Survey Parser #
#################


class SurveyParser:

    """ Collection of parsing methods and data frame headers
    to make life easier by using an instance of the SurveyParser """

    def __init__(self, path: Path):

        # Files to parse

        self.path = path

        # DataFrame columns: Surveys

        self.mlst_columns = [
            "file", "species", "sequence_type",
        ]
        self.mash_columns = [
            "file", "ref", "dist", "p-value", "match",
        ]
        self.kraken_report_columns = [
            "percent", "reads", "direct", "level", "taxid", "taxonomy",
        ]
        self.prokka_columns = [
            "locus_tag", "feature", "length", "gene", "ec", "cog", "product",
        ]

        self.sccion_columns = [
            'name', 'species', 'mlst', 'meca', 'pvl', 'spa', 'scc', 'resistance', 'plasmid'
        ]

        self.abricate_columns = [
            "file", "sequence", "start", "end", "gene", "coverage_bases",
            "coverage_map", "gaps",  "coverage", "identity", "database",
            "accession", "product",
        ]
        self.kleborate_columns = [
            "strain", "species", "st", "virulence_score", "resistance_score",
            "Yersiniabactin", "YbST", "Colibactin", "CbST", "Aerobactin",
            "AbST", "Salmochelin", "SmST", "rmpA", "rmpA2", "wzi",  "K_locus",
            "K_locus_confidence", "O_locus", "O_locus_confidence", "AGly",
            "Col", "Fcyn", "Flq", "Gly", "MLS", "Ntmdz", "Phe", "Rif", "Sul",
            "Tet", "Tmt", "Bla", "Bla_Carb", "Bla_ESBL", "Bla_ESBL_inhR",
            "Bla_broad", "Bla_broad_inhR",
        ]

    # Main access parsing method:

    def parse(self):

        """ Parse - Process - Clean """

        processes = {
            'mykrobe_lineage': SurveyProcessSetting(
                files=list(
                    (self.path / 'mykrobe').glob('*.json')
                ),
                parse_func=self.parse_mykrobe_lineage,
                remove='.json'
            ),
            'kraken': SurveyProcessSetting(
                files=list(
                    (self.path / 'kraken').glob('*.report')
                ),
                parse_func=self.parse_kraken,
                process_func=self.process_kraken,
                remove='.report'
            ),
            'mlst': SurveyProcessSetting(
                files=list(
                    (self.path / 'mlst').glob('*.tab')
                ),
                parse_func=self.parse_mlst,
                remove='.tab'
            ),
            'mash': SurveyProcessSetting(
                files=list(
                    (self.path / 'mash').glob('*.mash.tab')
                ),
                parse_func=self.parse_mash,
                remove='.mash.tab'
            ),
            'mykrobe_phenotype': SurveyProcessSetting(
                files=list(
                    (self.path / 'mykrobe').glob('*.json')
                ),
                parse_func=self.parse_mykrobe_susceptibility,
                remove='.json'
            ),
            'mykrobe_genotype': SurveyProcessSetting(
                files=list(
                    (self.path / 'mykrobe').glob('*.json')
                ),
                parse_func=self.parse_mykrobe_genotype,
                remove='.json'
            ),
            'abricate_resistance': SurveyProcessSetting(
                files=list(
                    (self.path / 'abricate' / 'resfinder').glob('*.tab')
                ),
                parse_func=self.parse_abricate,
                process_func=self.process_abricate,
                remove='.tab'
            ),
            'abricate_virulence': SurveyProcessSetting(
                files=list(
                    (self.path / 'abricate' / 'vfdb').glob('*.tab')
                ),
                parse_func=self.parse_abricate,
                process_func=self.process_abricate,
                remove='.tab'
            ),
            'abricate_plasmid': SurveyProcessSetting(
                files=list(
                    (self.path / 'abricate' / 'plasmidfinder').glob('*.tab')
                ),
                parse_func=self.parse_abricate,
                process_func=self.process_abricate,
                remove='.tab'
            ),
            'kleborate': SurveyProcessSetting(
                files=list(
                    (self.path / 'kleborate').glob('*.tab')
                ),
                parse_func=self.parse_kleborate,
                remove='.tab'
            ),
            'sccion': SurveyProcessSetting(
                files=list(
                    (self.path / 'sccion').glob('*.txt')
                ),
                parse_func=self.parse_sccion,
                remove='.txt'
            ),
        }

        for process, process_setting in processes.items():
            try:
                yield process, self._parse_and_process(
                    **dataclasses.asdict(process_setting),
                    process=process
                )
            except ValueError:
                print(f'Could not process: {process}')
                continue

    # Private helper methods for parsing

    def _parse_and_process(
        self,
        files: [],
        parse_func: Callable = None,
        process_func: Callable or None = None,
        remove: str or list = None,
        subset: str = None,
        process: str = None,
        *args, **kwargs
    ) -> dict or pandas.DataFrame:

        """Anonymous function to parse results based on a list of files,
        the name components to be stripped for extracting IDs, and the
        parse method of this class to return a dictionary of:

        :param files:
            list of input files from glob for parsing.
        :param remove:
            str or list of str, remove from filename for ID.
        :param parse_func:
            parse method of this class for parsing.
        :param process_func:
            process of data frame after reading
        :param subset:
            parse only files that contain substring

        :returns parsed results (dict) or aggregated results (dataframe)
        """

        files = list(files)

        if subset:
            files = retain_files(files, retain=subset)

        # Parse results from many files by ID:
        data = {
            get_id_from_fname(file, remove=remove): parse_func(file, *args)
            for file in tqdm(files, desc=f'{process}')
        }

        # Aggregate the data into DataFrames:
        data = self._aggregate(data)

        if process_func:
            data = process_func(data, **kwargs)

        return data

    def _aggregate(self, parsed_data):
        """ Aggregate results from the parsed data"""

        ids, dfs = self._split_data(parsed_data)

        df = pandas.concat(
            [
                df.assign(
                    id=[ids[i] for _ in range(len(df))]
                )
                for i, df in enumerate(dfs)
            ], sort=False).set_index('id')

        df.index.name = "ID"

        return df

    @staticmethod
    def _split_data(data) -> iter:
        """Split a data dictionary (keys: ids, values: data)
         into two id-sorted lists """

        return zip(*[(ide, data[ide]) for ide in sorted(data)])

    # Public helper methods for result file parsing

    def parse_kraken(self, file) -> pandas.DataFrame:
        """Parse single report output from Kraken"""
        return pandas.read_csv(
            file, header=None, sep="\t", names=self.kraken_report_columns
        )

    def parse_mash(self, file) -> pandas.DataFrame:
        """Parse single output from Mash"""
        return pandas.read_csv(
            file, sep="\t", header=None, names=self.mash_columns
        )

    def parse_abricate(self, file) -> pandas.DataFrame:
        """Parse single output from Abricate"""
        return pandas.read_csv(
            file, sep="\t", header=0, names=self.abricate_columns
        )

    def parse_prokka(self, file, gff=True) -> pandas.DataFrame:
        """Parse single output feature table from Prokka"""

        if gff:
            return self._parse_prokka_gff(file)
        else:
            return pandas.read_csv(
                file, sep="\t", header=0, names=self.prokka_columns
            )

    def parse_sccion(self, file) -> pandas.DataFrame:

        return pandas.read_csv(
            file, sep="\t", header=None, names=self.sccion_columns
        )

    @staticmethod
    def _parse_prokka_gff(file):

        data = []
        with open(file, 'r') as infile:
            for line in infile:
                if line.startswith('##'):
                    pass
                elif line.startswith('##FASTA'):
                    return 0
                else:
                    data.append(
                        line.strip().split('\t')
                    )

        return pandas.DataFrame(data, columns=[
            'contig', 'inference', 'type', 'start', 'stop', 'unknown_1',
            'strand', 'unknown_2', 'description'
        ])

    def parse_mlst(self, file) -> pandas.DataFrame:
        """Parse single output from mlst (github.com/tseemann/mlst)"""

        df = pandas.read_csv(file, sep="\t", header=None)

        df.columns = self.mlst_columns + [
            str(i) for i, col in enumerate(df.columns, 1)
        ][:-3]

        return df

    def parse_kleborate(self, file) -> pandas.DataFrame:
        """ Parse single file output from Kleborate """

        return pandas.read_csv(
            file, sep='\t', header=0
        )

    def parse_mykrobe_susceptibility(self, file) -> pandas.DataFrame:

        susceptibilities, _ = self._parse_mykrobe(file)

        return susceptibilities

    def parse_mykrobe_genotype(self, file) -> pandas.DataFrame:

        _, genotype = self._parse_mykrobe(file)

        return genotype

    def parse_mykrobe_lineage(self, file) -> pandas.DataFrame:

        return self._parse_mykrobe(file, lineage=True)

    @staticmethod
    def _parse_mykrobe(
        file, lineage=False
    ) -> (pandas.DataFrame, pandas.DataFrame) or pandas.DataFrame:
        """Parse single output from MykrobePredictor, only parse
        called resistances and the genes that called them.
        """

        with open(file, "r") as infile:
            mykrobe = json.load(infile)

            if lineage:
                # TODO make sure order highest if multiple?
                keys = [
                    k for k, v in list(
                        get_subdict("lineage", mykrobe)
                    )[0].items()
                ]
                return pandas.DataFrame.from_dict({'lineage': keys})
            else:
                susceptibility = list(
                    get_subdict("susceptibility", mykrobe)
                )[0]

                data = {
                    drug: {
                        "susceptibility": data["predict"],
                        "genotype": list(
                            data["called_by"].keys()
                        )
                        if "called_by" in data.keys() else None
                    }
                    for drug, data in susceptibility.items()
                }

                susceptibilities = pandas.DataFrame.from_dict(
                    {
                        drug: list(dat["susceptibility"])
                        for drug, dat in data.items()
                    }
                )

                genotypes = pandas.DataFrame.from_dict(
                    {
                        drug: [None] if not dat["genotype"]
                        else [",".join(dat["genotype"])]
                        for drug, dat in data.items()
                     }
                )

                return susceptibilities, genotypes

    # Public methods for processing

    @staticmethod
    def process_abricate(df):
        # Remove alleles from gene names:
        df.gene = df.gene.str.rsplit('_').str[0]
        return df

    @staticmethod
    def process_kraken(df):
        df.taxonomy = df.taxonomy.str.strip()
        return df


class AssemblyPipeline(PoreLogger):

    """ Collect summary data of hybrid assembly and genotyping pipeline """

    def __init__(self, path: Path, outdir: Path):

        PoreLogger.__init__(self, level=logging.INFO, name="AssemblyCollector")

        self.mp = MasterParser()

        self.path = path
        self.outdir = outdir

        self.outdir.mkdir(parents=True, exist_ok=True)

        self.components = {
            'illumina': self.path / 'illumina',
            'ont': self.path / 'ont',
            'hybrid': self.path / 'hybrid',
            'unicycler': self.path / 'unicycler'
        }

    def collect_dnadiff(self, exclude: list = None, string_sort: bool = False):

        dnadiff_path = self.components['hybrid'] / 'dnadiff'

        # Medaka polished ONT only assembly
        ont_medaka_vs_ref = ".medaka.report"
        # Hybrid assembly from polished assembly corrected with Illumina
        hybrid_medaka_vs_ref = ".medaka_hybrid.report"
        unicycler_vs_ref = ".unicycler.report"

        ont_files = dnadiff_path.glob("*"+ont_medaka_vs_ref)
        hybrid_medaka_files = dnadiff_path.glob("*"+hybrid_medaka_vs_ref)
        unicycler_files = dnadiff_path.glob("*"+unicycler_vs_ref)

        try:
            df_ont = pandas.concat(
                [
                    self.mp.parse_dnadiff(file=f, name_strip=ont_medaka_vs_ref)
                    for f in list(ont_files)
                ]
            )
        except ValueError:
            self.logger.info('Could not detect files for Dnadiff (ONT)')
            return

        try:
            df_hybrid_medaka = pandas.concat(
                [
                    self.mp.parse_dnadiff(file=f, name_strip=hybrid_medaka_vs_ref)
                    for f in list(hybrid_medaka_files)
                ]
            )
        except ValueError:
            self.logger.info('Could not detect files for Dnadiff (Hybrid)')
            return

        try:
            df_unicycler = pandas.concat(
                [
                    self.mp.parse_dnadiff(file=f, name_strip=unicycler_vs_ref)
                    for f in list(unicycler_files)
                ]
            )
        except ValueError:
            self.logger.info('Could not detect files for Dnadiff (Unicycler)')
            return

        if string_sort:
            df_ont = self.str_num_sort_col(df_ont, 'name')
            df_hybrid_medaka = self.str_num_sort_col(df_hybrid_medaka, 'name')
            df_unicycler = self.str_num_sort_col(df_unicycler, 'name')
        else:
            df_ont = df_ont.sort_values('name')
            df_hybrid_medaka = df_hybrid_medaka.sort_values('name')
            df_unicycler = df_unicycler.sort_values('name')

        df_ont['branch'] = [
            'ont_medaka' for _ in range(df_ont.__len__())
        ]
        df_hybrid_medaka['branch'] = [
            'hybrid_medaka' for _ in range(df_hybrid_medaka.__len__())
        ]
        df_unicycler['branch'] = [
            'hybrid_unicycler' for _ in range(df_unicycler.__len__())
        ]

        combined = pandas.concat(
            (df_ont, df_hybrid_medaka, df_unicycler), axis=0
        )

        if exclude:
            combined = combined.loc[~combined['name'].isin(exclude), :]

        combined.to_csv(self.outdir / 'dnadiff.tsv', sep='\t', index=False)

        return combined

    @staticmethod
    def str_num_sort(l: list):

        return sorted(l, key=lambda x: int(re.findall(r"\d+", x)[0]))

    def str_num_sort_col(self, df: pandas.DataFrame, col: str = None):

        try:
            names_sorted = sorted(
                df[col].tolist(), key=lambda x: int(
                    re.findall(r"\d+", x)[0]
                )
            )
            df = df.set_index('name') \
                .reindex(names_sorted).reset_index(level=0)
        except TypeError:
            self.logger.info(
                'Could not extract index using regex from'
                ' name string. Skipping sorting.'
            )

        return df

    def collect_genotypes(
        self, component: str = 'ont', name_index_regex: str = r"\d+", exclude: list = None
    ):

        genotype_path = self.components[component] / 'genotypes'

        if component == 'ont':
            fselect = "*.ont.tab"  # Medaka polished
            fstrip = [".medaka.fasta"]
        elif component == "hybrid":
            fselect = "*.hybrid.tab"  # Hybrid assembly without Medaka pre-polish
            fstrip = [".medaka.hybrid.fasta"]
        elif component == "illumina":
            fselect = "*.illumina.tab"  # Illumina assembly
            fstrip = [".assembly.fasta", ".fasta"]
        elif component == "unicycler":
            fselect = "*.unicycler.tab"  # Unicycler assembly
            fstrip = [".unicycler.fasta", ".fasta"]
        else:
            raise ValueError

        genotype_files = [f for f in genotype_path.glob(fselect)]
        resistance_files = [f for f in genotype_path.glob("*.json")]

        if len(genotype_files) >= 1:
            genotypes = pandas.concat(
                [self.mp.parse_sccion(file=f) for f in genotype_files]
            ).reset_index(drop=True)

            for fs in fstrip:
                genotypes['name'] = genotypes.name.str.replace(fs, "")

            try:
                names_sorted = sorted(
                    genotypes.name.tolist(), key=lambda x: int(
                        re.findall(name_index_regex, x)[0]
                    )
                )
                genotypes = genotypes.set_index('name')\
                    .reindex(names_sorted).reset_index(level=0)
            except IndexError:
                self.logger.info(
                    'Could not extract index using regex from name string'
                )
                genotypes = genotypes.sort_values('name').reset_index(drop=True)

            genotype_name = self.outdir / f'{component}_genotypes.tsv'

            if exclude is not None:
                genotypes = genotypes.drop(columns=exclude)

            genotypes.to_csv(genotype_name, sep='\t', index=False)

        else:
            self.logger.info(
                f'Could not detect genotypes at: {self.path}'
            )
            genotypes = None

        if len(resistance_files) >= 1:
            resistance_genotypes = pandas.concat(
                [
                    self.mp.parse_mykrobe_genotype(file=f)
                    for f in resistance_files
                ]
            )
            resistance_genotypes['name'] = \
                resistance_genotypes.name.str.replace(".illumina.json", "")

            try:
                names_sorted = sorted(
                    resistance_genotypes.name.tolist(), key=lambda x: int(
                        re.findall(name_index_regex, x)[0]
                    )
                )
                resistance_genotypes = resistance_genotypes.set_index('name')\
                    .reindex(names_sorted).reset_index(level=0)
            except TypeError:
                self.logger.info(
                    'Could not extract index using regex from resistance_genotypes'
                    ' name string. Skipping sorting.'
                )

            resistance_phenotypes = pandas.concat(
                [
                    self.mp.parse_mykrobe_susceptibility(file=f)
                    for f in resistance_files
                ]
            )
            resistance_phenotypes['name'] = \
                resistance_phenotypes.name.str.replace(".illumina.json", "")
            try:
                names_sorted = sorted(
                    resistance_phenotypes.name.tolist(), key=lambda x: int(
                        re.findall(name_index_regex, x)[0]
                    )
                )
                resistance_phenotypes = resistance_phenotypes.set_index('name') \
                    .reindex(names_sorted).reset_index(level=0)
            except TypeError:
                self.logger.info(
                    'Could not extract index using regex from resistance_phenotypes'
                    ' name string. Skipping sorting.'
                )

            resistance_genotypes_name = \
                self.outdir / f'{component}_resistance_genotypes.tsv'
            resistance_genotypes.to_csv(
                resistance_genotypes_name, sep='\t', index=False, na_rep='-'
            )

            resistance_phenotypes_name = \
                self.outdir / f'{component}_resistance_phenotypes.tsv'
            resistance_phenotypes.to_csv(
                resistance_phenotypes_name, sep='\t', index=False
            )
        else:
            self.logger.info(
                f'Could not detect Illumina resistance predictions '
                f'at: {self.path}'
            )

        return genotypes

    def plot_genotype_heatmap(
        self,
        reference: pandas.DataFrame,
        genotypes: {str: pandas.DataFrame},
        common_isolates: bool = True,
        exclude: list = None
    ):
        """ Plot lower and upper diagonal of genotype heatmaps """

        combined = {}
        for workflow, genotype in genotypes.items():
            # Drop plasmid comparison, not reliable
            if 'plasmid' in genotype.columns:
                genotype.drop(columns='plasmid', inplace=True)
            if 'plasmid' in reference.columns:
                reference.drop(columns='plasmid', inplace=True)

            # Only common names in both reference and genotype:
            common = genotype.merge(reference, on="name").name
            if common_isolates:
                genotype = genotype[genotype['name'].isin(common)]\
                    .reset_index(drop=True)
                reference = reference[reference['name'].isin(common)]\
                    .reset_index(drop=True)

            if not len(genotype) == len(reference):
                raise ValueError

            # Make sure all entries in columns are sorted
            for col in genotype.columns:
                interim = [
                    c.split(';') for c in genotype[col]
                ]
                sorted_entries = [
                    ";".join(sorted([e.strip() for e in tu])) for tu in interim
                ]
                genotype[col] = sorted_entries

            for col in reference.columns:
                interim = [
                    c.split(';') for c in reference[col]
                ]

                sorted_entries = [
                    ";".join(sorted([e.strip() for e in tu])) for tu in interim
                ]
                reference[col] = sorted_entries

            if exclude:
                reference = reference.loc[~reference['name'].isin(exclude), :].reset_index(drop=True)
                genotype = genotype.loc[~genotype['name'].isin(exclude), :].reset_index(drop=True)

            g = genotype.drop(columns='name')
            r = reference.drop(columns='name')

            match = g.eq(r)
            match_means = match.mean(axis=0)

            combined[workflow] = match_means

        df = pandas.DataFrame(combined)

        fig, ax = plt.subplots(
            nrows=1, ncols=1, figsize=(
                1 * 7, 1 * 4.5
            )
        )

        fig.subplots_adjust(hspace=0.8)

        sns.heatmap(
            df, linewidths=.5, cmap="Greens", ax=ax, annot=True, vmin=0, vmax=1
        )

        # plt.tight_layout()
        fig.savefig(self.outdir / 'genotype_reference_heatmap.png')

    def collect_statistics(self, mode: str = '.filtered') -> pandas.DataFrame or None:

        path = self.path / 'ont' / 'qc'
        files = path.glob(f"*{mode}.stats.txt")

        nanoq = []
        for f in files:
            iid = f.name.rstrip(f'{mode}.stats.txt')
            df = self.mp.parse_nanoq(file=f)
            df['name'] = iid
            nanoq.append(df)

        if nanoq:
            nanoq = pandas.concat(nanoq)
            nanoq = self.str_num_sort_col(nanoq, col='name')
        else:
            self.logger.info('Could not detect nanoq summary files')
            return

        files = path.glob(f"*.coverage.txt")

        coverm = []
        for f in files:
            iid = f.name.rstrip(f'.coverage.txt')
            df = self.mp.parse_coverm(file=f)
            df['name'] = iid
            coverm.append(df)

        if coverm:
            coverm = pandas.concat(coverm)
            coverm = self.str_num_sort_col(coverm, col='name')
            coverm = coverm.drop(
                columns=['genome', 'placeholder1', 'placeholder2']
            )
        else:
            self.logger.info('Could not detect coverage files')
            return

        return nanoq.merge(coverm, on='name', how='inner')


class SignalPipeline(PoreLogger):

    def __init__(self, workdir: Path = Path("/tmp/.signal_work")):

        PoreLogger.__init__(self, level=logging.INFO, name="SignalPipeline")

        self.workdir = workdir
        self.workdir.mkdir(parents=True, exist_ok=True)

    def clean(self):

         shutil.rmtree(self.workdir)



class PathogenPipeline(PoreLogger):

    def __init__(self, workdir: Path = Path("/tmp/.pathogen_work")):

        PoreLogger.__init__(self, level=logging.INFO, name="PathogenPipeline")

        self.database_reports = dict()
        self.database_reads = dict()

        self.workdir = workdir
        self.workdir.mkdir(parents=True, exist_ok=True)

    def clean(self):

         shutil.rmtree(self.workdir)

    def collect_results(
        self, path: Path, groupby: str = None, level: str = "S"
    ) -> (dict, dict):

        if not path.exists():
            self.logger.info(
                f'Result path for pathogen pipeline {path} is not populated yet'
            )
            return dict(), dict()

        databases = [
            db for db in path.glob("*") if db.is_dir()
        ]

        if not databases:
            self.logger.info('No database results detected yet')

        database_reports = {}
        database_reads = {}
        for db in databases:
            reports, total_reads = self.get_bracken_results(
                path=db, groupby=groupby, level=level
            )
            database_reports[db.name] = reports
            database_reads[db.name] = total_reads

        self.database_reports = database_reports
        self.database_reads = database_reads

        # For now just the Bracken results
        return database_reports, database_reads

    def get_server_data(
        self, host_label: str = "Homo sapiens",
        major_threshold: float = 0.5, minor_threshold: float = 0.1
    ) -> dict:

        # filtered for server data
        server_data = {}
        for db, reports in self.database_reports.items():
            server_data[db] = {}
            for g, r in reports.items():

                if minor_threshold < 1:
                    minor_limit_col = 'percent'
                    minor_threshold = float(minor_threshold)
                else:
                    minor_limit_col = 'reads'
                    minor_threshold = int(minor_threshold)

                if minor_threshold < 1:
                    major_limit_col = 'percent'
                    major_threshold = float(major_threshold)
                else:
                    major_limit_col = 'reads'
                    major_threshold = int(major_threshold)

                # Host
                if host_label in r.tax.tolist():
                    host_data = r[r['tax'] == host_label]
                else:
                    host_data = pandas.DataFrame()

                # Percent decimals have been increased in merge report function
                # so this should work for NovaSeq / Illumina
                major_taxa = r[
                    (r[major_limit_col] >= major_threshold) &
                    (r[minor_limit_col] >= minor_threshold) &
                    (r['tax'] != host_label)
                    ]
                minor_taxa = r[
                    (r[major_limit_col] < major_threshold) &
                    (r[minor_limit_col] >= minor_threshold) &
                    (r['tax'] != host_label)
                    ]

                minor_taxa = pandas.concat(
                    (minor_taxa, host_data)
                )

                server_data[db][g] = {
                    'major': self.report_to_server_data(major_taxa),
                    'minor': self.report_to_server_data(minor_taxa)
                }

        return server_data

    @staticmethod
    def report_to_server_data(report: pandas.DataFrame):

        reads = report.reads.tolist()
        percent = report.percent.tolist()
        labels = [
            label.strip() for label in report.tax.tolist()
        ]

        return [
            {
                'text': l,
                'values': [r],
                'percent': p,
                'background-color': "#737373",
                'detached': False,
                'assembly': [
                    {
                        'text': 'ctg_202',
                        'values': [202023],
                        'coverage': {
                            'windows': [0, 0, 0, 1, 4, 6, 3, 0, 0, 0],
                            'mean': 1.5
                        },
                    },
                    {
                        'text': 'ctg_202',
                        'values': [102736],
                        'coverage': {
                            'windows': [0, 0, 0, 1, 4, 6, 3, 0, 0, 0],
                            'mean': 1.5
                        }
                    }

                ]
            }
            for r, p, l in zip(reads, percent, labels)
        ]

    def get_bracken_results(
        self, path: Path, groupby: str = None, level: str = None
    ):

        # Bracken reports are the primary result for now
        bracken_files = list(
            path.glob("*.bracken.report")
        )

        # If a Bracken file is present a Kraken file must be present:
        kraken_files = list(
            path.glob("*.kraken.report")
        )

        # Groupby e.g. barcode regex, otherwise assume each file = sample

        bracken_reports = {}  # group id, pandas df
        reads = {}

        bracken_groups = self.group_files_by(
            files=bracken_files, regex=groupby, fext=".bracken.report"
        )
        kraken_groups = self.group_files_by(
            files=kraken_files, regex=groupby, fext=".kraken.report"
        )

        print(bracken_groups)
        print(kraken_groups)

        if set(bracken_groups.keys()) != set(kraken_groups.keys()):
            self.logger.error(
                'Bracken and Kraken files do not match. '
                'Something has gone wrong in the pipeline.'
            )

        # Merge reports of grouped:
        for group, bracken_files in bracken_groups.items():

            self.logger.debug(f'Merging reports from: {group}')
            # Merge Bracken reports, no output, single column
            bracken = merge_reports(
                r_files=bracken_files, out=None, c_only=True
            )
            # Merge Kraken reports to get accurate read counts

            if bracken is None:
                self.logger.error(
                    f'Failed merging of report files: {bracken_files}'
                )

            bracken_reports[group] = bracken

            reads[group] = self.get_merged_kraken_read_counts(
                report_files=kraken_groups[group]
            )

        if level:
            _level_reports = {}
            for group, df in bracken_reports.items():
                _level_reports[group] = df[
                    df['level'] == level
                    ].sort_values('percent', ascending=False)
            bracken_reports = _level_reports

        return bracken_reports, reads

    @staticmethod
    def get_merged_kraken_read_counts(
            report_files: [Path]
    ):

        unclassified = 0
        classified = 0
        total = 0
        for file in report_files:
            report = pandas.read_csv(
                file, sep='\t', comment='#', header=None,
                names=["percent", "reads", "direct", "level", "taxid", "tax"]
            )
            report.tax = [s.strip() for s in report['tax'].tolist()]

            u = int(report.at[0, 'reads'])
            root = int(report.at[1, 'reads'])

            unclassified += u
            classified += root
            total += root+u

        return {
            'unclassified': unclassified,
            'classified': classified,
            'total': total
        }

    def group_files_by(
        self, files: [Path],
        regex: str = r"barcode\d\d",
        fext: str = '.bracken.report'
    ) -> dict:

        grouped_files = dict()
        for file in files:
            if regex:
                try:
                    found = re.search(regex, file.name).group()  # Only expect one
                    if found not in grouped_files.keys():
                        grouped_files[found] = [file]
                    else:
                        grouped_files[found].append(file)
                except AttributeError:
                    self.logger.debug(
                        f"Could not find '{regex}' match in {file.name}"
                    )
                    pass
            else:
                # Bracken extension
                self.logger.info(
                    "No regex provided. Grouping by report filename."
                )
                file_name = file.name.replace(fext, "")
                grouped_files[file_name] = [file]

        return grouped_files

    def assemble(
        self,
        grouped_files: dict,
        assembler="raven",
        assembler_cpu: int = 2,
        threads: int = 2,
        raven_args: str = ""
    ):

        """ Multi-threaded Raven assembly of grouped files """

        assembly_methods = dict(
            raven=self._assemble_raven
        )

        with mp.Pool(processes=threads) as pool:
            pool.map_async(
                func=assembly_methods[assembler],
                iterable=[
                    list(data)+[assembler_cpu, raven_args]
                    for data in grouped_files.items()
                ]
            )

    def _assemble_raven(self, group_data: list):

        group, files, cpu, raven_args = group_data

        self.logger.info(
            f'Starting Raven assembly of group: {group} (n = {len(files)})'
        )
        cat_files = " ".join([str(f) for f in files])  # pathlib

        gfa_file = Path(f"{self.workdir / f'_{group}.gfa'}")
        fa_file = Path(f"{self.workdir / f'{group}.fa'}")
        fq_file = Path(f"{self.workdir / f'_{group}.fq'}")

        cat_string = f"cat {cat_files} > {fq_file}"
        raven_string = f"raven -t {cpu} {raven_args} --graphical-fragment-assembly {gfa_file} {fq_file}"
        awk_string = r"""awk '$1 ~/S/ {print ">"$2"\n"$3}'"""
        awk_string += f" {gfa_file} > {fa_file}"
        rm_string = f"rm {fq_file} {gfa_file} raven.cereal"

        cmd = " && ".join([cat_string, raven_string, awk_string, rm_string])

        try:
            run_cmd(cmd=cmd, shell=True)
        except:
            # always clean up
            if fq_file.exists():
                fq_file.unlink()
            if gfa_file.exists():
                gfa_file.unlink()

        n_contigs = 0
        bp = 0
        if fa_file.exists():
            for line in fa_file.open().read():
                if line.startswith('>'):
                    n_contigs += 1
                else:
                    bp += len(line.strip())

        msg = f'Assembled group: {group} [ contigs = {n_contigs}, bp = {bp} ]'

        self.logger.info(msg)

        return group, fa_file, n_contigs, bp


class MasterParser:

    def __init__(self):

        self.mlst_columns = [
            "file", "species", "sequence_type",
        ]
        self.mash_columns = [
            "file", "ref", "dist", "p-value", "match",
        ]
        self.kraken_report_columns = [
            "percent", "reads", "direct", "level", "taxid", "taxonomy",
        ]
        self.prokka_columns = [
            "locus_tag", "feature", "length", "gene", "ec", "cog", "product",
        ]

        self.sccion_columns = [
            'name', 'species', 'mlst', 'meca', 'pvl', 'spa', 'scc', 'resistance', 'plasmid'
        ]

        self.nanoq_columns = [
            'reads', 'bp', 'longest', 'shortest', 'mean_length',
            'median_length', 'mean_qscore', 'median_qscore'
        ]

        self.coverm_columns = [
            'genome', 'mean_coverage', 'placeholder1', 'placeholder2'
        ]

        self.abricate_columns = [
            "file", "sequence", "start", "end", "gene", "coverage_bases",
            "coverage_map", "gaps",  "coverage", "identity", "database",
            "accession", "product",
        ]
        self.kleborate_columns = [
            "strain", "species", "st", "virulence_score", "resistance_score",
            "Yersiniabactin", "YbST", "Colibactin", "CbST", "Aerobactin",
            "AbST", "Salmochelin", "SmST", "rmpA", "rmpA2", "wzi",  "K_locus",
            "K_locus_confidence", "O_locus", "O_locus_confidence", "AGly",
            "Col", "Fcyn", "Flq", "Gly", "MLS", "Ntmdz", "Phe", "Rif", "Sul",
            "Tet", "Tmt", "Bla", "Bla_Carb", "Bla_ESBL", "Bla_ESBL_inhR",
            "Bla_broad", "Bla_broad_inhR",
        ]

    # Public helper methods for result file parsing

    def parse_kraken(self, file) -> pandas.DataFrame:
        """Parse single report output from Kraken"""
        return pandas.read_csv(
            file, header=None, sep="\t", names=self.kraken_report_columns
        )

    def parse_mash(self, file) -> pandas.DataFrame:
        """Parse single output from Mash"""
        return pandas.read_csv(
            file, sep="\t", header=None, names=self.mash_columns
        )

    def parse_abricate(self, file) -> pandas.DataFrame:
        """Parse single output from Abricate"""
        return pandas.read_csv(
            file, sep="\t", header=0, names=self.abricate_columns
        )

    def parse_prokka(self, file, gff=True) -> pandas.DataFrame:
        """Parse single output feature table from Prokka"""

        if gff:
            return self._parse_prokka_gff(file)
        else:
            return pandas.read_csv(
                file, sep="\t", header=0, names=self.prokka_columns
            )

    def parse_sccion(self, file) -> pandas.DataFrame:

        return pandas.read_csv(
            file, sep="\t", header=None, names=self.sccion_columns
        )

    def parse_nanoq(self, file) -> pandas.DataFrame:

        return pandas.read_csv(
            file, sep=" ", header=None, names=self.nanoq_columns
        )

    def parse_coverm(self, file) -> pandas.DataFrame:

        return pandas.read_csv(
            file, sep='\t', header=None, skiprows=1, names=self.coverm_columns
        )

    @staticmethod
    def parse_dnadiff(file: Path, name_strip: str = None) -> pandas.DataFrame:

        contigs_query, contigs_ref = None, None
        snps, indels, aid = None, None, None
        one_to_one = False
        for line in file.open('r'):
            if line.startswith("1-to-1"):
                one_to_one = True
            if line.startswith('TotalSeqs'):
                contigs_query = [e for e in line.strip().split() if e][-1]
                contigs_ref = [e for e in line.strip().split() if e][-2]
            elif line.startswith('TotalSNPs'):
                snps = [e for e in line.strip().split() if e][-1]
            elif line.startswith('TotalIndels'):
                indels = [e for e in line.strip().split() if e][-1]
            elif line.startswith('AvgIdentity') and one_to_one:
                aid = [e for e in line.strip().split() if e][-1]
                one_to_one = False

        name = file.stem
        if name_strip:
            name = name.strip(name_strip)

        try:
            qual = 10*math.log10(100/(100-float(aid)))
        except ZeroDivisionError:
            qual = 70

        return pandas.DataFrame(
            data=[{
                'contigs_query': int(contigs_query),
                'contigs_ref': int(contigs_ref),
                'snps': int(snps),
                'indels': int(indels),
                'name': name,
                'identity': aid,
                'quality': float(qual)

            }],
        )

    @staticmethod
    def _parse_prokka_gff(file):

        data = []
        with open(file, 'r') as infile:
            for line in infile:
                if line.startswith('##'):
                    pass
                elif line.startswith('##FASTA'):
                    return 0
                else:
                    data.append(
                        line.strip().split('\t')
                    )

        return pandas.DataFrame(data, columns=[
            'contig', 'inference', 'type', 'start', 'stop', 'unknown_1',
            'strand', 'unknown_2', 'description'
        ])

    def parse_mlst(self, file) -> pandas.DataFrame:
        """Parse single output from mlst (github.com/tseemann/mlst)"""

        df = pandas.read_csv(file, sep="\t", header=None)

        df.columns = self.mlst_columns + [
            str(i) for i, col in enumerate(df.columns, 1)
        ][:-3]

        return df

    def parse_kleborate(self, file) -> pandas.DataFrame:
        """ Parse single file output from Kleborate """

        return pandas.read_csv(
            file, sep='\t', header=0
        )

    def parse_mykrobe_susceptibility(self, file) -> pandas.DataFrame:

        susceptibilities, _ = self._parse_mykrobe(file)

        return susceptibilities

    def parse_mykrobe_genotype(self, file) -> pandas.DataFrame:

        _, genotype = self._parse_mykrobe(file)

        return genotype

    def parse_mykrobe_lineage(self, file) -> pandas.DataFrame:

        return self._parse_mykrobe(file, lineage=True)

    @staticmethod
    def _parse_mykrobe(
        file, lineage=False
    ) -> (pandas.DataFrame, pandas.DataFrame) or pandas.DataFrame:
        """Parse single output from MykrobePredictor, only parse
        called resistances and the genes that called them.
        """

        with open(file, "r") as infile:
            mykrobe = json.load(infile)

            if lineage:
                # TODO make sure order highest if multiple?
                keys = [
                    k for k, v in list(
                        get_subdict("lineage", mykrobe)
                    )[0].items()
                ]
                lineage_data = pandas.DataFrame.from_dict({'lineage': keys})
                lineage_data['name'] = file.name
                return lineage_data
            else:
                susceptibility = list(
                    get_subdict("susceptibility", mykrobe)
                )[0]

                data = {
                    drug: {
                        "susceptibility": data["predict"],
                        "genotype": list(
                            data["called_by"].keys()
                        )
                        if "called_by" in data.keys() else None
                    }
                    for drug, data in susceptibility.items()
                }

                susceptibilities = pandas.DataFrame.from_dict(
                    {
                        drug: list(dat["susceptibility"])
                        for drug, dat in data.items()
                    }
                )

                genotypes = pandas.DataFrame.from_dict(
                    {
                        drug: [None] if not dat["genotype"]
                        else [",".join(dat["genotype"])]
                        for drug, dat in data.items()
                     }
                )
                genotypes['name'] = file.name
                susceptibilities['name'] = file.name

                return susceptibilities, genotypes
