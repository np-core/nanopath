import uuid
import pandas
import shlex
import subprocess
import urllib.request
import requests

from io import StringIO
from ratelimiter import RateLimiter

from tqdm import tqdm
from io import StringIO
from pathlib import Path
from pandas.errors import EmptyDataError
from nanopath.pathfinder import get_aspera_key
from nanopath.pathfinder import get_genome_sizes


class MiniAspera:

    def __init__(self, force=False):

        self.port = 33001
        self.limit = 1024

        self.force = force

        self.ascp = 'ascp'
        self.key = get_aspera_key()

        self.fasp = "era-fasp@fasp.sra.ebi.ac.uk:"

    def download_batch(
        self,
        data: pandas.DataFrame,
        outdir: Path = Path.cwd(),
        ftp: bool = False,
    ):

        outdir.mkdir(exist_ok=True, parents=True)

        with tqdm(total=len(data)) as pbar:
            pbar.set_description("Downloading batch")

            for i, fastq in data.iterrows():

                if ftp:
                    fq1_address = fastq["ftp_1"]
                else:
                    fq1_address = fastq["ftp_1"].replace(
                        "ftp.sra.ebi.ac.uk", self.fasp
                    )

                fq1_path = outdir / Path(fq1_address).name

                self.download(
                    address=fq1_address,
                    outfile=fq1_path,
                    force=self.force,
                    ftp=ftp
                )

                if fastq["ftp_2"] and not isinstance(
                    fastq["ftp_2"], float
                ):
                    if ftp:
                        fq2_address = fastq["ftp_2"]
                    else:
                        fq2_address = fastq["ftp_2"].replace(
                            "ftp.sra.ebi.ac.uk", self.fasp
                        )
                    fq2_path = outdir / Path(fq2_address).name
                    self.download(
                        address=fq2_address,
                        outfile=fq2_path,
                        force=self.force,
                        ftp=ftp
                    )

                pbar.update(1)

    def download_raw(
        self,
        data: pandas.DataFrame,
        submitted: bool = False,
        outdir: Path = Path.cwd(),
        ftp: bool = False,
    ):

        outdir.mkdir(exist_ok=True, parents=True)

        with tqdm(total=len(data)) as pbar:
            pbar.set_description("Downloading from raw table")

            for i, row in data.iterrows():

                if submitted:
                    ftp_data = row['submitted_ftp']
                else:
                    ftp_data = row['fastq_ftp']

                files = ftp_data.split(";")

                for fq in files:
                    if ftp:
                        fq_address = fq
                    else:
                        fq_address = fq.replace(
                            "ftp.sra.ebi.ac.uk", self.fasp
                        )

                    fq_path = outdir / Path(fq_address).name

                    self.download(
                        address=fq_address,
                        outfile=fq_path,
                        force=self.force,
                        ftp=ftp
                    )

                pbar.update(1)

    def download(self, address, outfile, force=False, ftp=False):

        # Skip existing files
        if not force and outfile.exists():
            print(f"File exists: {outfile}")
            return

        if ftp:
            cmd = shlex.split(f"wget {address} -O {outfile}")
        else:
            cmd = shlex.split(
                f"{self.ascp} -QT -l {str(self.limit)}m"
                f" -P{str(self.port)} -i {self.key} -q "
                f"{address} {outfile}"
            )
        try:
            subprocess.call(cmd)
        except subprocess.CalledProcessError:
            print("Error in subprocess.")
            raise  # handle errors in the called executable
        except OSError:
            print("Executable not found.")
            raise  # executable not found

    @staticmethod
    def read_batch(file):

        return pandas.read_csv(file, sep='\t', header=0)


class Survey:

    """
    Simple accessor class to pull short-read data from the ENA.
    """

    def __init__(self, outdir: Path = None, fields: str = "all"):

        self.outdir = outdir

        if fields == "all":
            self.fields = "all"
        else:
            self.fields = "study_accession,sample_accession,run_accession,tax_id,scientific_name," \
                "fastq_ftp,fastq_bytes,submitted_ftp,submitted_bytes," \
                "read_count,base_count," \
                "instrument_platform,instrument_model," \
                "library_layout,library_source,library_strategy," \
                "location,country,collection_date"

        self.url = f'https://www.ebi.ac.uk/ena/portal/api/search?result=read_run' \
                   f'&fields={self.fields}&query="'

        self.term_chunk_size = 200  # REST API limited by 6000 characters

        self.results = dict()
        self.query = pandas.DataFrame()

    @staticmethod
    def write_query_file(csv_file="query.csv", query_results=None):

        query_results.to_csv(csv_file)

    def read_query_file(self, file):

        self.query = pandas.read_csv(file, sep='\t', header=0)

        return self.query

    def parse_biosample(self):

        pass

    def display(self):

        pass

    def filter_query(self, ops: str):
        """ Filter the query dataframe by
        applying pandas query (boolean
        filter operations) string """
        self.query = self.query.query(ops)

    @staticmethod
    def batch_output(batches, outdir="batches", exist_ok=True):

        outdir = Path.cwd() / outdir
        for i, batch in enumerate(batches):
            batch_dir = outdir / f"batch_{i}"
            batch_dir.mkdir(parents=True, exist_ok=exist_ok)

            yield batch_dir, batch

    @staticmethod
    def batch(query_result: pandas.DataFrame, batch_size=None, max_gb=None):

        if max_gb:
            running_size = []
            batches = []
            start = 0
            for i, entry in query_result.iterrows():
                running_size.append(float(entry["size"])/1000)
                if sum(running_size) > max_gb:
                    batches += [query_result[start:i]]
                    start = i
                    running_size = []
            return batches
        elif batch_size:
            return [
                query_result[i:i + batch_size] for i in range(
                    0, query_result.shape[0], batch_size
                )
            ]
        else:
            raise ValueError(
                "Either maximum gigabytes or batch size must be set."
            )

    def query_ena(
        self, species="Staphylococcus aureus", scheme="illumina",
        study: str = None, sample: list or str = None, term: str = None,
        submitted_fastq: bool = False, allow_merged_fastq: bool = False,
        allow_missing: bool = False,
    ) -> dict:

        """
        Search the ENA warehouse for raw sequence reads with the
        following default parameters to conform to current implementations
        of WGS analysis pipelines.
        """

        # Format queries correctly:
        if isinstance(scheme, str):
            scheme = scheme.lower()

        # At the moment, restrict queries to schemes:
        if scheme == "illumina":

            platform = "ILLUMINA"
            source = "GENOMIC"
            layout = "PAIRED"
            strategy = "WGS"

        elif scheme == "nanopore":

            platform = "OXFORD_NANOPORE"
            source = "GENOMIC"
            layout = "SINGLE"
            strategy = "WGS"

        else:
            platform = ""
            source = ""
            layout = ""
            strategy = ""

        if species:
            terms = [self._construct_species_query(
                species, platform, source, layout, strategy
            )]
        elif study:
            terms = [self._construct_study_query(study)]
        elif term:
            terms = [term]
        elif sample:
            print(f'Splitting query into chunks of {self.term_chunk_size}')
            sample_chunks = self.chunks(sample, self.term_chunk_size)
            terms = [self._construct_sample_query(sc) for sc in sample_chunks]
        else:
            raise ValueError("Need to specify either species, study accession, or custom search term for Survey.")

        queries = {}
        for i, t in enumerate(terms):
            print(f'Submitting query: {i}')
            url = self.url + t
            url = url + '"'
            url = url.replace(" ", "%20")
            print(url)
            df = self._query(url)

            query_results = self._sanitize_ena_query(
                df, url, submitted_fastq=submitted_fastq,
                allow_merged_fastq=allow_merged_fastq,
                allow_missing=allow_missing
            )

            queries[i] = {'term': t, 'url': url, 'results': query_results}

        return queries

    @staticmethod
    def _query(url) -> pandas.DataFrame:

        query_results = StringIO(urllib.request.urlopen(url)
                                 .read().decode('utf-8'))

        return pandas.read_csv(query_results, sep="\t")

    @staticmethod
    def _sanitize_ena_query(
        df, url, submitted_fastq, allow_merged_fastq, allow_missing
    ) -> pandas.DataFrame:

        if not allow_missing:
            # Drop rows with missing FTP links:
            df = df.dropna(subset=["fastq_ftp"])

        if df.empty:
            raise ValueError(
                f"Query results are empty, check your "
                f"query string or confirm search manually "
                f"at ENA: {url}"
            )

        try:
            genome_sizes = get_genome_sizes()
        except FileNotFoundError:
            genome_sizes = None

        sanitized_dict = {}
        # Iterate over rows of dataframe to sanitize query entries:
        for index, entry in df.iterrows():

            # FTP Links
            if not submitted_fastq:
                ftp_links = str(entry["fastq_ftp"]).strip(";").split(";")
                ftp_sizes = str(entry["fastq_bytes"]).strip(";").split(";")
            else:
                ftp_links = str(entry["submitted_ftp"]).strip(";").split(";")
                ftp_sizes = str(entry["submitted_bytes"]).strip(";").split(";")

            if entry["library_layout"] == "PAIRED":

                # TODO: split FASTQ into forward and reverse after download
                if len(ftp_links) < 2 and allow_merged_fastq:
                    ftp_1, ftp_2 = ftp_links[0], None
                else:
                    # Get last two links and file sizes,
                    # which should be forward + reverse,
                    # check that they are conforming to
                    # pattern _1.fastq.gz and _2.fastq.gz:
                    ftp_links = ftp_links[-2:]
                    ftp_1, ftp_2 = ftp_links
                    ftp_sizes = ftp_sizes[-2:]

            elif entry["library_layout"] == "SINGLE":

                # If there are fewer than one or more
                # than two links, ignore accession:
                if len(ftp_links) != 1:
                    continue
                else:
                    ftp_1 = ftp_links[0]
                    ftp_2 = None
            else:
                raise ValueError("Layout must be either SINGLE or PAIRED")

            # Convert to MB:
            try:
                size = sum([int(byte)/1024/1024 for byte in ftp_sizes])
            except ValueError:
                size = None

            try:
                reads = int(entry["read_count"])
            except ValueError:
                reads = None

            try:
                bases = int(entry["base_count"])
            except ValueError:
                bases = None

            if bases:
                if genome_sizes is None:
                    coverage = None
                else:
                    try:
                        coverage = bases/(float(
                            genome_sizes.loc[entry["tax_id"], "size"]
                        )*1000000)
                    except KeyError or ValueError or ZeroDivisionError:
                        coverage = None
            else:
                coverage = None

            name = Path(ftp_1).stem.split('_')

            try:
                name = '_'.join(name[:-1])  # to get rid of _1.fq.gz for example
            except IndexError:
                name = None

            entry_dict = {
                "id": str(uuid.uuid4()),
                "accession": entry['run_accession'],
                "name": name,
                "ftp_1": ftp_1,
                "ftp_2": ftp_2,
                "size": size,
                "reads": reads,
                "bases": bases,
                "coverage": coverage,
                "layout": entry["library_layout"],
                "platform": entry["instrument_platform"],
                "model": entry["instrument_model"],
                "source": entry["library_source"],
                "strategy": entry["library_strategy"],
                "tax_id": entry["tax_id"],
                "sample": entry["sample_accession"],
                "location": entry["location"],
                "country": entry["country"],
                "collection_date": entry["collection_date"],
                "scientific_name": entry['scientific_name']
            }

            sanitized_dict[
                entry["run_accession"]
            ] = entry_dict

        try:
            df = pandas.DataFrame(sanitized_dict).T
        except EmptyDataError:
            raise ValueError(f"No results were returned for query: {url}")

        return df

    @staticmethod
    def chunks(lst, n):
        """Yield successive n-sized chunks from lst."""
        for i in range(0, len(lst), n):
            yield lst[i:i + n]

    @staticmethod
    def _construct_species_query(
        species, platform, source, layout, strategy
    ) -> str:
        """ Construct a specific ENA query string to obtain paired-end
        Illumina reads or Nanopore reads from the DataWarehouse. """

        q = ''
        if species:
            q += f'tax_name("{species}")'
        if platform:
            q += f' AND instrument_platform={platform}'
        if platform:
            q += f' AND library_layout={layout}'
        if source:
            q += f' AND library_source={source}'
        if strategy:
            q += f' AND library_strategy={strategy}'

        return q

    @staticmethod
    def _construct_study_query(study: list or str):
        """Construct a specific ENA query
        string to obtain reads from a Study Accession"""

        if isinstance(study, str):
            return f'"study_accession={study}"'
        else:
            return "%20OR%20".join(
                f'study_accession={s}' for s in study
            )

    @staticmethod
    def _construct_sample_query(sample: list or str):
        if isinstance(sample, str):
            return f'"run_accession={sample}"'
        else:
            return "%20OR%20".join(
                f'run_accession={s}' for s in sample
            )


class BioSampler:

    """ Class to query BioSample at ENA """

    def __init__(self):

        self.url = "https://www.ebi.ac.uk/ena/portal/api/search?result=read_run" \
                   "&fields=sample_accession,location,country,collection_date,scientific_name&query="

    def process_list(self, file: Path, query_column: str, sep: str, output: Path):

        df = self.read(file, sep=sep)
        data_queries = [f"run_accession={run_accession}" for run_accession in df[query_column]]

        data_query_chunks = self.chunks(data_queries, 200)  # < max length of query string!

        query_data = []
        for query_chunks in data_query_chunks:
            query_string = " OR ".join(query_chunks)
            response_text = self.query_url(query=query_string)
            df = pandas.read_csv(
                StringIO(response_text), sep='\t', header=0
            )
            print(df)
            query_data.append(df)

        df = pandas.concat(query_data)
        df.to_csv(output, sep="\t", header=True, index=False)

    @staticmethod
    def chunks(lst, n):
        """Yield successive n-sized chunks from lst."""
        for i in range(0, len(lst), n):
            yield lst[i:i + n]

    @staticmethod
    def read(file: Path, sep="\t") -> pandas.DataFrame:

        """ Read table of data and extract the column containing BioSample IDs """

        return pandas.read_csv(file, sep=sep, header=0)

    @RateLimiter(max_calls=3, period=1)
    def query_json(self, query: str):

        print(f'Query: {query}')
        try:
            response = requests.get(
                self.url + f"{query}", headers={"Accept": "application/json"}
            )
            response.raise_for_status()
        except requests.exceptions.HTTPError as http_err:
            print(f'HTTP error occurred: {http_err}')  # Python 3.6
        except Exception as err:
            print(f'Other error occurred: {err}')  # Python 3.6
        else:
            print(f'Success: {query}')

            return response.json()

    @RateLimiter(max_calls=3, period=1)
    def query_url(self, query: str = ""):

        try:
            url = self.url + f"{query}"
            url = url.replace(" ", "%20")
            response = requests.get(url)
            response.raise_for_status()
        except requests.exceptions.HTTPError as http_err:
            print(f'HTTP error occurred: {http_err}')  # Python 3.6
        except Exception as err:
            print(f'Other error occurred: {err}')  # Python 3.6
        else:
            print('Success!')
            return response.text
