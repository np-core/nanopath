""" Old code """

import os
import bs4
import json
import time
import pandas
import shutil
import requests
import textwrap
import urllib

from numpy import int64
from subprocess import call, CalledProcessError

from nanopath.utils import PoreLogger

from numpy import nan

from BCBio import GFF
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# Somewhere in report summary a SettingWithCopyWarning is thrown.
# Disabled warning for now, but need to fix.

pandas.options.mode.chained_assignment = None

FASTA_EXTENSIONS = (".fa", ".fasta")


class NCTC3000(PoreLogger):
    """
    Class for parsing HTML table of closed bacterial genomes from NCTC3000. Creates data table for local copy
    of all (or a subset by species) genome assemblies from NCTC3000. Generates .fasta and .gbk files, then
    analyses with mlst (https://github.com/tseemann/mlst) and abricate (https://github.com/tseemann/abricate) for
    virulence and antibiotic resistance prediction. Summarize output and update projects and species.
    """

    def __init__(self, path=os.getcwd(), species="Staphylococcus aureus", force=False, verbose=True):

        PoreLogger.__init__(self, name="NCTC")

        self.date = time.strftime("%d-%m-%Y-%H-%M")

        self.path = os.path.abspath(path)
        self.project = True

        self.species = species
        self.species_name = self.species.replace(" ", "_").lower()
        self.species_path = os.path.join(self.path, self.species_name)

        self.force = force
        self.verbose = verbose

        # Infrastructure

        self.config_path = os.path.join(self.path, ".nctc")

        # Analysis with NCTC Programs and Snakemake

        self.nctc = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))

        self.abricate = os.path.join(self.nctc, "bin", "abricate", "bin", "abricate")

        self.env_file = os.path.join(self.nctc, "env", "nctc.json")
        self.snake_file = os.path.join(self.nctc, "pipe", "nctc.snk")
        self.config_file = os.path.join(self.nctc, "pipe", "nctc.json")

        self.config_run_file = None

        # HTML NCTC3000

        self.url = "http://www.sanger.ac.uk/resources/downloads/bacteria/nctc/"

        self.tables = dict()

        self.files = dict()
        self.report_files = dict()

    def _load_project(self):

        files = [file for file in os.listdir(self.config_path) if file.endswith(".tab")]

        table_history = self._load_history(files)

        try:
            last = sorted(table_history.keys(), key=lambda x: time.strptime(x, "%d-%m-%Y-%H-%M"))[-1]
        except IndexError:
            raise ValueError("There is no history to load in", self.config_path)

        last_table = table_history[last]

        self.tables[last_table["name"]] = NCTCTable(
            name=last_table["name"], dataframe=last_table["dataframe"], links=last_table["links"]
        )

    def _load_history(self, files):

        table_history = dict()

        for file in files:
            components = file.replace(".tab", "").split("_")

            date = components[-2]
            table_name = components[0]
            table_type = components[-1]

            if date not in table_history.keys():
                table_history[date] = {"name": table_name, "created": date}

            if table_type == "links":
                table_history[date]["links"] = os.path.join(self.config_path, file)
            else:
                table_history[date]["dataframe"] = os.path.join(self.config_path, file)

        return table_history

    def _setup_project(self):

        """
        Called upon initialisation of class NCTC3000
        1. If self.path exists and contains path/.nctc, create a new project directory.
        2. If self.path exists and does not contain path/.nctc, it is a user directory.
        3. If self.path does not exist, create new project directory.
        """

        if os.path.exists(self.path) and os.path.exists(self.config_path):
            self.logger.info("Operating on project at", self.path)
            self._load_project()
        elif os.path.exists(self.path) and not os.path.exists(self.config_path):
            self.logger.info("Could not find project at", self.path)
            exit(1)
        elif not os.path.exists(self.path):
            self.logger.info("Creating project at", self.path)
            self._create_project()
        else:
            self.logger.info("Hmmmmm...")

    def make(self, strict=True, force=False):

        """
        Main function for creating a new project and maintaining it.
        """

        if self.verbose:
            self._print_parameters(locals())

        self._setup_project()

        self.parse_website()
        self.parse_species(strict=strict, force=force)

    def _create_project(self):

        gff_path = os.path.join(self.species_path, "gff")
        fasta_path = os.path.join(self.species_path, "fasta")
        genbank_path = os.path.join(self.species_path, "genbank")

        for path in (self.path, self.species_path, self.config_path,
                     gff_path, fasta_path, genbank_path):
            os.makedirs(path)

    def parse_website(self):

        self.logger.info("Parsing NCTC300 website at Sanger Institute")

        site = requests.get(self.url)
        html = bs4.BeautifulSoup(site.text, "lxml")

        for table in html.find_all('table'):
            nctc_table = NCTCTable(name="nctc_assemblies", table=table)

            self.tables[nctc_table.name] = nctc_table

            nctc_table.summarise()

            self.logger.info("Storing table", nctc_table.name, "from date", self.date, "at", self.config_path)

            nctc_table.store(self.config_path)

    def parse_species(self, name="genomes", assembly="manual", strict=True, force=False):

        self.logger.info("Fetching assemblies for", self.species, "from NCTC3000")

        bacteria = self.tables[name]

        if strict:
            species_indices = bacteria.dataframe[(bacteria.dataframe["Species"] == self.species) &
                                                 (bacteria.dataframe["Chromosomes"] > 0)].index
        else:
            species_indices = bacteria.dataframe[(bacteria.dataframe["Species"] == self.species)].index

        gff_files = download_assemblies(gff_path=os.path.join(self.species_path, "gff"),
                                        links=get_links(bacteria, species_indices, assembly), force=force)

        gff_to_fasta(gff_files, fasta_path=os.path.join(self.species_path, "fasta"),
                     genbank_path=os.path.join(self.species_path, "genbank"), verbose=self.verbose)

    def _print_parameters(self, parameters):

        """
        http://stackoverflow.com/questions/10724495/getting-all-arguments-and-values-passed-to-a-python-function
        :param parameters: within function call to local()
        :return:
        """

        parameters.pop("self")

        for key, value in sorted(parameters.items()):
            self.logger.info("{key} = {value}".format(key=key, value=value))

    def collect(self, out_path, user_path=None, merge="outer", resistance_db="resfinder", virulence_db="vfdb",
                cnv=True, cnv_mincov=95, single=False, csv=True, force=False):

        out_path = os.path.abspath(out_path)

        if self.verbose:
            self._print_parameters(locals())

        if os.path.exists(out_path) and force:
            shutil.rmtree(out_path)

        os.makedirs(out_path)

        _, _, mlst_path, res_path, vir_path = self._get_analysis_paths(resistance_db, virulence_db, mode="collect")

        frames = parse_reports(mlst_path, res_path, vir_path, resistance_db=resistance_db, virulence_db=virulence_db)

        cleaned_frames = clean_reports(frames)
        merged_frames = merge_frames(cleaned_frames, merge=merge)

        if cnv:
            cnv_frames = get_copy_numbers(cleaned_frames, cnv_mincov)
            cnv_merged = merge_frames(cnv_frames, left=cleaned_frames["mlst"], merge=merge)

            cleaned_frames.update(cnv_frames)
            merged_frames.update(cnv_merged)

        if single:
            write_frames = cleaned_frames
        else:
            write_frames = merged_frames

        for output, frame in write_frames.items():
            if csv:
                sep = ","
                file_path = os.path.join(out_path, output + ".csv")
            else:
                sep = "\t"
                file_path = os.path.join(out_path, output + ".tab")

            frame.to_csv(file_path, sep=sep, na_rep=".")
            self.logger.info("Written summary to file:", file_path)

    def type(self, resistance_db="resfinder", virulence_db="vfdb", minid=90, mincov=95,
             cluster=True, mlst=True, resistance=True, virulence=True, force=False):

        """
        Analyse a collection of species reference genomes using
            mlst: https://github.com/tseemann/mlst
            abricate: https://github.com/tseemann/abricate
        If executed in a cluster environment the function uses Snakemake and the wrapper class Snek.
        Recommend local runs as they are quick enough for single species, cluster if all are required.
        """

        if self.verbose:
            self._print_parameters(locals())

        analysis_path, fasta_path, mlst_path, res_path, vir_path = self._get_analysis_paths(
            resistance_db, virulence_db, mode="type"
        )

        if force:
            self.logger.info("May the force be with you.")

        self.logger.info("Running analyses on files (.fasta) in", fasta_path)


        if mlst:
            mlst = MLST(target_path=fasta_path, out_dir=mlst_path, exec_path="mlst", force=force)

            mlst.run(min_id=minid, min_cov=mincov)

        if resistance:
            res = Abricate(target_path=fasta_path, out_dir=res_path, exec_path=self.abricate, force=force)
            res.run(db=resistance_db, min_id=minid)

        if virulence:
            vir = Abricate(target_path=fasta_path, out_dir=vir_path, exec_path=self.abricate, force=force)
            vir.run(db=virulence_db, min_id=minid)

    def _get_analysis_paths(self, res_db, vir_db, mode="type"):

        """
        Called by analyse and summarise to define paths available for analysis and paths containing analysis reports.
        If user_path is not defined, return path to species and fasta directories of a project.
        :param user_path: path to analysis directory, requires files at user_path/fasta
        :return: user_path, fasta_path
        """

        analysis_path = None
        fasta_path = None

        if os.path.exists(self.path) and os.path.exists(self.config_path):
            self.logger.info("Operating on species", self.species, "in project", self.path)
            analysis_path = os.path.join(self.species_path, "analysis")
            fasta_path = os.path.join(self.species_path, "fasta")

        elif os.path.exists(self.path) and not os.path.exists(self.config_path):
            self.logger.info("Could not find project at", self.path, "- assuming user path.")

            analysis_path = os.path.join(self.path, "analysis")
            fasta_path = self.path

            if mode == "type":
                fasta_files = [file for file in os.listdir(self.path) if file.endswith(FASTA_EXTENSIONS)]

                if len(fasta_files) == 0:
                    self.logger.info("Could not detect any files (.fasta, .fa) in path", self.path)
                    exit(1)
                else:
                    self.logger.info('Found files (.fasta, .fa) in path', self.path)
            elif mode == "collect":
                if not os.path.exists(analysis_path):
                    self.logger.info("Could not find analysis directory at", analysis_path)
                    self.logger.info("Perhaps you need to run typing first?")
                    exit(1)

        elif not os.path.exists(self.path):
            self.logger.info("Could not find path", self.path)
            exit(1)
        else:
            self.logger.info("Hmmmmm...")
            exit(1)

        # Generate paths to analysis results

        mlst_path = os.path.join(analysis_path, "mlst")
        res_path = os.path.join(analysis_path, "resistance", res_db)
        vir_path = os.path.join(analysis_path, "virulence", vir_db)

        # Return user or species analysis and fasta paths

        return analysis_path, fasta_path, mlst_path, res_path, vir_path

    def _modify_config(self, vf_db, res_db, min_id, min_cov, force):

        with open(self.config_file, "r") as infile:
            config = json.load(infile)

        config["snek"]["force"] = force
        config["abricate_path"] = self.abricate
        config["env"] = self.env_file
        config["abricate_vir_db"] = vf_db
        config["abricate_res_db"] = res_db
        config["minid"] = min_id
        config["mlst_mincov"] = min_cov

        self.config_run_file = os.path.join(self.config_path, "run_" + os.path.basename(self.config_file))

        with open(self.config_run_file, "w") as outfile:
            json.dump(config, outfile)


class NCTCTable:
    """ Class parsing a HTML Soup table into a Pandas dataframe """

    def __init__(self, name, table=None, dataframe=pandas.DataFrame(), links=pandas.DataFrame()):

        """
        Initialize class NCTCTable.
        Calls:
        _transform()
        :param table: HTML soup string containing table section of NCTC3000.
        """

        self.name = name
        self.table = table

        self.dataframe = dataframe
        self.links = links

        if self.table is not None:
            self._transform()

    def store(self, path):

        filename = os.path.join(path, self.name + "_" + time.strftime("%d-%m-%Y-%H-%M"))

        filename_table = filename + "_data.tab"
        filename_links = filename + "_links.tab"

        self.dataframe.to_csv(filename_table, sep="\t")
        self.links.to_csv(filename_links, sep="\t")

    def _transform(self):

        """
        Private method called upon initiation of NCTCTable.
        Assigns the HTML table descriptor ('summary') as name to the class (self.name)
        and extracts data from the HTML table (self.dataframe), including FTP links to
        assembled genomes (self.links).
        Calls:
        _transform_links()
        _clean_tables()
        """

        self.name = self.table["summary"]

        if self.name == "NCTC genomes":
            self.name = "genomes"

        data = [[td.text for td in row.findAll("td")]
                for row in self.table.findAll("tr")]

        ftp_links = [[a["href"] for a in row.find_all("a", href=True)
                      if a["href"].startswith("ftp")]
                     for row in self.table.findAll("tr")]

        links = [self._transform_links(links) for links in ftp_links]

        self.dataframe = pandas.DataFrame(data[1:], columns=[tag.text for tag in self.table.findAll("th")])

        self.links = pandas.DataFrame(links[1:], columns=links[0])

        self._clean_tables()

        assert len(self.dataframe) == len(self.links)

    def summarise(self, top=10):

        """
        Prints a short summary of the current NCTC3000.
        :param top: Number of species and their genome counts to show in summary.
        :return: message: Summary message as formatted string.
        """

        stamp("Here is a summary of the project:")

        message = textwrap.dedent("""
        NCTC3000, {date}
        -------------------------------------------
        Database: {db_name}
        Species:                       {n_species}
        Entries:                       {n_entries}
        Entries with Assembly:         {n_assemblies}
        Entries with Chromosomes:      {n_chromosomes}
        Top {top_number}:
        {top_species_counts}
        -------------------------------------------
        """).format(date=time.strftime("%d-%m-%Y %H:%M"),
                    db_name=self.name,
                    n_species=len(self.dataframe["Species"].unique()),
                    n_entries=len(self.dataframe),
                    n_assemblies=self.dataframe["Manual Assembly"].count(),
                    n_chromosomes=len(self.dataframe[self.dataframe["Chromosomes"]
                                      .astype(float) > 0]),
                    top_number=top,
                    top_species_counts=self.dataframe["Species"].value_counts()[:top].to_string())

        print(message)

        return message

    def _clean_tables(self):

        """
        Need better cleaning of Tables...
        """

        self.dataframe = self.dataframe.replace("Pending", nan)

        self.dataframe.columns = ["Species", "Strain", "Sample", "Runs", "Automated Assembly", "Manual Assembly",
                                  "Chromosomes", "Plasmids", "Unidentified"]

        self.dataframe[["Chromosomes", "Plasmids", "Unidentified"]] = \
            self.dataframe[["Chromosomes", "Plasmids", "Unidentified"]].astype(float)

    @staticmethod
    def _transform_links(links):

        """
        Private static method to transform link data from HTML table into a Pandas dataframe.
        :param links: List of links extracted from HTML table on NCTC3000.
        :return: Dictionary of links from HTML table, missing links are None.
        """

        link_dict = {"manual_gff": nan,
                     "manual_embl": nan,
                     "automatic_gff": nan,
                     "automatic_embl": nan}

        if not links:
            return link_dict
        else:
            for link in links:
                if "manual" and "gff" in link:
                    link_dict["manual_gff"] = link
                if "manual" and "embl" in link:
                    link_dict["manual_embl"] = link
                if "automatic" and "gff" in link:
                    link_dict["automatic_gff"] = link
                if "automatic" and "embl" in link:
                    link_dict["automatic_embl"] = link

        return link_dict


class NCTCProgram:
    """
    Superclass NCTCProgram - minimal wrappers for Abricate and MLST
            htpps://github.com/tseemann/mlst
            htpps://github.com/tseemann/abricate
    """

    def __init__(self, target_path=os.getcwd(), out_dir="mlst", exec_path="mlst", force=False):

        self.out_dir = out_dir
        self.exec_path = exec_path
        self.target_path = target_path

        self.output = None

        if not self.target_path.endswith("/"):
            self.target_path += "/"

        if not self.out_dir.endswith("/"):
            self.out_dir += "/"

        self.go = True

        if os.path.exists(out_dir):
            stamp("Analysis exists at", self.out_dir)
            if force:
                stamp("Force activated, removing directory", self.out_dir)
                shutil.rmtree(self.out_dir)
            else:
                self.go = False

        if self.go:
            stamp("Creating analysis directory", self.out_dir)
            os.makedirs(out_dir)


class MLST(NCTCProgram):

    def run(self, min_cov=80, min_id=90):

        if self.go:
            output = os.path.join(self.out_dir, "mlst.report")

            cmd = " ".join([self.exec_path, "--mincov", str(min_cov), "--minid", str(min_id),
                            self.target_path + "*", ">", output])

            try:
                stamp("Running mlst with database", "on", self.target_path)
                call(cmd, shell=True, stdout=open(os.devnull, 'w'), stderr=open(os.devnull, 'w'))
            except CalledProcessError:
                shutil.rmtree(self.out_dir)
                stamp("Could not run MLST. Removing directory", self.out_dir)
            except KeyboardInterrupt:
                stamp("Interrupted by user, removing analysis directory", self.out_dir)
                shutil.rmtree(self.out_dir, ignore_errors=True)

            return output


class Abricate(NCTCProgram):

    def _get_file_paths(self, file, db):

        file_name, extension = os.path.splitext(file)
        file_path = os.path.join(self.target_path, file)
        output = os.path.join(self.out_dir, file_name + "_" + db + ".tab")

        return file_path, output

    def run(self, db="resfinder", min_id=80):

        if self.go:
            files = [file for file in os.listdir(self.target_path) if file.endswith(FASTA_EXTENSIONS)]

            fail = []
            for file in files:

                file_path, output = self._get_file_paths(file, db)

                cmd = " ".join([self.exec_path, "--db", db, "--minid", str(min_id), "--nopath", file_path, ">", output])

                try:
                    stamp("Running Abricate with database", db, "on", file)
                    call(cmd, shell=True, stdout=open(os.devnull, 'w'), stderr=open(os.devnull, 'w'))
                except CalledProcessError:
                    stamp("Could not run Abricate for file", file)
                    fail.append(file)
                    continue
                except KeyboardInterrupt:
                    stamp("Interrupted by user, removing analysis directory", self.out_dir)
                    shutil.rmtree(self.out_dir, ignore_errors=True)
                finally:
                    if len(fail) == len(files):
                        stamp("Failed to run all files in Abricate, removing directory", self.out_dir)
                        shutil.rmtree(self.out_dir, ignore_errors=True)

            summary_file = self.summarize(db=db)

            return summary_file

    def summarize(self, db):

        if self.go:
            summary = os.path.join(self.out_dir, db + ".report")

            cmd = " ".join([self.exec_path, "--summary", self.out_dir + "*.tab", ">", summary])

            try:
                stamp("Writing results from database", db, "to", summary)
                call(cmd, shell=True, stdout=open(os.devnull, 'w'), stderr=open(os.devnull, 'w'))
            except CalledProcessError:
                raise

            return summary


def stamp(*args):
    print(str(time.strftime("[%H:%M:%S]")) + " " + " ".join([str(arg) for arg in args]))


def get_copy_numbers(frames, cnv_mincov):
    cnv_frames = {}

    for output, dataframe in frames.items():
        if "virulence" in output or "resistance" in output:
            cnv_frame = {}
            for column in dataframe.columns:
                if column not in ("ID", "Count"):
                    cnv_frame[column] = _get_cnv_column(output, dataframe, column, cnv_mincov)
                else:
                    cnv_frame[column] = dataframe[column]

            cnv_frames[output + "_cnv"] = pandas.DataFrame(cnv_frame)

    # Return dictionary without key "mlst", just copy number variation for Abricate
    return cnv_frames


def _get_cnv_column(output, dataframe, column, cnv_mincov):
    cnv_column = []
    for entry in dataframe[column].str.split(";"):
        if len(entry) == 1 and entry[0] == ".":
            cnv_column.append(nan)
        else:
            cov_filtered = []
            for e in entry:
                try:
                    e = float(e)
                    if e >= cnv_mincov:
                        cov_filtered.append(e)
                except ValueError:
                    print("Could not convert coverage in dataframe for", output, "to float.")

            cnv_column.append(len(cov_filtered))

    return cnv_column


def clean_reports(frames):
    clean = {}

    for output, frame in frames.items():
        ids = frame.iloc[:, 0]
        cleaned = frame[ids.str.endswith(FASTA_EXTENSIONS)]
        cleaned_ids = cleaned.iloc[:, 0].map(lambda x: _extract_file_name(x))
        cleaned = _clean_output_reports(output, cleaned)
        cleaned.insert(0, "ID", cleaned_ids)
        cleaned.reindex()

        clean[output] = cleaned

    return clean


def _clean_output_reports(output, dataframe):
    if "virulence" in output or "resistance" in output:
        keep = ~dataframe.columns.str.contains("|".join(['Unnamed:', '#FILE']))
        cleaned = dataframe.iloc[:, keep]
        cleaned.rename(columns={'NUM_FOUND': 'Count'}, inplace=True)

    elif output == "mlst":
        cleaned = dataframe.iloc[:, 1:]
        cleaned.rename(columns={1: 'Species', 2: 'MLST'}, inplace=True)
        cleaned.rename(columns=_get_mlst_header(cleaned), inplace=True)
    else:
        cleaned = None

    return cleaned


def _get_mlst_header(df):
    allele_cols = {}
    i = 1
    for col in df.columns:
        if type(col) is not int64:
            allele_cols[col] = col
        else:
            allele_cols[col] = "MLST_Allele_" + str(i)
            i += 1

    return allele_cols


def _extract_file_name(entry):
    file_name = os.path.basename(entry)
    name, extension = os.path.splitext(file_name)

    return name


def merge_frames(frames, left=None, merge="outer"):
    merged = {}

    for output in frames.keys():
        if "virulence" in output or "resistance" in output:
            if left is None:
                merged_df = pandas.merge(left=frames["mlst"], right=frames[output], left_on="ID", right_on="ID",
                                         how=merge)
            else:
                merged_df = pandas.merge(left=left, right=frames[output], left_on="ID", right_on="ID",
                                         how=merge)

            merged[output] = merged_df

    return merged


def parse_reports(mlst_path, res_path, vir_path, resistance_db="resfinder", virulence_db="vfdb"):
    files = {"mlst": os.path.join(mlst_path, "mlst.report"),
             "resistance_" + resistance_db: os.path.join(res_path, resistance_db + ".report"),
             "virulence_" + virulence_db: os.path.join(vir_path, virulence_db + ".report")}

    frames = {}

    for output, file in files.items():
        if os.path.exists(file):
            stamp("Found output file", file)
            if output == "mlst":
                frames[output] = pandas.read_csv(file, sep="\t", header=None)
            else:
                frames[output] = pandas.read_csv(file, sep="\t")
        else:
            stamp("Could not find output file, skipping", file)

    return frames


def download_assemblies(gff_path, links, force=False):
    exist = 0
    downloaded = 0
    out_paths = []

    for entry in links.itertuples():
        file, url = entry.File, entry.Link

        out_path = os.path.join(gff_path, file)

        if os.path.exists(out_path) and not force:
            exist += 1
        else:
            stamp("Downloading", file, "from", url)
            try:
                urllib.request.urlretrieve(url, out_path)
            except KeyboardInterrupt:
                stamp("Interrupted by user, deleting last file at", out_path)
                if os.path.exists(out_path):
                    os.remove(out_path)
            downloaded += 1

        out_paths.append(out_path)

    stamp("Found", exist, "files and downloaded", downloaded, "new assemblies to", gff_path)

    return out_paths


def gff_to_fasta(gff_files, fasta_path, genbank_path, verbose=True):
    files = {}

    for file in gff_files:

        if verbose:
            stamp("Converting file", file)

        file_name, ext = os.path.splitext(file)
        base_name = os.path.basename(file_name)

        fasta = os.path.join(fasta_path, base_name + ".fasta")
        genbank = os.path.join(genbank_path, base_name + ".gbk")

        if not os.path.exists(fasta) or not os.path.exists(genbank):

            records = _parse_gff(file, file_name)

            if not os.path.exists(fasta):
                SeqIO.write(records, fasta, format="fasta")
            if not not os.path.exists(genbank):
                SeqIO.write(records, genbank, format="genbank")

        files[base_name] = {"gff": os.path.realpath(file),
                            "fasta": os.path.realpath(fasta),
                            "genbank": os.path.realpath(genbank)}

    return files


def _parse_gff(file, file_name):
    records = []

    with open(file, "r") as infile:

        for i, rec in enumerate(GFF.parse(infile)):

            # Enumerates the contigs (can be chromosome, plasmid and unidentified)
            # based on total number of contigs (not type)
            rec_id = rec.id + "_" + str(i + 1)

            if len(rec_id) > 15:
                rec_id = "contig_" + "_" + str(i + 1)

            seq_record = SeqRecord(Seq(str(rec.seq)), id=rec_id,
                                   description=os.path.basename(file_name),
                                   features=rec.features)

            records.append(seq_record)

    return records


def get_links(bacteria, species_indices, assembly):
    links = bacteria.links.loc[species_indices, assembly + "_gff"].dropna()

    file_names = get_file_names(bacteria.dataframe.loc[links.index])

    files = pandas.concat([file_names, links], axis=1)
    files.columns = ["File", "Link"]

    return files


def get_file_names(link_entries):
    species_names = link_entries["Species"].str.replace(" ", "_")
    strain_names = species_names.str.cat([link_entries["Strain"]], sep="_")

    return strain_names.str.cat([".gff" for _ in range(len(link_entries))])