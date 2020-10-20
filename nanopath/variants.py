import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt
import collections
from sklearn.model_selection import train_test_split
from sklearn.metrics import roc_auc_score
from sklearn.ensemble import RandomForestClassifier
from itertools import permutations
import pickle
import vcf
import pyfastx
import logging
from sklearn.metrics import roc_curve
from sklearn.metrics import precision_recall_curve
from sklearn.metrics import average_precision_score
from pysam import VariantFile

from nanopath.utils import PoreLogger, run_cmd
from nanopath.processors import SnippySample, ClairSample, MedakaSample, ForestSample
from pathlib import Path

basesN = ['A', 'C', 'G', 'T']
perm = permutations(basesN, 2)
N, perms = 0, {}
for p in perm:
    perms.setdefault(p[0], {}).setdefault(p[1], N)
    N += 1


def totalReads(row):
    l = row['INFO'].split(';')
    return int(l[2].replace('TotalReads=', ''))


def add_base_change_n(row):
    r = row['ref']
    a = row['alt']
    return perms[r][a]


def alpha_beta(row):
    if row['ref'] == 'A' and row['alt'] == 'G':
        return 0
    elif row['ref'] == 'A' and row['alt'] == 'T':
        return 1
    elif row['ref'] == 'A' and row['alt'] == 'C':
        return 1
    elif row['ref'] == 'C' and row['alt'] == 'A':
        return 1
    elif row['ref'] == 'C' and row['alt'] == 'G':
        return 1
    elif row['ref'] == 'C' and row['alt'] == 'T':
        return 0
    elif row['ref'] == 'G' and row['alt'] == 'A':
        return 0
    elif row['ref'] == 'G' and row['alt'] == 'C':
        return 1
    elif row['ref'] == 'G' and row['alt'] == 'T':
        return 1
    elif row['ref'] == 'T' and row['alt'] == 'A':
        return 1
    elif row['ref'] == 'T' and row['alt'] == 'C':
        return 0
    elif row['ref'] == 'T' and row['alt'] == 'G':
        return 1


def ps(row):
    if row['alt'] == 'N':
        return None
    else:
        return float(row[row['alt'] + '_ps'])

def read_fasta(fasta: Path) -> dict:
    return {
        name: seq.upper() for name, seq in
        pyfastx.Fasta(str(fasta), build_index=False)  # capital bases
    }

class CoreGenome:

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

        self.samples = []  # all samples

    def parse_snp_vcf(self, path: Path, vcf_glob: str = "*.vcf"):

        for vcf in sorted([
            f for f in path.glob(vcf_glob)
        ]):
            # class used after filtering with the classifiers to construct core genome,
            # simply considers snps only regardless of variant caller (default)

            sample = ForestSample(vcf)

            self.logger.info(
                f'Parsed {len(sample.data)} SNPs from sample {sample.name}'
            )
            self.samples.append(sample)

    def core_genome(
        self,
        allow_missing: float = 0.2,
        missing_char: str = "N",
        snp_limit: int = 0,
        include_reference: bool = True
    ):

        """ Determine core variants from Snippy and Medaka samples  """

        # Compute the common core genome SNPs across all VCFs and form into alignment:

        snp_counts = {}
        for sample in self.samples:
            if snp_limit > 0 and len(sample.data) > snp_limit:
                continue

            chrom_snps = {}
            for chrom_name, data in sample.data.groupby('chromosome'):
                chrom_snp_counter = {}
                for _, row in data.iterrows():
                    try:
                        chrom_position = row['position']
                    except (KeyError):
                        raise ValueError(
                            f"Could not parse SNP (index: {_}) on chromosome {chrom_name} in sample {sample.name}"
                        )
                    if chrom_position not in chrom_snp_counter.keys():
                        chrom_snp_counter[chrom_position] = 1
                    else:
                        chrom_snp_counter[chrom_position] += 1

                snp_counter = collections.Counter(chrom_snp_counter)
                chrom_snps[chrom_name] = snp_counter
                self.logger.debug(
                    f"Sample: {sample.name} Chromosome: {chrom_name}, SNPs counted {len(snp_counter)}"
                )

            snp_counts[sample.name] = chrom_snps

        self.logger.info(f"Parsed {len(snp_counts)} samples for core genome alignment")

        if len(snp_counts) == 0:
            raise ValueError(f'Could not parse SNPs with limit: {snp_limit}')


        # Sort by chromosome across samples:
        core_chromosomes = {}
        for sample_name, chrom_snp_counters in snp_counts.items():
            for chrom_name, counts in chrom_snp_counters.items():
                sample_snps = {'sample': sample_name, 'snp_count': counts}
                if chrom_name not in core_chromosomes.keys():
                    core_chromosomes[chrom_name] = [sample_snps]
                else:
                    core_chromosomes[chrom_name].append(sample_snps)

        # Check number of samples and compute SNP counts for each chromosome:
        self.logger.info(f"Parsed {len(core_chromosomes)} chromosomes for core genome alignment")
        core_snp_counts = {}
        core_sample_names = {}
        for chrom_name, sample_snps in core_chromosomes.items():
            self.logger.info(f"Chromosome: {chrom_name} --> {len(sample_snps)} samples")
            if len(sample_snps) != len(snp_counts):
                self.logger.warning(
                    f"There are {len(sample_snps)} / {len(snp_counts)} total samples for chromosome {chrom_name}"
                )
            if chrom_name not in core_snp_counts.keys():
                core_snp_counts[chrom_name] = collections.Counter(
                    sample_snps[0]['snp_count']
                )
                for sample_data in sample_snps[1:]:
                    core_snp_counts[chrom_name].update(
                        collections.Counter(sample_data['snp_count'])
                    )

            core_sample_names[chrom_name] = [sample_data['sample'] for sample_data in sample_snps]

            print(
                f"Unique SNPs across all samples ({len(core_sample_names[chrom_name])}) "
                f"on chromosome {chrom_name}: {len(core_snp_counts[chrom_name])}"
            )

        for chrom_name, core_snp_count in core_snp_counts.items():
            allow_missing_snps = int(len(core_sample_names[chrom_name]) * allow_missing)
            core_threshold_snps = len(core_sample_names[chrom_name]) - allow_missing_snps
            self.logger.info(
                f"Allowing fo {allow_missing*100}% missing SNPs ({allow_missing_snps}) across samples"
            )
            self.logger.info(
                f"Core genome alignment threshold for inclusion: {core_threshold_snps} SNPs)"
            )
            core_snps = collections.Counter({k: c for k, c in core_snp_count.items() if c >= core_threshold_snps})
            print(core_snps)


class RandomForestFilter(PoreLogger):

    """ Train a random-forest classifier to filter variant calls (SNP)

    0. Use ST93 genomes against JKD reference for testing, use all genomes with > 100x coverage later

    1. Subsample training genomes (ONT) that have matching Illumina data to various
        levels of coverage and call SNPs with Medaka and Clair (np-core/np-variants)
    3. Call SNPs on the corresponding Illumina genomes with Snippy, and subset VCF for
        SNPs and broken MNP/COMPLEX (np-core/np-variants)
    4. Compare SNP calls and create a SNP validation file for each genome (based on
        Snippy reference VCF of Illumina calls)
    5.

    """

    def __init__(self, outdir: Path):

        PoreLogger.__init__(self, name="RandomForest", level=logging.INFO)

        self.outdir = outdir
        self.outdir.mkdir(parents=True, exist_ok=True)
        self.model_dir = self.outdir / 'models'
        self.model_dir.mkdir(parents=True, exist_ok=True)
        self.training_dir = self.outdir / 'training'
        self.training_dir.mkdir(parents=True, exist_ok=True)
        self.evaluation_dir = self.outdir / 'evaluation'
        self.evaluation_dir.mkdir(parents=True, exist_ok=True)

        self.ont_calls = dict()  # key: sample name, ont calls with truth column [subclass of Sample]
        self.features = dict()  # key: sample name, feature data frame for ont calls

        self.features_combined = pd.DataFrame()
        self.feature_combinations = {
            'qual': ['quality'],
            'composite': [
                'quality', 'reads_all', 'proximity', 'base_change', 'majority_base_percent',
                'top_base_matches_variant_call', 'deletions_percent', 'insertions_percent'
            ]
        }

        self.probability = 0

        self.vcf_header = ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', 'EXTRA']

    def evaluate_model(
        self,
        model_file: Path,
        vcf_snippy: Path = None,
        vcf_ont: Path = None,
        stats_ont: Path = None,
        dir_snippy: Path = None,
        dir_ont: Path = None,
        caller: str = 'clair',
        prefix: str = "prefix",
        break_complex: bool = True,
        mask_weak: float = 0.8
    ):

        """ Evaluate RF filter model on Snippy reference VCFs and Clair / Medaka variant VCFs """

        self.logger.info(f"Load random forest model: {model_file}")
        with model_file.open('rb') as model_in:
            model = pickle.load(model_in)

        try:
            feature_combo = model_file.name.strip(".sav").split(".")[1]
            self.logger.info(f'Using model feature combination: {feature_combo}')
        except IndexError:
            self.logger.info(f'Could not extract feature combination identifier from: {model_file.name}')
            raise

        try:
            use_features = self.feature_combinations[feature_combo]
        except KeyError:
            self.logger.info(f'Could not find {feature_combo} in native model feature combinations')
            raise

        if dir_snippy and dir_ont:
            comparisons = self.get_evaluation_comparisons(dir_snippy=dir_snippy, dir_ont=dir_ont)
            stats_ont = None
        else:
            comparisons = [vcf_snippy, vcf_ont]
            stats_ont = stats_ont

        self.logger.info(f"Reading files from reference (Snippy) and variant (ONT) callers")
        ont_with_truth, snippies = self.get_data_from_comparisons(
            comparisons=comparisons, caller=caller, break_complex=break_complex, outdir=self.evaluation_dir,
            prefix=prefix, stats=stats_ont  # none if from directory
        )

        self.logger.info(f"Parsing features from variant calls for prediction")
        _, ont_with_features = self.parse_features(ont_calls=ont_with_truth)  # same order as snippy_samples

        classifier_truth_summaries = []
        application_truth_summaries = []
        for i, ont in enumerate(ont_with_features):
            snippy = snippies[i]
            self.logger.info(
                f"Predict SNP validity using the >{feature_combo}< classifier on sample: {ont.name}"
            )

            snp_prediction_features = np.array(
                ont.features[use_features]
            )

            # Prediction with Random Forest Model
            snp_prediction = model.predict(snp_prediction_features)
            snp_prediction_probabilities = model.predict_proba(snp_prediction_features)
            # column header were switched in original code, but no effect of it downstream
            snp_probabilities = pd.DataFrame(snp_prediction_probabilities, columns=[False, True])
            snp_probability = snp_probabilities.max(axis=1)  # extract the corresponding probability

            ont.features['probability'] = snp_probability
            ont.features['prediction'] = snp_prediction  # here the snp prediction from classifier added
            # classify the SNP prediction in comparison to the truth value of the (ONT) SNP call
            self.logger.info(
                f"Evaluate classifier performance to sample {ont.name} vs Snippy reference {snippy.name}"
            )

            if mask_weak > 0:
                print("MASKING WEAK POSITIONS")
                ont_features = ont.features.copy()
                # ont call or N, handled in get_truth_summary
                ont_features['mask_weak'] = [
                    irow['alt'] if irow['ps'] >= mask_weak else 'N'
                    for i, irow in ont_features.iterrows()
                ]
            else:
                ont_features = ont.features.copy()

            # Evaluation processing, weak mask not in snp classifier evaluation!

            ont_features['classifier_evaluation'] = ont_features.apply(self.classify_snp_prediction, axis=1)

            classifier_prediction_evaluations = ont_features.classifier_evaluation.value_counts()

            classifier_truth_summary = self.get_truth_summary(
                true_positives=classifier_prediction_evaluations.get('true_positive'),
                true_negatives=classifier_prediction_evaluations.get('true_negative'),
                false_positives=classifier_prediction_evaluations.get('false_positive'),
                false_negatives=classifier_prediction_evaluations.get('false_negative'),
                snippy=None, ont_data=ont_features, name=ont.name
            )
            classifier_truth_summaries.append(classifier_truth_summary)

            self.logger.info(
                f"Evaluate classifier application to sample {ont.name} vs Snippy reference {snippy.name}"
            )

            # Subset the SNPs by feature prediction as if classifier was applied to sample
            ont.filtered = ont.features[ont.features.prediction == True]

            # Recompute truth against Snippy reference SNPs
            ont, app_summary = self.find_true_snps(
                snippy=snippy, ont=ont, caller=caller, filtered=True
            )  # after RFF, calls get truth summary internally

            application_truth_summaries.append(app_summary)

        classifier_truth_all = pd.DataFrame(classifier_truth_summaries)\
            .set_index('name').sort_values(by=['name', 'coverage'])
        application_truth_all = pd.DataFrame(application_truth_summaries)\
            .set_index('name').sort_values(by=['name', 'coverage'])

        print(classifier_truth_all)
        print(application_truth_all)

        classifier_truth_all.to_csv(self.evaluation_dir / f"{prefix}_classifier_truth.tsv", sep="\t")
        application_truth_all.to_csv(self.evaluation_dir / f"{prefix}_application_truth.tsv", sep="\t")

    def mask_probabilites(self, r):
        if r['probs'] >= self.probability:
            return False
        else:
            return True

    @staticmethod
    def classify_snp_prediction(row):

        """ Classify the prediction of the Illumina SNPs with the RF classifier (i.e. performance of classifier) """

        if row['true_snp'] == True and row['prediction'] == True:
            return 'true_positive'
        elif row['true_snp'] == False and row['prediction'] == True:
            return 'false_positive'
        elif row['true_snp'] == False and row['prediction'] == False:
            return 'true_negative'
        else:
            return 'false_negative'  # can occur here because we check against prediction, not calls this time

    def prepare_training_data(
        self, dir_snippy: Path, dir_ont: Path, caller: str = 'clair',
        break_complex: bool = True, snippy_ext: str = ".ref.vcf"
    ):

        """ Prepare SNP validation data from Snippy reference VCFs and Clair / Medaka variant VCFs """

        comparisons = self.get_coverage_comparisons(dir_snippy=dir_snippy, dir_ont=dir_ont, snippy_ext=snippy_ext)

        ont_with_truth, snippies = self.get_data_from_comparisons(
            comparisons=comparisons, caller=caller, break_complex=break_complex, outdir=self.training_dir
        )

        features, _ = self.parse_features(ont_calls=ont_with_truth)

        # Combined features for training
        self.features_combined = pd.concat(features)  # combined feature frames
        self.features_combined = self.features_combined.reset_index(drop=True)
        self.features_combined.to_csv(self.training_dir / 'training_features.tsv', sep='\t', index=False)

    def train_models(self, test_size: float = 0.3, model_prefix: str = "forest_model"):

        training_features = self.features_combined

        # make sure no missing values (should not be the case)
        training_features = training_features[~training_features['true_snp'].isna()]

        results = {}
        models = {}
        for fc, features in self.feature_combinations.items():
            self.logger.info(f'Training model on feature combination: {fc}')

            # Get features for training and add the validation column
            feature_data_columns = features+['true_snp']
            feature_data = training_features[feature_data_columns]

            # Remove any feature with a missing value
            feature_data = feature_data.dropna()

            # Set up input features (X) and validation (Y)
            x = np.array(feature_data[features])
            y = np.array(feature_data['true_snp'])

            # Generate training and testing data
            self.logger.info(f'Split data into test and train: {fc}')
            train_x, test_x, train_y, test_y = \
                train_test_split(x, y, test_size=test_size, random_state=1)  # TODO random_state

            self.logger.info(f'Training random forest model: {fc}')
            # Train random forest model
            model = RandomForestClassifier()
            model.fit(train_x, train_y)

            self.logger.info(f'Applying model to test data: {fc}')

            # Apply random forest model to testing data
            test_probabilities = model.predict_proba(test_x)

            # Create model evaluation scores
            auc = roc_auc_score(
                y_true=test_y, y_score=test_probabilities[:, 1]
            )
            false_positive_rate, true_positive_rate, thresholds = roc_curve(
                y_true=test_y, y_score=test_probabilities[:, 1]
            )

            model_file = self.model_dir / f"{model_prefix}.{fc}.sav"
            self.logger.info(f'Writing model to file: {model_file}')
            with model_file.open("wb") as model_out:
                pickle.dump(model, model_out)

            self.logger.info(f'Applying model to all data: {fc}')
            # Apply random forest model to all data
            model_predictions = model.predict(x)
            model_probabilities = model.predict_proba(x)

            self.logger.info(f'Assessing feature importances: {fc}')
            # Get and plot feature importances for this model
            importances = model.feature_importances_
            std = np.std([tree.feature_importances_ for tree in model.estimators_], axis=0)
            indices = np.argsort(importances)[::-1]

            self.logger.info(f'Plotting feature importances: {fc}')
            self.plot_feature_importances(importances, std, indices, features, fc)

            results[fc] = {
                'false_positive_rate': false_positive_rate,
                'true_positive_rate': true_positive_rate,
                'AUC': auc,
                'model_probabilities': model_probabilities,  # on all data, train + test
                'model_predictions': model_predictions,  # on all data, train + test
                'test_y': test_y,
                'score_y': test_probabilities[:, 1],
                'model': model
            }
            models[fc] = model

        self.plot_roc_curve(model_results=results)
        self.plot_recall_precision(model_results=results)

        return results, models

    def plot_roc_curve(self, model_results):

        for model_type, results in model_results.items():
            plt.plot(
                results['false_positive_rate'],
                results['true_positive_rate'],
                label=f'{model_type}, AUC: {results["AUC"]:.3f}'
            )
        plt.plot([0, 1], [0, 1], color='darkblue', linestyle='--')
        plt.xlabel('False Positive Rate')
        plt.ylabel('True Positive Rate')
        plt.title('Receiver Operating Characteristic (ROC) Curve')
        plt.legend()
        plt.tight_layout()
        plt.savefig(self.training_dir / 'roc.pdf')
        plt.savefig(self.training_dir / 'roc.png')
        plt.savefig(self.training_dir / 'roc.svg')
        plt.clf()

    def plot_recall_precision(self, model_results):

        for model_type, results in model_results.items():
            average_precision = average_precision_score(results['test_y'], results['score_y'])
            precision, recall, _ = precision_recall_curve(results['test_y'], results['score_y'])
            plt.step(recall, precision, where='post', label=f'{model_type} AP: {average_precision:.3f}')

        plt.xlabel('Recall')
        plt.ylabel('Precision')
        plt.ylim([0.0, 1.05])
        plt.xlim([0.0, 1.0])
        plt.title('Precision-Recall curve')
        plt.legend()
        plt.tight_layout()
        plt.savefig(self.training_dir / 'pr.pdf')
        plt.savefig(self.training_dir / 'pr.png')
        plt.savefig(self.training_dir / 'pr.svg')
        plt.clf()

    def plot_feature_importances(self, importances, std, indices, features, fc):

        plt.figure()
        plt.title("Relative importance of features")

        df = pd.DataFrame(
            {'Importance': importances, 'std': std, 'indices': indices, 'feature': features}
        )
        df.to_csv(self.training_dir / f'{fc}_feat_importance.tsv')

        sns.barplot(y='feature', x="Importance", xerr=df['std'], capsize=.2, data=df)

        plt.tight_layout()

        plt.savefig(self.training_dir / f'{fc}_feature_importance.pdf')
        plt.savefig(self.training_dir / f'{fc}_feature_importance.png')
        plt.savefig(self.training_dir / f'{fc}_feature_importance.svg')
        plt.clf()

    @staticmethod
    def get_vcf_files(dir_snippy, dir_ont, snippy_ext: str = ".ref.vcf", ont_ext: str = ".vcf"):

        if not dir_snippy.exists():
            raise ValueError('Could not find Snippy VCF directory')
        if not dir_ont.exists():
            raise ValueError('Could not find ONT VCF directory')

        snippy_files = list(dir_snippy.glob(f'*{snippy_ext}'))
        ont_files = list(dir_ont.glob(f'*{ont_ext}'))  # output from coverage subsamples in training is format {name}_{coverage}.vcf

        if not snippy_files:
            raise ValueError('Could not find Snippy VCF')
        if not ont_files:
            raise ValueError('Could not find ONT VCF')

        return snippy_files, ont_files

    def get_evaluation_comparisons(self, dir_snippy, dir_ont, snippy_ext: str = ".ref.vcf"):

        snippy_files, ont_files = self.get_vcf_files(dir_snippy=dir_snippy, dir_ont=dir_ont, snippy_ext=snippy_ext)

        print(snippy_files, ont_files)

        snippy_file_names = [f.name.replace(snippy_ext, "") for f in snippy_files]
        ont_file_names = [f.stem for f in ont_files]  # .vcf

        print(snippy_file_names, ont_file_names)

        comparisons = [
            (dir_snippy / f"{name}{snippy_ext}", dir_ont / f"{name}.vcf")
            for name in snippy_file_names if name in ont_file_names
        ]

        if not comparisons:
            raise ValueError('Could not find matching files to compare')

        self.logger.info(
            f'Detected {len(comparisons)} matching sample VCFs (Snippy - ONT)'
        )

        return comparisons

    def get_coverage_comparisons(self, dir_snippy, dir_ont, snippy_ext):

        snippy_files, ont_files = self.get_vcf_files(dir_snippy=dir_snippy, dir_ont=dir_ont, snippy_ext=snippy_ext)

        # Get ONT call sample names and coverage
        ont_coverage = {}
        for file in ont_files:
            name, coverage = file.stem.split("_")
            if name not in ont_coverage.keys():
                ont_coverage[name] = [coverage]
            else:
                ont_coverage[name].append(coverage)

        # Check if all coverages are matching
        coverage_lengths = {
            name: len(coverages) for name, coverages in ont_coverage.items()
        }
        if len(set(coverage_lengths.values())) != 1:
            raise ValueError(
                f'Could not detect matching coverage for each file: {coverage_lengths}'
            )

        snippy_ref_names = set(
            [f.stem.replace(".ref", "").split("_")[0] for f in snippy_files]
        )

        for name in snippy_ref_names:
            if name not in ont_coverage.keys():
                raise ValueError(f'Could not detect matching ONT calls for Snippy reference sample: {name}')

        self.logger.info(
            f'Detected {len(snippy_ref_names)} matching sample VCFs (Snippy - ONT)'
        )

        comparisons = []
        for name in snippy_ref_names:
            ont_covs = ont_coverage[name]
            for cov in ont_covs:
                comparisons.append(
                    (dir_snippy / f"{name}_{cov}.ref.vcf", dir_ont / f"{name}_{cov}.vcf")
                )

        if not comparisons:
            raise ValueError('Could not find matching files to compare')

        self.logger.info(
            f"{len(comparisons)} subsampled VCFs (ONT) matched {len(snippy_ref_names)} reference VCFs (Snippy)"
        )

        return comparisons

    def get_data_from_comparisons(
        self, comparisons, caller: str, break_complex: bool, outdir: Path, prefix: str = "base", stats: Path = None
    ):

        truth_summaries = []
        ont_with_truth = []
        snippies = []
        for comparison in comparisons:
            snippy, ont = self.read_samples_from_comparison(
                comparison=comparison, caller=caller, break_complex=break_complex, stats=stats
            )
            ont, truth_summary = self.find_true_snps(snippy=snippy, ont=ont, caller=caller)

            ont_with_truth.append(ont)
            snippies.append(snippy)
            truth_summaries.append(truth_summary)

        truth_comparison = pd.DataFrame(truth_summaries).set_index('name').sort_values(by=['name', 'coverage'])

        print(truth_comparison)

        truth_comparison.to_csv(outdir / f'{prefix}_{caller}_truth.tsv', sep='\t', index=True)

        return ont_with_truth, snippies

    @staticmethod
    def read_samples_from_comparison(comparison: tuple, caller: str = "clair", break_complex: bool = True, stats: Path = None):

        snippy_vcf, ont_vcf = comparison

        snippy = SnippySample(vcf=snippy_vcf, break_complex=break_complex)

        if stats is None:
            pysamstats = ont_vcf.parent / f"{ont_vcf.stem}.txt"
        else:
            pysamstats = stats

        if not pysamstats.exists():
            raise ValueError(
                f'Could not find associated stats file ({pysamstats}) for calls: {ont_vcf.name}'
            )

        if caller == "medaka":
            ont = MedakaSample(vcf=ont_vcf, stats=pysamstats)
        elif caller == "clair":
            ont = ClairSample(vcf=ont_vcf, stats=pysamstats)
        else:
            raise ValueError('Caller argument must be one of: clair or medaka')

        return snippy, ont

    def find_true_snps(self, snippy, ont, caller, filtered: bool = False):

        if filtered:
            data = ont.filtered
        else:
            data = ont.data

        if "mask_weak" in data.columns.tolist():
            weak = len(data[data['mask_weak'] == "N"])
            self.logger.info(f"Found {weak} weak positions in {ont.name}")
            data = data[data['mask_weak'] != "N"]

        true_snps = snippy.data.merge(
            data, on=['chromosome', 'position'], how='inner',
            suffixes=["_snippy", f"_{caller}_filtered" if filtered else f"_{caller}"]
        )

        true_indices = [
            (row['chromosome'], row['position']) for _, row in true_snps.iterrows()
        ]

        truth_column = []
        for i, row in data.iterrows():
            idx = (row['chromosome'], row['position'])
            if idx in true_indices:
                truth_column.append(True)
            else:
                truth_column.append(False)

        data['true_snp_filtered' if filtered else 'true_snp'] = truth_column
        data['name_filtered' if filtered else 'name'] = [ont.name for _ in truth_column]

        true_positives = len(true_snps)
        false_positives = len(data) - len(true_snps)
        false_negatives = len(snippy.data) - len(true_snps)
        true_negatives = 0  # always 0, no way to validate false Snippy calls with nanopore calls

        if filtered:
            ont.filtered = data
        else:
            ont.data = data

        truth_summary = self.get_truth_summary(
            true_positives=true_positives,
            false_positives=false_positives,
            false_negatives=false_negatives,
            true_negatives=true_negatives,
            snippy=snippy, ont_data=data, name=ont.name,  # rff filtered or unfiltered data from ONT calls
        )

        return ont, truth_summary

    def get_truth_summary(
        self, true_positives, false_positives, false_negatives, true_negatives, snippy, ont_data, name: str
    ):

        """ Get comparison statistics on call or prediction performance: accuracy, precision, recall, f1"""

        # If extracted from dictionary with dict.get()
        if true_positives is None:
            true_positives = 0
        if false_positives is None:
            false_positives = 0
        if false_negatives is None:
            false_negatives = 0
        if true_negatives is None:
            true_negatives = 0

        accuracy = (true_positives + true_negatives) / \
                   (true_positives + true_negatives + false_positives + false_negatives)

        precision = true_positives / (true_positives + false_positives)
        recall = true_positives / (true_positives + false_negatives)

        f1 = 2 * (recall * precision) / (recall + precision)

        try:
            _, cov = name.strip(".vcf").split("_")
        except ValueError:
            cov = None

        if cov is not None:
            try:
                cov = int(cov)
            except TypeError:
                self.logger.debug("Could not convert parsed coverage value to integer")
                raise

        results = {
            'name': name,
            'coverage': cov,
            'ont': len(ont_data),
            'true_positives': true_positives,
            'false_positives': false_positives,
            'false_negatives': false_negatives,
            'true_negatives': true_negatives,
            'accuracy': accuracy,
            'precision': precision,
            'recall': recall,
            'f1': f1
        }

        if snippy is not None:
            results['snippy'] = len(snippy.data)

        return results

    @staticmethod
    def parse_features(ont_calls: list = None):

        # False positives, no chained assignments, therefore disable warning
        pd.options.mode.chained_assignment = None

        features = []
        sample_features = []
        for sample in ont_calls:

            data = sample.data

            # Must compute proximity by chromosome! 222073

            chrdf = []
            for chromosome, df in data.groupby('chromosome'):
                df['alt_length'] = df.alt.map(len)
                df['ref_length'] = df.ref.map(len)
                df['indel_length'] = df['alt_length'] - df['ref_length']

                df['5prime'] = df['position'] - df['position'].shift(1)
                df['3prime'] = df['position'] - df['position'].shift(-1)
                df['3prime'] = df['3prime'].abs()

                df['proximity'] = df[['5prime', '3prime']].min(axis=1)
                df['proximity'].fillna(5000, inplace=True)  # not sure what this does

                # Ensures that only SNPs are processed
                # should be the case from parsing anyway
                df = df[df.ref_length == 1]
                df = df[df.alt_length == 1]

                df['alpha_beta'] = df.apply(alpha_beta, axis=1)  # ref <-> alt transition or transversion
                df['base_change'] = df.apply(add_base_change_n, axis=1)  # hmmmmm

                # Pysamstats features added from sample calls object
                df = df.merge(
                    sample.stats, left_on=['chromosome', 'position'],
                    right_on=['chrom', 'pos'], how='left', suffixes=[None, "_stats"]
                )

                df['total'] = df['A'] + df['C'] + df['G'] + df['T']
                for b in ('A', 'T', 'C', 'G', 'insertions', 'deletions'):
                    df[f'{b}_percent'] = (df[b] / df['reads_all'])*100
                    df[f'{b}_ps'] = df[b] / df['total']

                df['reads_all'] = df[['A', 'T', 'C', 'G']].sum(axis=1)
                df['top_base'] = df[['A', 'T', 'C', 'G']].max(axis=1)
                df['top_base_seq'] = df[['A', 'T', 'C', 'G']].idxmax(axis=1)

                df['majority_base_percent'] = (df['top_base'] / df['reads_all'])
                df['top_base_matches_variant_call'] = np.where(df.alt == df.top_base_seq, 1, 0)

                df = df[df['majority_base_percent'].notna()]

                if len(df) > 0:
                    df['ps'] = df.apply(ps, axis=1)
                else:
                    raise ValueError(f"No remaining features in {sample.name} on chromosome {chromosome}")

                chrdf.append(df)

            chromosome_features = pd.concat(chrdf)

            sample.features = chromosome_features  # store features in ont calls object
            sample_features.append(sample)

            features.append(chromosome_features)

        return features, sample_features


class ForestClassifier:
    def __init__(
        self, vcf: Path = None, stats: Path = None, output: Path = None, model_file: Path = None,
        caller: str = None, mask_weak: bool = False, probability: float = 0
    ):

        self.vcf = vcf
        self.stats = stats
        self.output = output

        self.model_file = model_file

        if model_file:
            self.model = self.load_model()

        self.caller = caller

        self.probability = probability
        self.mask_weak = mask_weak

        if self.vcf:
            self.chromosome = self.get_chromosome()  # TODO: single chromosome what if multiple?

        self.keep = set()

        self.vcf_header = ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', 'EXTRA']

    def get_chromosome(self):

        with self.vcf.open("r") as fi:
            return [
                ln.split(",")[0].replace('##contig=<ID=', '')
                for ln in fi if ln.startswith('##contig=<ID=')
            ][0]

    def load_model(self):
        return pickle.load(
            self.model_file.open('rb')
        )

    def mask_probabilites(self, r):
        if r['probs'] >= self.probability:
            return False
        else:
            return True

    def ps(self, row):
        if row['ALT'] == 'N':
            x = None
        else:
            x = float(row[row['ALT'] + '_PS'])
        return x


    def read_vcf(self):

        df = pd.read_csv(self.vcf, sep='\t', comment='#', names=self.vcf_header)

        if len(df) == 0:
            raise ValueError("Could not read data from VCF")

        df['ALT_len'] = df.ALT.map(len)
        df['REF_len'] = df.REF.map(len)
        df['INDEL length'] = df['ALT_len'] - df['REF_len']

        df['5prime proximity'] = df['POS'] - df['POS'].shift(1)
        df['3prime proximity'] = df['POS'] - df['POS'].shift(-1)

        df['3prime proximity'] = df['3prime proximity'].abs()
        df['proximty'] = df[['5prime proximity', '3prime proximity']].min(axis=1)
        df['proximty'].fillna(5000, inplace=True)
        df['alphabeta'] = df.apply(alpha_beta, axis=1)

        df = df[df.REF_len == 1]
        df = df[df.ALT_len == 1]
        if len(df) > 0:
            df['baseChange'] = df.apply(add_base_change_n, axis=1)

        psdf = pd.read_csv(self.stats, sep='\t')
        df = df.merge(psdf, left_on=['CHROM', 'POS'], right_on=['chrom', 'pos'], how='left')

        df['total'] = df['A'] + df['C'] + df['G'] + df['T']
        bases = ['A', 'T', 'C', 'G', 'insertions', 'deletions']
        for b in bases:
            df['{0} %'.format(b)] = (df[b] / df['reads_all']) * 100
            df['{0}_PS'.format(b)] = (df[b] / df['total'])
        df['reads_all'] = df[['A', 'T', 'C', 'G']].sum(axis=1)
        df['top_base'] = df[['A', 'T', 'C', 'G']].max(axis=1)
        df['top_base_seq'] = df[['A', 'T', 'C', 'G']].idxmax(axis=1)
        df['majority base %'] = (df['top_base'] / df['reads_all'])
        df['Top Base matches Nanopolish'] = np.where(df.ALT == df.top_base_seq, 1, 0)
        df = df[df['majority base %'].notna()]
        if len(df) > 0:
            df['ps'] = df.apply(self.ps, axis=1)


        self.SNPs = df
        return len(df)

    def classify(self):
        feature_combinations = {
            'medaka': ['QUAL', 'reads_all', 'proximty', 'baseChange', 'majority base %',
                       'Top Base matches Nanopolish', 'deletions %', 'insertions %'],
            'clair': ['QUAL', 'reads_all', 'proximty', 'baseChange', 'majority base %',
                      'Top Base matches Nanopolish', 'deletions %', 'insertions %']
        }
        features = feature_combinations[self.caller]
        print(self.SNPs[features])
        X = np.array(self.SNPs[features])
        preds = self.model.predict(X)
        probs = self.model.predict_proba(X)
        probs = pd.DataFrame(probs, columns=[True, False])
        p = probs.max(axis=1)
        self.SNPs['preds'] = preds
        self.SNPs['probs'] = p
        #        self.SNPs['preds']=self.SNPs.apply(self.maskProbs,axis=1)
        self.SNPs['mask'] = self.SNPs.apply(self.mask_probabilites, axis=1)
        keep = self.SNPs[self.SNPs.preds == True]
        s = set(keep.POS)
        self.keep.update(s)
        mask = self.SNPs[self.SNPs['mask'] == True]
        psmask = self.SNPs[self.SNPs['ps'] < 0.8]
        self.mask = set(mask['POS'])
        print(len(self.mask))
        self.mask.update(list(psmask['POS']))
        print(len(self.mask))

    def _vcf_reader(self):
        vcf_reader = vcf.Reader(self.vcf.open('r'))
        for record in vcf_reader:
            yield record

    def filter(self):
        vcf_reader = self._vcf_reader()
        vcf_oneread = vcf.Reader(self.vcf.open('r'))
        vcf_writer = vcf.Writer(self.output.open('w'), vcf_oneread)
        for record in vcf_reader:
            if record.POS not in self.keep:
                continue
            else:
                if self.mask_weak:
                    if record.POS in self.mask:
                        record.ALT = 'N'
                        print('Masking weak position:', record.POS)
                vcf_writer.write_record(record)

        vcf_writer.close()
