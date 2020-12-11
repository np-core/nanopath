""" Beastling v0.1

Wrapper for constrained models of bacterial transmission dynamics in BEAST
from SNP alignments produced with Snippy and Gubbins in PhyBeast Nextflow

Authors: Eike Steinig, Michael Meehan, Sebastian Duchene

"""

import jinja2
import pandas
import pyfastx
import logging
import textwrap
import yaml
import math
import uuid

from pathlib import Path
from nanopath.utils import PoreLogger


class BeastlingError(Exception):
    pass


class Beastling(PoreLogger):

    """ Beastling configures XML files for constrained bacterial models """

    def __init__(
        self,
        alignment: Path, data: Path,
        clock_model: str = 'strict',
        sample_every: int = 1000,
        chain_length: int = 100000,
        chain_type: str = 'default',
        chain_number: int = 4,
        prefix: str = "bdss",
        keep_missing: bool = True,
        sample_prior: bool = False
    ):

        """ Initiate Beastling

        :param alignment: core genome SNP alignment FASTA

        :param data: tab-delimited sample data with
            minimum headers: name, date

        :param clock_model: clock model: one of the following
            strict, relaxed_lognormal, relaxed_exponential

        :param keep_missing: keep sequences that have no date meta data and
            set their dates to missing (0) - otherwise remove from alignment

        """

        PoreLogger.__init__(self, level=logging.INFO, name='Beastling')

        self.version = "2.6"
        self.logger.info(f'Beastling for BEAST v{self.version}')

        # Input data and alignment
        self.data = self.read_data(data=data)
        self.alignment = self.read_fasta(fasta=alignment)

        # Substitution model
        self.substitution_model = 'general_time_reversible'
        self.gamma_categories = 4

        # Clock model
        self.clock_model = clock_model
        self.clock = None  # configured with priors later

        if sample_prior:
            self.sample_from_prior = "true"
        else:
            self.sample_from_prior = "false"

        # Demographic model
        self.model_priors = dict()  # configured in model subclasses

        # MCMC settings
        self.chain_type = chain_type
        self.chain_length = chain_length
        self.chain_number = chain_number
        self.sample_every = sample_every

        # Loggers
        self.posterior_log = f"{prefix}.log"
        self.posterior_every = sample_every  # same as sample frequency
        self.tree_log = f"{prefix}.trees"
        self.tree_every = sample_every  # same as sample frequency

        self.prefix = prefix

        # Housekeeping
        self.keep_missing = keep_missing  # keep samples with missing dates

        # Information
        self.alignment_file = alignment
        self.data_file = data

        # Constraints
        self.allowed_distributions = (
            'lognormal', 'exponential', 'beta', 'gamma', 'uniform'
        )
        self.allowed_gamma_modes = {
            'shape_mean': 'ShapeMean',
            'shape_scale': 'ShapeScale',
            'shape_rate': 'ShapeRate',
            'one_parameter': 'OneParameter'
        }
        self.allowed_clock_priors = (
            'rate', 'uced', 'ucld_mean', 'ucld_sd'
        )
        self.allowed_interval_params = (
            'sampling_proportion',
            'reproductive_number',
            'become_uninfectious',
            'rho'
        )
        self.required_fields = {
            'all': ('initial', 'lower', 'upper', 'dimension'),
            'interval': ('initial',),
            'lognormal': ('mean', 'sd', 'real_space'),
            'exponential': ('mean',),
            'gamma': ('alpha', 'beta', 'mode'),
            'beta': ('alpha', 'beta'),
            'uniform': ()
        }

        self.allowed_model_priors = ()  # defined in subclasses

    # Data input methods

    @staticmethod
    def read_fasta(fasta: Path) -> dict:
        return {
            name: seq.upper() for name, seq in
            pyfastx.Fasta(str(fasta), build_index=False)  # capital bases
        }

    @staticmethod
    def read_data(data: Path):

        """ Parse sample meta data file, missing dates are set to 0

        :param data: path to tab-delimited sample meta data: name, date, type
        :return: dataframe of sample meta data
        """

        df = pandas.read_csv(
            data, sep='\t', header=0, na_values=['-']
        )

        # TODO: ad-hoc fix - need to assign individual priors to missing dates
        df['date'] = df['date'].fillna(
            max(df.date)
        )

        return df

    @staticmethod
    def read_config(file: Path) -> dict:

        """ Read prior settings from YAML """

        with file.open() as fh:
            return yaml.load(fh, Loader=yaml.FullLoader)

    @staticmethod
    def load_template(file: str):

        """ Load the template XML file for BDSS """

        template_loader = jinja2.FileSystemLoader(
            searchpath=f"{Path(__file__).parent / 'templates'}"
        )
        template_env = jinja2.Environment(loader=template_loader)
        template = template_env.get_template(file)

        return template

    def print_configuration(self):

        """ Print initial configuration to terminal """

        self.logger.info(f'Input alignment: {self.alignment_file}')
        self.logger.info(f'Input sample data: {self.data_file}')
        self.logger.info(
            f'Substitution model: {self.substitution_model} '
            f'({self.gamma_categories})'
        )
        self.logger.info(f'Molecular clock model: {self.clock_model}')

        if self.chain_type == 'coupled':
            self.logger.info(f'MCMCMC chain length: {self.chain_length}')
            self.logger.info(f'Active hot chains: {self.chain_number-1}')
        else:
            self.logger.info(f'MCMC chain length: {self.chain_length}')

        self.logger.info(f'Posterior log: {self.posterior_log}')
        self.logger.info(f'Log posterior every: {self.posterior_every}')
        self.logger.info(f'Tree log: {self.tree_log}')
        self.logger.info(f'Log trees every: {self.tree_every}')

    def check_configuration(self):

        """ Some input data and configuration checks """

        self.logger.info('Checking input and configurations')

        # Check clock model setting
        if self.clock_model not in (
            'strict', 'relaxed_exponential', 'relaxed_lognormal'
        ):
            raise BeastlingError(
                'Clock model must be one of: strict, '
                'relaxed_exponential, relaxed_lognormal'
            )

        # Check meta data headers: name and date
        if not any(x in ['date'] for x in self.data.columns):
            raise BeastlingError(
                'Sample data file must contain one column '
                'with header: date'
            )
        if not any(x in ['name'] for x in self.data.columns):
            raise BeastlingError(
                'Sample data file must contain one column '
                'with header: name'
            )

        # Check consistent sequence lengths in alignment
        sequence_lengths = [len(seq) for _, seq in self.alignment.items()]
        if len(set(sequence_lengths)) != 1:
            raise BeastlingError(
                'Sequences in input alignment are not the same length'
            )

        alignment_names = [name for name, _ in self.alignment.items()]

        # Subset the meta data table by sequences present in alignment
        # so that we don't have to worry when using it later
        self.data = self.data[self.data['name'].isin(alignment_names)]

        self.logger.info(
            f'Sequences in alignment: {len(alignment_names)}'
        )
        self.logger.info(
            f'Entries in sample data: {len(self.data)}'
        )

    # Configuration of chain, date and data blocks in XML

    def get_chain_xml(self) -> str:

        """ Construct an XML for the selected type of Markov Chain """

        self.logger.info('XML: setting chain parameters')

        if self.chain_type == 'coupled':
            return f'<run id="mcmc" spec="beast.coupledMCMC.CoupledMCMC" ' \
                f'chainLength="{self.chain_length}" ' \
                f'chains="{self.chain_number}" ' \
                f'target="0.234" logHeatedChains="false" deltaTemperature="0.1" ' \
                f'optimise="true" resampleEvery="{self.sample_every}" >'
        else:
            return f'<run id="mcmc" spec="MCMC" chainLength="{self.chain_length}" sampleFromPrior="{self.sample_from_prior}">'

    def get_date_xml(self) -> str:

        """ Construct the date string for the state XML """

        # TODO: does this need to be integers in the string?

        self.logger.info('XML: generating date string for state block')

        date_strings = [
            f'{row["name"]}={row["date"]}'
            for i, row in self.data.iterrows()
        ]

        return ",".join(date_strings)

    def get_data_xml(self) -> str:

        """ Construct a XML string of the input sequence data block """

        self.logger.info('XML: generating sequence data block')

        data_block = ""
        for name, seq in self.alignment.items():

            # Check how many unique nucleotides in sequence
            # and fail if they are not exclusively ACTG
            bases = set(seq)
            for base in bases:
                if base not in ('A', 'C', 'T', 'G', 'N'):
                    raise BeastlingError(
                        f'Sequence for {name} contains base '
                        f'other than ACTGN: {bases}'
                    )

            # TODO: check if total count parameter is right or always four
            data_block += f'<sequence ' \
                f'id="seq_{name}" ' \
                f'spec="Sequence" ' \
                f'taxon="{name}" ' \
                f'totalcount="{len(bases)}" ' \
                f'value="{seq}"/>\n'

        return data_block

    # Data integrity and config checks

    def check_prior_config(self, prior_config: dict, distribution: bool = True):

        for prior, prior_data in prior_config.items():

            # if prior not in self.allowed_model_priors:
            #     raise BeastlingError(
            #         f'Prior `{prior}` is not allowed for this model'
            #     )

            # Some models (Coalescent Bayesian Skyline) do not need
            # a distribution field specified
            if distribution:
                try:
                    distr = prior_data['distribution']
                except KeyError:
                    raise BeastlingError(
                        'Prior configuration requires a `distribution` field'
                    )

                # Check allowed distributions
                if distr not in self.allowed_distributions:
                    raise BeastlingError(
                        f'Distribution `{distr}` not supported.'
                    )

                # Check specific required fields:
                for field in self.required_fields[distr]:
                    if field not in prior_data.keys():
                        raise BeastlingError(
                            f'Field `{field}` is missing in `{prior}`'
                            f'prior for distribution `{distr}`'
                        )

            # Check general required fields for model and clock priors only
            if not prior.startswith('interval'):
                required_fields = 'all'
            else:
                required_fields = 'interval'

            for field in self.required_fields[required_fields]:
                if field not in prior_data.keys():
                    raise BeastlingError(
                        f'Field `{field}` is missing in `{prior}` prior'
                    )

            # Check allowed values for Gamma distribution field `mode`
            if 'mode' in prior_data.keys():
                if prior_data.get('mode') not in self.allowed_gamma_modes.keys():
                    raise BeastlingError(
                        f'Allowed values for `mode` are: '
                        f'{self.allowed_gamma_modes.keys()}'
                    )

        self.logger.info('Prior configuration: OK')

    def check_clock_priors(self, prior_config: dict):

        for p in prior_config.keys():
            if p not in self.allowed_clock_priors:
                raise BeastlingError(
                    f'Prior `{p}` is not a supported prior for clock model'
                )

        self.check_prior_config(prior_config=prior_config, distribution=True)
        self.logger.info('Selected clock priors: OK')

    def check_model_priors(self, prior_config: dict, distribution: bool = True):

        for p in prior_config.keys():
            if p not in self.model_priors.keys():
                raise BeastlingError(
                    f'Prior `{p}` is not a supported prior for demographic model'
                )

        self.check_prior_config(
            prior_config=prior_config, distribution=distribution
        )
        self.logger.info('Selected model priors: OK')

    # Configuration methods for model and clock priors

    def set_model_priors(self, prior_config: dict, distribution: bool = True):

        self.check_model_priors(
            prior_config=prior_config, distribution=distribution
        )

        if distribution:
            distributions = self.load_distributions(prior_config=prior_config)
        else:
            distributions = dict()  # e.g. Coalescent Bayesian Skyline

        loaded = []
        for p, prior_data in prior_config.items():
            self.model_priors[p].distribution = distributions.get(p)
            self.model_priors[p].initial = prior_data.get('initial')
            self.model_priors[p].lower = prior_data.get('lower')
            self.model_priors[p].upper = prior_data.get('upper')
            self.model_priors[p].dimension = prior_data.get('dimension')

            self.model_priors[p].check_init()

            loaded.append(p)

        required = set(self.model_priors.keys())
        loaded = set(loaded)

        if loaded != required:
            for p in loaded.symmetric_difference(required):
                self.logger.info(
                    f'Could not detect prior `{p}` in configuration'
                )
            raise BeastlingError(
                f'Could not detect required priors in configuration (logged)'
            )

        self.logger.info('Configured model priors: OK')

    def set_clock(self, prior_config: dict):

        self.check_clock_priors(prior_config=prior_config)
        distributions = self.load_distributions(prior_config=prior_config)

        if self.clock_model == 'strict':
            clock = Strict(
                priors=[Rate()]
            )
        elif self.clock_model == 'relaxed_exponential':
            clock = RelaxedExponential(
                priors=[UCED()]
            )
        elif self.clock_model == 'relaxed_lognormal':
            clock = RelaxedLogNormal(
                priors=[UCLDMean(), UCLDSD()]
            )
        else:
            raise BeastlingError(
                f'Clock model `{self.clock_model}` not supported'
            )

        for prior in clock.priors:
            name = prior.name
            try:
                config = prior_config[name]
            except KeyError:
                raise BeastlingError(
                    f'Could not detect clock prior setting ({name}) '
                    f'for clock model ({self.clock_model})'
                )
            try:
                distribution = distributions[name]
            except KeyError:
                raise BeastlingError(
                    f'Could not detect clock distribution ({name}) '
                    f'for clock model ({self.clock_model})'
                )

            prior.distribution = distribution
            prior.initial = config.get('initial')
            prior.lower = config.get('lower')
            prior.upper = config.get('upper')
            prior.dimension = config.get('dimension')

            if 'fixed' in config.keys():
                # Set the whole clock to a fixed value (no operators)
                clock.fixed = config.get('fixed', False)
                self.logger.info(f'Set a fixed clock at rate: {prior.initial}')

            prior.check_init()

        self.clock = clock
        self.logger.info('Configured clock priors: OK')

    # Distribution helper method

    def load_distributions(self, prior_config: dict) -> dict:

        distributions = {}
        for prior, prior_data in prior_config.items():
            distribution = prior_data['distribution']

            # Distribution selection for parameters of prior
            if distribution == 'lognormal':
                distributions[prior] = LogNormal(
                    mean=prior_data.get('mean'),
                    sd=prior_data.get('sd'),
                    real_space=prior_data.get('real_space')
                )
            elif distribution == 'exponential':
                distributions[prior] = Exponential(
                    mean=prior_data.get('mean')
                )
            elif distribution == 'gamma':
                distributions[prior] = Gamma(
                    alpha=prior_data.get('alpha'),
                    beta=prior_data.get('beta'),
                    mode=self.allowed_gamma_modes[
                        prior_data.get('mode')
                    ]
                )
            elif distribution == 'beta':
                distributions[prior] = Beta(
                    alpha=prior_data.get('alpha'),
                    beta=prior_data.get('beta'),
                )
            elif distribution == 'uniform':
                distributions[prior] = Uniform()
            else:
                raise BeastlingError(
                    f'Distribution {distribution} is not supported'
                )

        self.logger.info('Distribution configuration: OK')

        return distributions


class CoalescentBayesianSkyline(Beastling):

    def __init__(
        self,
        alignment: Path,
        data: Path,
        clock_model: str,
        chain_type: str,
        chain_length: int,
        chain_number: int,
        prefix: str
    ):
        Beastling.__init__(
            self,
            alignment=alignment,
            data=data,
            clock_model=clock_model,
            chain_type=chain_type,
            chain_length=chain_length,
            chain_number=chain_number,
            prefix=prefix
        )

        self.allowed_model_priors = ('pop_size', 'group_size')

        self.logger.info('Coalescent Bayesian Skyline')

        # Model priors
        self.model_priors = dict(
            pop_size=PopulationSize(),
            group_size=GroupSize()
        )

    def construct_template(self, xml: Path):

        template = self.load_template(file='cosky.xml')

        self.check_model_dimensions()

        ps_prior = self.model_priors.get('pop_size')
        gs_prior = self.model_priors.get('group_size')

        rendered = template.render(
            data_xml=self.get_data_xml(),
            date_xml=self.get_date_xml(),
            mcmc_xml=self.get_chain_xml(),
            tree_log=self.tree_log,
            tree_every=self.tree_every,
            posterior_log=self.posterior_log,
            posterior_every=self.posterior_every,
            pop_size_param=ps_prior.get_parameter_xml(),
            group_size_state_node=gs_prior.get_group_size_state_node(),
            clock_param=self.clock.get_parameter_xml(),
            clock_prior=self.clock.get_prior_xml(),
            clock_state_node=self.clock.get_state_node_xml(),
            clock_branch_rate=self.clock.get_branch_rate_xml(),
            clock_scale_operator=self.clock.get_scale_operator_xml(),
            clock_updown_operator=self.clock.get_updown_operator_xml(),
            clock_logger=self.clock.get_logger_xml(),
        )

        with xml.open('w') as fh:
            fh.write(rendered)

    def check_model_dimensions(self):

        """ Model dimensions for group size and pop size need to be the same """

        if self.model_priors['pop_size'].dimension != self.model_priors['group_size'].dimension:
            raise BeastlingError(
                'Dimensions of `pop_size` and `group_size` need to be the same'
            )


class BirthDeathSkylineSerial(Beastling):

    """ Birth Death Skyline Serial (BDSS)

    Specifically, the model can be used to estimate:

        - date of origin
        - reproductive number
        - sampling proportion
        - become uninfectious rate

    """

    def __init__(
        self,
        alignment: Path,
        data: Path,
        clock_model: str,
        chain_type: str,
        chain_length: int,
        chain_number: int,
        prefix: str,
        sample_prior: bool
    ):

        Beastling.__init__(
            self,
            alignment=alignment,
            data=data,
            clock_model=clock_model,
            chain_type=chain_type,
            chain_length=chain_length,
            chain_number=chain_number,
            prefix=prefix,
            sample_prior=sample_prior
        )

        self.allowed_model_priors = [
            'sampling_proportion',
            'reproductive_number',
            'become_uninfectious',
            'origin'
        ]

        self.logger.info('Birth-Death Skyline Serial')

        # Model priors
        self.model_priors = dict(
            origin=Origin(),
            reproductive_number=ReproductiveNumberBDSS(),
            sampling_proportion=SamplingProportionBDSS(),
            become_uninfectious=BecomeUninfectiousBDSS()
        )

    def construct_template(self, xml: Path, r_dimension: int = None):

        template = self.load_template(file='bdss.xml')

        og_prior = self.model_priors.get('origin')
        rn_prior = self.model_priors.get('reproductive_number')
        sp_prior = self.model_priors.get('sampling_proportion')
        bu_prior = self.model_priors.get('become_uninfectious')

        slice_function_xml, slice_rate_xml, slice_logger_xml \
            = self.get_slice_xml()

        rendered = template.render(
            data_xml=self.get_data_xml(),
            date_xml=self.get_date_xml(),
            mcmc_xml=self.get_chain_xml(),
            tree_log=self.tree_log if r_dimension is None else f"{self.prefix}_d{r_dimension}.trees",
            tree_every=self.tree_every,
            posterior_log=self.posterior_log if r_dimension is None else f"{self.prefix}_d{r_dimension}.log",
            posterior_every=self.posterior_every,
            origin_param=og_prior.get_parameter_xml(),
            origin_prior=og_prior.get_prior_xml(),
            reproductive_number_param=rn_prior.get_parameter_xml(d=r_dimension),
            reproductive_number_prior=rn_prior.get_prior_xml(),
            sampling_proportion_param=sp_prior.get_parameter_xml(),
            sampling_proportion_prior=sp_prior.get_prior_xml(),
            become_uninfectious_param=bu_prior.get_parameter_xml(),
            become_uninfectious_prior=bu_prior.get_prior_xml(),
            clock_param=self.clock.get_parameter_xml(),
            clock_prior=self.clock.get_prior_xml(),
            clock_state_node=self.clock.get_state_node_xml(),
            clock_branch_rate=self.clock.get_branch_rate_xml(),
            clock_scale_operator=self.clock.get_scale_operator_xml(),
            clock_updown_operator=self.clock.get_updown_operator_xml(),
            clock_logger=self.clock.get_logger_xml(),
            slice_functions=slice_function_xml,
            slice_rates=slice_rate_xml,
            slice_loggers=slice_logger_xml
        )

        with xml.open('w') as fh:
            fh.write(rendered)

    def get_slice_xml(self) -> (str, str, str):
        """ Helper function to get XML for slices from all configured priors """

        slice_priors = [
            self.model_priors.get('reproductive_number'),
            self.model_priors.get('become_uninfectious'),
            self.model_priors.get('sampling_proportion'),
        ]

        reverse_time_array = \
            '<reverseTimeArrays spec="beast.core.parameter.BooleanParameter" ' \
            'value="{0} {1} {2} false false"/>'.format(
                str(slice_priors[0].sliced).lower(),
                str(slice_priors[1].sliced).lower(),
                str(slice_priors[2].sliced).lower()
            )

        slice_function_xml, slice_rate_xml, slice_logger_xml = "", "", ""
        for p in slice_priors:
            slice_function_xml += p.get_slice_function_xml()
            slice_rate_xml += p.get_slice_rate_xml() + '\n'
            slice_logger_xml += p.get_slice_logger_xml()

        if len(slice_function_xml) == 0:
            reverse_time_array = ""

        slice_rate_xml += reverse_time_array

        return slice_function_xml, slice_rate_xml, slice_logger_xml

    # Slice configuration methods

    def check_slice_config(self, slice_config: dict):

        for interval_parameter in slice_config.keys():
            if interval_parameter not in self.allowed_interval_params:
                raise BeastlingError(
                    'Parameter for intervals must be one of: '
                    'sampling_proportion, reproductive_number or '
                    'become_uninfectious'
                )

        for interval_parameter, parameter_data in slice_config.items():

            if 'change_times' not in parameter_data.keys():
                raise BeastlingError(
                    f'Setting `change_times` is missing in '
                    f'interval parameter {interval_parameter}'
                )
            else:
                if not isinstance(parameter_data.get('change_times'), list):
                    raise BeastlingError(
                        f'Interval change times for {interval_parameter} '
                        f'must be specified as list'
                    )
                else:
                    number_changes = len(
                        parameter_data.get('change_times')
                    )

            if 'intervals' not in parameter_data.keys():
                raise BeastlingError(
                    f'Setting `intervals` is missing in interval '
                    f'parameter {interval_parameter}'
                )
            else:
                if not isinstance(parameter_data.get('intervals'), list):
                    raise BeastlingError(
                        'Interval change times must be specified as list'
                    )
                else:
                    number_intervals = len(
                        parameter_data.get('intervals')
                    )

            if number_changes != number_intervals:
                raise BeastlingError(
                    'Number of change times is not equal to number of intervals'
                )

            slices = parameter_data.get('intervals')
            self.check_prior_config(prior_config={
                f'interval_{i}': pd for i, pd in enumerate(slices)
            }, distribution=True)  # Make this a format for configuration checks

        self.logger.info('Interval configuration: OK')

    def set_slices(self, slice_config: dict):

        """ Should always be called AFTER setting priors """

        self.check_slice_config(slice_config=slice_config)

        for parameter, data in slice_config.items():

            if parameter == 'sampling_proportion':
                prior = SamplingProportionBDSS()
                prior.lower = 0.
                prior.upper = 1.0
            elif parameter == 'become_uninfectious':
                prior = BecomeUninfectiousBDSS()
                prior.lower = 0.
                prior.upper = math.inf
            elif parameter == 'reproductive_number':
                prior = ReproductiveNumberBDSS()
                prior.lower = 0.
                prior.upper = math.inf
            else:
                raise BeastlingError(
                    'Parameter for interval slicing not set correctly.'
                )

            distributions = self.load_distributions(prior_config={
                f'interval_{i}': pd for i, pd in enumerate(
                    data.get('intervals')
                )
            })  # Make this a format for distribution loading, returns dict

            prior.sliced = True
            prior.intervals = data.get('change_times')
            prior.dimension = len(prior.intervals)
            prior.initial = [pd.get('initial') for pd in data.get('intervals')]

            # Because of distribution dict, make sure keys are sorted
            sorted_distribution_keys = sorted(
                distributions.keys(), key=lambda x: x.split('_')[-1]
            )

            try:
                prior.distribution = [
                    distributions[p] for p in sorted_distribution_keys
                ]
            except IndexError:
                raise BeastlingError(
                    "Could not find entry in adjusted prior keys - "
                    "something is very wrong"
                )

            prior.check_init()
            self.model_priors[parameter] = prior


class BirthDeathSkylineContemporary(Beastling):

    """ Birth Death Skyline Contemporary (BDSC)

    Specifically, the model can be used to estimate:

        - date of origin
        - reproductive number
        - sampling probability (rho)
        - become uninfectious rate

    """

    def __init__(
        self,
        alignment: Path,
        data: Path,
        clock_model: str,
        chain_type: str,
        chain_length: int,
        chain_number: int,
        prefix: str,
        sample_prior: bool
    ):

        Beastling.__init__(
            self,
            alignment=alignment,
            data=data,
            clock_model=clock_model,
            chain_type=chain_type,
            chain_length=chain_length,
            chain_number=chain_number,
            prefix=prefix,
            sample_prior=sample_prior
        )

        self.allowed_model_priors = [
            'rho',
            'reproductive_number',
            'become_uninfectious',
            'origin'
        ]

        self.logger.info('Birth-Death Skyline Contemporary')

        # Model priors
        self.model_priors = dict(
            origin=OriginBDSC(),
            reproductive_number=ReproductiveNumberBDSC(),
            rho=RhoBDSC(),
            become_uninfectious=BecomeUninfectiousBDSC()
        )

    def construct_template(self, xml: Path):

        template = self.load_template(file='bdsc.xml')

        og_prior = self.model_priors.get('origin')
        rn_prior = self.model_priors.get('reproductive_number')
        sp_prior = self.model_priors.get('rho')
        bu_prior = self.model_priors.get('become_uninfectious')

        slice_function_xml, slice_rate_xml, slice_logger_xml \
            = self.get_slice_xml()

        rendered = template.render(
            data_xml=self.get_data_xml(),
            date_xml=self.get_date_xml(),
            mcmc_xml=self.get_chain_xml(),
            tree_log=self.tree_log,
            tree_every=self.tree_every,
            posterior_log=self.posterior_log,
            posterior_every=self.posterior_every,
            origin_param=og_prior.get_parameter_xml(),
            origin_prior=og_prior.get_prior_xml(),
            reproductive_number_param=rn_prior.get_parameter_xml(),
            reproductive_number_prior=rn_prior.get_prior_xml(),
            rho_param=sp_prior.get_parameter_xml(),
            rho_prior=sp_prior.get_prior_xml(),
            become_uninfectious_param=bu_prior.get_parameter_xml(),
            become_uninfectious_prior=bu_prior.get_prior_xml(),
            clock_param=self.clock.get_parameter_xml(),
            clock_prior=self.clock.get_prior_xml(),
            clock_state_node=self.clock.get_state_node_xml(),
            clock_branch_rate=self.clock.get_branch_rate_xml(),
            clock_scale_operator=self.clock.get_scale_operator_xml(),
            clock_updown_operator=self.clock.get_updown_operator_xml(),
            clock_logger=self.clock.get_logger_xml(),
            slice_functions=slice_function_xml,
            slice_rates=slice_rate_xml,
            slice_loggers=slice_logger_xml
        )

        with xml.open('w') as fh:
            fh.write(rendered)

    def get_slice_xml(self) -> (str, str, str):
        """ Helper function to get XML for slices from all configured priors """

        slice_priors = [
            self.model_priors.get('reproductive_number'),
            self.model_priors.get('become_uninfectious'),
            self.model_priors.get('rho'),
        ]

        reverse_time_array = \
            '<reverseTimeArrays spec="beast.core.parameter.BooleanParameter" ' \
            'value="{0} {1} {2} false false"/>'.format(
                str(slice_priors[0].sliced).lower(),
                str(slice_priors[1].sliced).lower(),
                str(slice_priors[2].sliced).lower()
            )

        slice_function_xml, slice_rate_xml, slice_logger_xml = "", "", ""
        for p in slice_priors:
            slice_function_xml += p.get_slice_function_xml()
            slice_rate_xml += p.get_slice_rate_xml() + '\n'
            slice_logger_xml += p.get_slice_logger_xml()

        if len(slice_function_xml) == 0:
            reverse_time_array = ""

        slice_rate_xml += reverse_time_array

        return slice_function_xml, slice_rate_xml, slice_logger_xml

    # Slice configuration methods

    def check_slice_config(self, slice_config: dict):

        for interval_parameter in slice_config.keys():
            if interval_parameter not in self.allowed_interval_params:
                raise BeastlingError(
                    'Parameter for intervals must be one of: '
                    'sampling_proportion, reproductive_number or '
                    'become_uninfectious'
                )

        for interval_parameter, parameter_data in slice_config.items():

            if 'change_times' not in parameter_data.keys():
                raise BeastlingError(
                    f'Setting `change_times` is missing in '
                    f'interval parameter {interval_parameter}'
                )
            else:
                if not isinstance(parameter_data.get('change_times'), list):
                    raise BeastlingError(
                        f'Interval change times for {interval_parameter} '
                        f'must be specified as list'
                    )
                else:
                    number_changes = len(
                        parameter_data.get('change_times')
                    )

            if 'intervals' not in parameter_data.keys():
                raise BeastlingError(
                    f'Setting `intervals` is missing in interval '
                    f'parameter {interval_parameter}'
                )
            else:
                if not isinstance(parameter_data.get('intervals'), list):
                    raise BeastlingError(
                        'Interval change times must be specified as list'
                    )
                else:
                    number_intervals = len(
                        parameter_data.get('intervals')
                    )

            if number_changes != number_intervals:
                raise BeastlingError(
                    'Number of change times is not equal to number of intervals'
                )

            slices = parameter_data.get('intervals')
            self.check_prior_config(prior_config={
                f'interval_{i}': pd for i, pd in enumerate(slices)
            }, distribution=True)  # Make this a format for configuration checks

        self.logger.info('Interval configuration: OK')

    def set_slices(self, slice_config: dict):

        """ Should always be called AFTER setting priors """

        self.check_slice_config(slice_config=slice_config)

        for parameter, data in slice_config.items():

            if parameter == 'rho':
                prior = RhoBDSC()
                prior.lower = 0.
                prior.upper = 1.0
            elif parameter == 'become_uninfectious':
                prior = BecomeUninfectiousBDSC()
                prior.lower = 0.
                prior.upper = math.inf
            elif parameter == 'reproductive_number':
                prior = BecomeUninfectiousBDSC()
                prior.lower = 0.
                prior.upper = math.inf
            else:
                raise BeastlingError(
                    'Parameter for interval slicing not set correctly.'
                )

            distributions = self.load_distributions(prior_config={
                f'interval_{i}': pd for i, pd in enumerate(
                    data.get('intervals')
                )
            })  # Make this a format for distribution loading, returns dict

            prior.sliced = True
            prior.intervals = data.get('change_times')
            prior.dimension = len(prior.intervals)
            prior.initial = [pd.get('initial') for pd in data.get('intervals')]

            # Because of distribution dict, make sure keys are sorted
            sorted_distribution_keys = sorted(
                distributions.keys(), key=lambda x: x.split('_')[-1]
            )

            try:
                prior.distribution = [
                    distributions[p] for p in sorted_distribution_keys
                ]
            except IndexError:
                raise BeastlingError(
                    "Could not find entry in adjusted prior keys - "
                    "something is very wrong"
                )

            prior.check_init()
            self.model_priors[parameter] = prior


class MultiTypeBirthDeath(Beastling):

    def __init__(
        self,
        alignment: Path,
        data: Path,
        clock_model: str,
        chain_type: str,
        chain_length: int,
        chain_number: int,
        prefix: str,
        trait: str,
    ):

        Beastling.__init__(
            self,
            alignment=alignment,
            data=data,
            clock_model=clock_model,
            chain_type=chain_type,
            chain_length=chain_length,
            chain_number=chain_number,
            prefix=prefix
        )

        self.allowed_model_priors = [
            'sampling_proportion',
            'reproductive_number',
            'become_uninfectious',
            'rate_matrix'
        ]

        self.logger.info('Multi-Type Birth-Death')

        self.trait = trait
        self.demes: int = 0
        self.change_times: list = []
        self.migration_rates: list = []

        self.map_tree_log = Path(self.tree_log).stem + "_map.trees"
        self.typed_node_tree_log = Path(self.tree_log).stem + "_typed_node.trees"

        # Model priors
        self.model_priors = dict(
            rate_matrix=RateMatrix(),
            reproductive_number=ReproductiveNumberMTBD(),
            sampling_proportion=SamplingProportionMTBD(),
            become_uninfectious=BecomeUninfectiousMTBD()
        )

    def get_type_xml(self) -> str:

        """ Construct the type string for the state XML in MTBD """

        self.logger.info('XML: generating type string for state block')

        date_strings = [
            f'{row["name"]}={row[self.trait]}'
            for i, row in self.data.iterrows()
        ]

        return ",".join(date_strings)

    # Configuration methods

    def set_model_config(self, model_config: dict):

        demes, change_times, migration_rates = \
            self.check_model_config(model_config=model_config)

        self.demes = demes
        self.change_times = change_times
        self.migration_rates = migration_rates

    @staticmethod
    def check_model_config(model_config: dict):

        try:
            demes = model_config['demes']
        except KeyError:
            raise BeastlingError(
                "MTBD model requires setting of model config: `demes`"
            )

        try:
            change_times = model_config['change_times']
        except KeyError:
            raise BeastlingError(
                "MTBD model requires setting of model config: `change_times`"
            )

        try:
            migration_rates = model_config['migration_rates']
        except KeyError:
            raise BeastlingError(
                "MTBD model requires setting of model config: `migration_rates`"
            )

        if not isinstance(change_times, list):
            raise BeastlingError('Model parameter `change_times` must be a list')
        if not isinstance(migration_rates, list):
            raise BeastlingError('Model parameter `migration_rates` must be a list')

        expected_rates = demes*(demes-1)
        if len(migration_rates) != expected_rates:
            raise BeastlingError(
                f'Number of migration rates ({migration_rates}) must be '
                f'equivalent to: number of demes * number of demes -1 '
                f'({expected_rates})'
            )

        if not len(change_times) > 1:
            raise BeastlingError(
                'Number of change times must be > 1'
            )

        return demes, change_times, migration_rates

    def check_deme_config(self):

        """ Are demes and associated parameters are properly configured? """

        nsp = len(self.model_priors['sampling_proportion'].initial)
        nct = len(self.change_times)

        if nsp != self.demes*nct:
            raise BeastlingError(
                f'Number of sampling proportion initial values must be '
                f'of length: {self.demes*nct}'
            )

        # Make sure dimensions are set correctly
        for mp, prior in self.model_priors.items():
            if mp in (
                'reproductive_number',
                'sampling_proportion',
                'become_uninfectious'
            ):

                if not prior.dimension == len(prior.initial):
                    self.logger.warning(
                        'Dimension must be the same length '
                        'as the number of initial values!'
                    )
                    self.logger.warning(
                        f'Changing dimension of prior `{prior}` from '
                        f'{prior.dimension} to {len(prior.initial)}'
                    )
                    prior.dimension = len(prior.initial)

        types = self.data[self.trait].unique().tolist()

        if len(types) != self.demes:
            raise BeastlingError(
                f'Number of types specified in data file ({len(types)}) does'
                f'not match the number of configured demes ({self.demes})'
            )

    def get_sampling_slice_xml(self):

        """ Default sampling slices in multi-type birth-death model """

        try:
            change_times = [float(t) for t in self.change_times]
        except TypeError:
            raise BeastlingError(
                'Change times must be floats'
            )
        dim = len(change_times)
        change_times.reverse()
        change_times = " ".join([str(t) for t in change_times])

        return f"""
        <parameter id="samplingRateChangeTimes.s:snp" spec="parameter.RealParameter" dimension="{dim}" name="samplingRateChangeTimes">{change_times}</parameter>
        <reverseTimeArrays id="BooleanParameter.0" spec="parameter.BooleanParameter" dimension="6">false false true false false false</reverseTimeArrays>
        """

    def get_geofreq_xml(self):

        # TODO: what's going on with dimensions and name in case demes > 2 - see BEAUTi XMLs

        values = " ".join(
            [format(1/v, ".5f") for v in self.demes*[self.demes]]
        )

        return f"""
        <parameter id="geo-frequencies.t:snp" spec="parameter.RealParameter" dimension="2" lower="0.0" name="stateNode" upper="1.0">{values}</parameter>
        """

    def construct_template(self, xml: Path):

        # Note the default template gamma shape prior was set with limits to 0
        # (upper + lower) in the standard output of the model in BEAUTi -
        # changed to Infinity but make sure this is intended

        # Similarly the rateMatrix prior is supposed to be bounded 0 - 100
        # but in the default output is set to Infinity

        template = self.load_template(file='mtbd.xml')

        rm_prior = self.model_priors.get('rate_matrix')
        rn_prior = self.model_priors.get('reproductive_number')
        sp_prior = self.model_priors.get('sampling_proportion')
        bu_prior = self.model_priors.get('become_uninfectious')

        rendered = template.render(
            data_xml=self.get_data_xml(),
            date_xml=self.get_date_xml(),
            type_xml=self.get_type_xml(),
            mcmc_xml=self.get_chain_xml(),
            tree_log=self.tree_log,
            map_tree_log=self.map_tree_log,
            typed_node_tree_log=self.typed_node_tree_log,
            tree_every=self.tree_every,
            posterior_log=self.posterior_log,
            posterior_every=self.posterior_every,
            geo_frequency_param=self.get_geofreq_xml(),
            reproductive_number_param=rn_prior.get_parameter_xml(),
            reproductive_number_prior=rn_prior.get_prior_xml(),
            sampling_proportion_param=sp_prior.get_parameter_xml(),
            sampling_proportion_prior=sp_prior.get_prior_xml(),
            become_uninfectious_param=bu_prior.get_parameter_xml(),
            become_uninfectious_prior=bu_prior.get_prior_xml(),
            rate_matrix_param=rm_prior.get_parameter_xml(),
            rate_matrix_prior=rm_prior.get_prior_xml(),
            clock_param=self.clock.get_parameter_xml(),
            clock_prior=self.clock.get_prior_xml(),
            clock_state_node=self.clock.get_state_node_xml(),
            clock_branch_rate=self.clock.get_branch_rate_xml(),
            clock_scale_operator=self.clock.get_scale_operator_xml(),
            clock_updown_operator=self.clock.get_updown_operator_xml(),
            clock_logger=self.clock.get_logger_xml(),
            sampling_rate_slice=self.get_sampling_slice_xml()
        )

        with xml.open('w') as fh:
            fh.write(rendered)

# Distribution classes


class Distribution:
    """ Distribution base class type checking and later extensions """
    pass


class Exponential(Distribution):
    def __init__(
        self, mean: float = 10.0
    ):
        self.mean = mean

    def get_xml(self):

        return f"""
            <Exponential id="Exponential.{uuid.uuid4()}" name="distr">
                <parameter id="RealParameter.{uuid.uuid4()}" spec="parameter.RealParameter" estimate="false" name="mean">{self.mean}</parameter>
            </Exponential>
        """


class Uniform(Distribution):
    def __init__(self):
        pass

    @staticmethod
    def get_xml():
        return f"""
            <Uniform id="Uniform.{uuid.uuid4()}" name="distr"/>
        """


class LogNormal(Distribution):
    def __init__(
        self, mean: float = 1.0, sd: float = 1.25, real_space: bool = False
    ):
        self.real_space = str(real_space).lower()
        self.mean = mean
        self.sd = sd

    def get_xml(self):

        return f"""
            <LogNormal id="LogNormalDistributionModel.{uuid.uuid4()}" meanInRealSpace="{self.real_space}" name="distr">
                <parameter id="RealParameter.{uuid.uuid4()}" spec="parameter.RealParameter" estimate="false" name="M">{self.mean}</parameter>
                <parameter id="RealParameter.{uuid.uuid4()}" spec="parameter.RealParameter" estimate="false" name="S">{self.sd}</parameter>
            </LogNormal>
        """


class Gamma(Distribution):
    def __init__(
        self, alpha: float = 2.0, beta: float = 2.0, mode: str = 'ShapeMean',
    ):
        self.alpha = alpha
        self.beta = beta
        self.mode = mode

    def get_xml(self):

        return f"""
            <Gamma id="Gamma.{uuid.uuid4()}" mode="{self.mode}" name="distr">
                <parameter id="RealParameter.{uuid.uuid4()}" spec="parameter.RealParameter" estimate="false" name="alpha">{self.alpha}</parameter>
                <parameter id="RealParameter.{uuid.uuid4()}" spec="parameter.RealParameter" estimate="false" name="beta">{self.beta}</parameter>
            </Gamma>
        """


class Beta(Distribution):
    def __init__(
        self, alpha: float = 1.0, beta: float = 1.0
    ):
        self.alpha = alpha
        self.beta = beta

    def get_xml(self):

        return f"""
            <Beta id="Beta.{uuid.uuid4()}" name="distr">
                <parameter id="RealParameter.{uuid.uuid4()}" spec="parameter.RealParameter" estimate="false" name="alpha">{self.alpha}</parameter>
                <parameter id="RealParameter.{uuid.uuid4()}" spec="parameter.RealParameter" estimate="false" name="beta">{self.beta}</parameter>
            </Beta>
        """


# Prior base class


class Prior:

    """ Base class for priors """

    def __init__(
        self,
        distribution: Distribution or [Distribution] = None,
        initial: float or list = None,
        lower: float or str = None,
        upper: float or str = None,
        dimension: int = 1
    ):

        self.distribution = distribution  # load a distribution subclass

        self.initial = initial
        self.lower = lower
        self.upper = upper

        self.dimension = dimension

        self.sliced = False  # for sliced SamplingProportion prior
        self.intervals = []  # for sliced SamplingProportion prior

        self.name = None  # prior identifier name defined in model subclasses
        self.idx = None  # prior identifier defined in all prior subclasses
        self.x = None  # param identifier defined in all prior subclasses

        self.scw = None  # scale operator weight defined in clock subclasses
        self.scx = None  # scale identifier defined in clock subclasses

        self.allowed_negative_infinity = ('-infinity', '-inf', '-Infinity')
        self.allowed_positive_infinity = ('infinity', 'inf', 'Infinity')

        self.param_spec = "parameter.RealParameter"

    def get_bounds(self):

        if self.lower == -math.inf:
            lower = 'lower="-Infinity"'
        elif self.lower is None:
            lower = ''
        elif self.lower == math.inf:
            raise BeastlingError(
                'Lower bound cannot be positive infinity'
            )
        else:
            lower = f'lower="{self.lower}"'

        if self.upper == math.inf:
            upper = 'upper="Infinity"'
        elif self.upper is None:
            upper = ''
        elif self.lower == -math.inf:
            raise BeastlingError(
                'Upper bound cannot be negative infinity'
            )
        else:
            upper = f'upper="{self.upper}"'

        return lower, upper

    def check_init(self):

        # Initial value checks
        if self.initial in self.allowed_negative_infinity \
                or self.initial in self.allowed_positive_infinity:
            raise BeastlingError(
                f'Initial value cannot be: {self.initial}'
            )

        if isinstance(self.initial, list):
            try:
                self.initial = [float(i) for i in self.initial]
            except TypeError:
                raise BeastlingError(
                    f'Initial values in list must be valid floats'
                )
        else:
            try:
                if self.param_spec == "parameter.RealParameter":
                    self.initial = float(self.initial)
                elif self.param_spec == "parameter.IntegerParameter":
                    self.initial = int(self.initial)
                else:
                    raise BeastlingError(
                        f'Parameter specification: `{self.param_spec}` not found'
                    )
            except TypeError:
                raise BeastlingError(
                    f'Initial value must be a valid float'
                )

        try:
            self.intervals = [float(i) for i in self.intervals]
        except TypeError:
            raise BeastlingError(
                f'Interval values in list must be valid floats'
            )

        if self.lower is not None or self.upper is not None:
            # Upper / lower infinity checks
            if self.lower in self.allowed_positive_infinity:
                raise BeastlingError(
                    f'Lower bound cannot be positive infinity'
                )

            if self.upper in self.allowed_negative_infinity:
                raise BeastlingError(
                    f'Upper bound cannot be negative infinity'
                )

            # Upper / lower value checks and conversion
            if self.lower in self.allowed_negative_infinity:
                self.lower = -math.inf
            else:
                try:
                    self.lower = float(self.lower)
                except TypeError:
                    raise BeastlingError(
                        'Lower bound value must be valid float or infinity string'
                    )

            if self.upper in self.allowed_positive_infinity:
                self.upper = math.inf
            else:
                try:
                    self.upper = float(self.upper)
                except TypeError:
                    raise BeastlingError(
                        'Upper bound value must be valid float or infinity string'
                    )

            # Initital bounds check
            if isinstance(self.initial, list):
                check_initial = self.initial
            else:
                check_initial = [self.initial]

            for i in check_initial:
                if not self.lower <= i <= self.upper:
                    raise BeastlingError(
                        f'Initial value {i} must be within '
                        f'bounds: {self.lower} and {self.upper}'
                    )

        if self.sliced:
            if not isinstance(self.initial, list):
                raise BeastlingError(
                    'Something went wrong, slices set but '
                    'initial values is not a list'
                )

            if self.dimension != len(self.initial) or \
                    self.dimension != len(self.intervals):
                raise BeastlingError(
                    'In sliced SamplingProportion, the number of dimensions '
                    'must match the number of initial values and intervals'
                )

            if len(self.initial) != len(self.intervals):
                raise BeastlingError(
                    'Length of initial values must be the same length of '
                    'intervals when setting intervals of sampling proportion'
                )

    def get_parameter_xml(self, d: int = None) -> str:

        lower, upper = self.get_bounds()

        # Allow for higher dimensions using slices
        if isinstance(self.initial, list):
            initial = " ".join(str(i) for i in self.initial)
        else:
            initial = self.initial

        return textwrap.dedent(
            f'<parameter id="{self.x}:snp" spec="{self.param_spec}" '
            f'dimension="{d if d else self.dimension}" {lower} {upper} '
            f'name="stateNode">{initial}</parameter>'
        )

    def get_prior_xml(self) -> str:

        if not self.sliced:
            # Normal singular prior distribution
            return f"""
                <prior id="{self.idx}:snp" name="distribution" x="@{self.x}:snp">
                    {self.distribution.get_xml()}
                </prior>
            """
        else:
            # Sliced sampling proportion distribution per interval
            sliced_priors = ''
            for i, distribution in enumerate(self.distribution):
                sliced_priors += f"""
                    <prior id="{self.name}Prior{i+1}" name="distribution" x="@{self.name}{i+1}">
                        {distribution.get_xml()}
                    </prior>
                """
            return sliced_priors

    def get_logger_xml(self) -> str:

        return f'<log idref="{self.x}:snp"/>'

    # Clock scale operator for priors
    def get_scale_operator_xml(self):
        pass

    # Sliced priors
    def get_slice_function_xml(self) -> str:

        if not self.sliced:
            return ''
        else:
            xml = ''
            for i, interval in enumerate(self.initial):
                xml += f'<function spec="beast.core.util.Slice" id="{self.name}{i+1}" ' \
                    f'arg="@{self.x}:snp" index="{i}" count="1"/>\n'
            return xml

    def get_slice_rate_xml(self) -> str:

        if not self.sliced:
            return ''
        else:
            intervals = " ".join(str(i) for i in self.intervals)
            if self.name == 'samplingProportion':
                rate_change_times = 'samplingRateChangeTimes'
            elif self.name == 'rho':
                rate_change_times = 'samplingRateChangeTimes'
            elif self.name == 'reproductiveNumber':
                rate_change_times = 'birthRateChangeTimes'
            elif self.name == 'becomeUninfectious':
                rate_change_times = 'deathRateChangeTimes'
            else:
                raise BeastlingError(
                    'Rate change times are only defined for: '
                    'rho, samplingProportion, reproductiveNumber and becomeUninfectious'
                )

            return f'<{rate_change_times} spec="beast.core.parameter.RealParameter" value="{intervals}"/>'

    def get_slice_logger_xml(self) -> str:

        if not self.sliced:
            return ''
        else:
            loggers = ''
            for i, value in enumerate(self.initial):
                loggers += f'<log idref="{self.name}{i+1}"/>\n'
            return loggers

# Clock model priors


class Rate(Prior):
    def __init__(self):
        Prior.__init__(self)
        self.name = 'rate'
        self.idx = "ClockPrior.c"
        self.x = "clockRate.c"


class UCED(Prior):
    def __init__(self):
        Prior.__init__(self)
        self.name = 'uced'
        self.idx = "UCMeanRatePrior.c"
        self.x = "ucedMean.c"


class UCLDMean(Prior):
    def __init__(self):
        Prior.__init__(self)
        self.name = 'ucld_mean'
        self.idx = "MeanRatePrior.c"
        self.x = "ucldMean.c"


class UCLDSD(Prior):
    def __init__(self):
        Prior.__init__(self)
        self.name = 'ucld_sd'
        self.idx = "ucldStdevPrior.c"
        self.x = "ucldStdev.c"


# Model priors for Birth-Death Skyline Serial

class Origin(Prior):
    def __init__(self):
        Prior.__init__(self)
        self.name = "origin"
        self.idx = f"originPrior_BDSKY_Serial.t"
        self.x = f"origin_BDSKY_Serial.t"


class ReproductiveNumberBDSS(Prior):
    def __init__(self):
        Prior.__init__(self)
        self.name = "reproductiveNumber"
        self.idx = "reproductiveNumberPrior_BDSKY_Serial.t"
        self.x = "reproductiveNumber_BDSKY_Serial.t"


class SamplingProportionBDSS(Prior):
    def __init__(self):
        Prior.__init__(self)
        self.name = 'samplingProportion'
        self.idx = "samplingProportionPrior_BDSKY_Serial.t"
        self.x = "samplingProportion_BDSKY_Serial.t"


class BecomeUninfectiousBDSS(Prior):
    def __init__(self):
        Prior.__init__(self)
        self.name = 'becomeUninfectious'
        self.idx = "becomeUninfectiousRatePrior_BDSKY_Serial.t"
        self.x = "becomeUninfectiousRate_BDSKY_Serial.t"

# Contemp priors


class OriginBDSC(Prior):
    def __init__(self):
        Prior.__init__(self)
        self.name = "origin"
        self.idx = f"originPrior_BDSKY_Contempt"
        self.x = f"origin_BDSKY_Contemp.t"


class ReproductiveNumberBDSC(Prior):
    def __init__(self):
        Prior.__init__(self)
        self.name = "reproductiveNumber"
        self.idx = "reproductiveNumberPrior_BDSKY_Contemp.t"
        self.x = "reproductiveNumber_BDSKY_Contemp.t"


class RhoBDSC(Prior):
    def __init__(self):
        Prior.__init__(self)
        self.name = 'rho'
        self.idx = "rhoPrior_BDSKY_Contemp.t"
        self.x = "rho_BDSKY_Contemp.t"


class BecomeUninfectiousBDSC(Prior):
    def __init__(self):
        Prior.__init__(self)
        self.name = 'becomeUninfectious'
        self.idx = "becomeUninfectiousRatePrior_BDSKY_Contemp.t"
        self.x = "becomeUninfectiousRate_BDSKY_Contemp.t"


# MultiType BirthDeath Priors


class ReproductiveNumberMTBD(Prior):
    def __init__(self):
        Prior.__init__(self)
        self.idx = "RPrior.t"
        self.x = "R0.t"


class SamplingProportionMTBD(Prior):
    def __init__(self):
        Prior.__init__(self)
        self.idx = "samplingProportionPrior.t"
        self.x = "samplingProportion.t"

    # Using a distribution component for prior here, not sure why:
    def get_prior_xml(self) -> str:

        dim, incl = self.get_include_string()

        return f"""
            <distribution id="{self.idx}:snp" spec="multitypetree.distributions.ExcludablePrior" x="@{self.x}:snp">
                <xInclude id="samplingProportionXInclude.t:snp" spec="parameter.BooleanParameter" dimension="{dim}">{incl}</xInclude>
                {self.distribution.get_xml()}
            </distribution>
        """

    def get_include_string(self) -> (int, str):

        incl = []
        for v in self.initial:
            if v != 0:
                incl.append('true')
            else:
                incl.append('false')

        return len(self.initial), " ".join(incl)


class BecomeUninfectiousMTBD(Prior):
    def __init__(self):
        Prior.__init__(self)
        self.idx = "becomeUninfectiousRatePrior.t"
        self.x = "becomeUninfectiousRate.t"


class RateMatrix(Prior):
    def __init__(self):
        Prior.__init__(self)
        self.idx = "rateMatrixPrior.t"
        self.x = "rateMatrix.t"


# Model priors for Coalescent Bayesian Skyline

class PopulationSize(Prior):
    def __init__(self):
        Prior.__init__(self)
        self.name = 'bPopSizes'
        self.idx = None  # no XML prior necessary
        self.x = "bPopSizes.t"


class GroupSize(Prior):
    def __init__(self):
        Prior.__init__(self)
        self.name = 'bGroupSizes'
        self.idx = None  # no XML prior necessary
        self.x = "bGroupSizes.t"

        self.param_spec = "parameter.IntegerParameter"

    def get_group_size_state_node(self):

        return f'<stateNode id="bGroupSizes.t:snp" spec="parameter.IntegerParameter" ' \
            f'dimension="{self.dimension}">{self.initial}</stateNode>'

# MultiType Priors


# Clock model base class

class Clock:

    def __init__(
            self, priors: [Prior]
    ):
        self.priors = priors
        self.fixed: bool or None = None
        self.state_node: str = ''  # defined in some clock subclasses
        self.branch_rate_model: str = self.get_branch_rate_xml()
        self.updown_operator: str = self.get_updown_operator_xml()

    # Combining prior XMLs

    def get_parameter_xml(self):
        return "\n".join([p.get_parameter_xml() for p in self.priors])

    def get_logger_xml(self):
        return "\n".join([p.get_logger_xml() for p in self.priors])

    def get_prior_xml(self):
        return "\n".join([p.get_prior_xml() for p in self.priors])

    def get_state_node_xml(self):
        return self.state_node

    # Defined in subclasses

    def get_scale_operator_xml(self) -> str:
        pass

    def get_branch_rate_xml(self) -> str:
        pass

    def get_updown_operator_xml(self) -> str:
        pass


# Clock model subclasses


class Strict(Clock):
    def __init__(self, priors):
        Clock.__init__(self, priors=priors)

    def get_branch_rate_xml(self) -> str:
        return f'<branchRateModel id="StrictClock.c:snp" spec="beast.evolution.branchratemodel.StrictClockModel" clock.rate="@clockRate.c:snp"/>'

    def get_scale_operator_xml(self):
        if self.fixed:
            return ""
        else:
            return f'<operator id="StrictClockRateScaler.c:snp" spec="ScaleOperator" parameter="@clockRate.c:snp" scaleFactor="0.5" weight="3.0"/>'

    def get_updown_operator_xml(self):
        if self.fixed:
            return ""
        else:
            return textwrap.dedent(f"""
                <operator id="strictClockUpDownOperator.c:snp" spec="UpDownOperator" scaleFactor="0.75" weight="3.0">
                    <up idref="clockRate.c:snp"/>
                    <down idref="Tree.t:snp"/>
                </operator>
            """)


class RelaxedExponential(Clock):
    def __init__(self, priors):
        Clock.__init__(self, priors=priors)
        self.state_node = f'<stateNode id="expRateCategories.c:snp" spec="parameter.IntegerParameter" dimension="718">1</stateNode>'

    def get_branch_rate_xml(self) -> str:
        return textwrap.dedent(f""" 
            <branchRateModel id="ExponentialRelaxedClock.c:snp" spec="beast.evolution.branchratemodel.UCRelaxedClockModel" clock.rate="@ucedMean.c:snp" rateCategories="@expRateCategories.c:snp" tree="@Tree.t:snp">
                <Exponential id="Exponential.c:snp" name="distr">
                    <parameter id="UCExpLambda.c:snp" spec="parameter.RealParameter" name="mean">1.0</parameter>
                </Exponential>
            </branchRateModel>
        """)

    def get_scale_operator_xml(self):
        if self.fixed:
            return ''
        else:
            return textwrap.dedent(f""" 
                <operator id="ucedMeanScaler.c:snp" spec="ScaleOperator" parameter="@ucedMean.c:snp" scaleFactor="0.5" weight="1.0"/>
                <operator id="ExpCategoriesRandomWalk.c:snp" spec="IntRandomWalkOperator" parameter="@expRateCategories.c:snp" weight="10.0" windowSize="1"/>
                <operator id="ExpCategoriesSwapOperator.c:snp" spec="SwapOperator" intparameter="@expRateCategories.c:snp" weight="10.0"/>
                <operator id="ExpCategoriesUniform.c:snp" spec="UniformOperator" parameter="@expRateCategories.c:snp" weight="10.0"/>
            """)

    def get_updown_operator_xml(self):
        if self.fixed:
            return ''
        else:
            return textwrap.dedent(f"""
                <operator id="relaxedUpDownOperatorExp.c:snp" spec="UpDownOperator" scaleFactor="0.75" weight="3.0">
                    <up idref="ucedMean.c:snp"/>
                    <down idref="Tree.t:snp"/>
                </operator>
            """)


class RelaxedLogNormal(Clock):
    def __init__(self, priors):
        Clock.__init__(self, priors=priors)
        self.state_node = f'<stateNode id="rateCategories.c:snp" spec="parameter.IntegerParameter" dimension="718">1</stateNode>'

    def get_branch_rate_xml(self) -> str:
        return textwrap.dedent(f""" 
            <branchRateModel id="RelaxedClock.c:snp" spec="beast.evolution.branchratemodel.UCRelaxedClockModel" clock.rate="@ucldMean.c:snp" rateCategories="@rateCategories.c:snp" tree="@Tree.t:snp">
                <LogNormal id="LogNormalDistributionModel.c:snp" S="@ucldStdev.c:snp" meanInRealSpace="true" name="distr">
                    <parameter id="RealParameter.{uuid.uuid4()}" spec="parameter.RealParameter" estimate="false" lower="0.0" name="M" upper="1.0">1.0</parameter>
                </LogNormal>
            </branchRateModel>
        """)

    def get_scale_operator_xml(self):
        if self.fixed:
            return ""
        else:
            return textwrap.dedent(f""" 
                <operator id="ucldMeanScaler.c:snp" spec="ScaleOperator" parameter="@ucldMean.c:snp" scaleFactor="0.5" weight="1.0"/>
                <operator id="ucldStdevScaler.c:snp" spec="ScaleOperator" parameter="@ucldStdev.c:snp" scaleFactor="0.5" weight="3.0"/>
                <operator id="CategoriesRandomWalk.c:snp" spec="IntRandomWalkOperator" parameter="@rateCategories.c:snp" weight="10.0" windowSize="1"/>
                <operator id="CategoriesSwapOperator.c:snp" spec="SwapOperator" intparameter="@rateCategories.c:snp" weight="10.0"/>
                <operator id="CategoriesUniform.c:snp" spec="UniformOperator" parameter="@rateCategories.c:snp" weight="10.0"/>
            """)

    def get_updown_operator_xml(self):
        if self.fixed:
            return ""
        else:
            return textwrap.dedent(f"""
                <operator id="relaxedUpDownOperator.c:snp" spec="UpDownOperator" scaleFactor="0.75" weight="3.0">
                    <up idref="ucldMean.c:snp"/>
                    <down idref="Tree.t:snp"/>
                </operator>
            """)
