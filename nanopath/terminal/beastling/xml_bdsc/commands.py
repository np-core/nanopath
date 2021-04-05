import click

from pathlib import Path
from nanopath.beastling import BirthDeathSkylineContemporary
from nanopath.utils import modify_model_priors

@click.command()
@click.option(
    "--alignment", "-a", required=True, type=Path,
    help="Variable site alignment of non-recombinant core-genome SNPs",
)
@click.option(
    "--data", "-d", required=True, type=Path,
    help="Data file for samples in the alignment with headers: name, date"
)
@click.option(
    "--yaml", "-y", required=False, type=Path, default=None,
    help="YAML configuration file including prior settings [None]"
)
@click.option(
    "--yaml_dir", "-yd", required=False, type=Path, default=None,
    help="YAML configuration file directory to process all YAML files from [None]"
)
@click.option(
    "--yaml_glob", "-yg", required=False, type=str, default="*.yaml",
    help="Glob on the directory to subset config files matching alignment and data input ['*.yaml']"
)
@click.option(
    "--clock", "-c", required=False, type=str, default='strict',
    help="Molecular clock model: strict, relaxed_exponential, relaxed_lognormal [strict]"
)
@click.option(
    "--mcmc", "-m", required=False, type=str, default='default',
    help="Chain type: MCMC (default) or Coupled MCMC (coupled) [default]"
)
@click.option(
    "--length", "-l", required=False, type=int, default=1e7,
    help="Number of steps in Monte Carlo chain [1e7]"
)
@click.option(
    "--hot", "-h", required=False, type=int, default=3,
    help="Number of hot chains in Coupled Monte Carlo chain [3]"
)
@click.option(
    "--intervals", "-i", is_flag=True,
    help="Use sampling proportion intervals from configuration YAML [false]"
)
@click.option(
    "--sample_prior", "-s", is_flag=True,
    help="Sample from prior [false]"
)
@click.option(
    "--prefix", "-p", required=False, type=str, default='bdss',
    help="Prefix for sample logs from BEAST [bdss]"
)
@click.option(
    "--outdir", "-o", required=False, type=Path, default=Path('bdss'),
    help="Outdir for XML files [$PWD/bdss]"
)
@click.option(
    "--model_prior", "-pr", required=False, type=str, default=None, multiple=True,
    help="One or multiple args setting the replacement prior value in the YAML file with keys in string [:]"
)
@click.option(
    "--tag", "-t", is_flag=True,
    help="If modified on the fly attach the key path and setting to the output prefix [false]"
)
def xml_bdsc(alignment, data, outdir, yaml, yaml_dir, yaml_glob, clock, mcmc, length, hot, intervals, prefix, sample_prior, model_prior, tag):

    """ Pre-configured Birth-Death Skyline Contemporary XML """

    if yaml_dir is not None:
        yaml_files = {f.stem: f for f in yaml_dir.glob(f"{yaml_glob}")}
    else:
        yaml_files = {prefix: yaml}

    outdir.mkdir(parents=True, exist_ok=True)

    for prefix, y in yaml_files.items():

        bdsc = BirthDeathSkylineContemporary(
            alignment=alignment,
            data=data,
            clock_model=clock,
            chain_type=mcmc,
            chain_length=length,
            chain_number=hot+1,
            prefix=prefix,
            sample_prior=sample_prior
        )

        bdsc.print_configuration()
        bdsc.check_configuration()

        config = bdsc.read_config(file=yaml)

        # Set model prior configuration
        model_priors = config.get('priors').get('model')
        # Modify the model prior configs if settings are passed
        if model_prior:
            model_priors = modify_model_priors(model_priors, model_prior, tag, prefix)

        bdsc.set_model_priors(prior_config=model_priors, distribution=True)

        # Set clock prior configuration
        clock_priors = config.get('priors').get('clock')
        bdsc.set_clock(prior_config=clock_priors)

        if intervals:
            # Set slice configurations and overwrite associated priors
            slice_config = config.get('priors').get('intervals')
            bdsc.set_slices(slice_config=slice_config)

        bdsc.construct_template(
            xml=outdir / f'{prefix}.xml'
        )
