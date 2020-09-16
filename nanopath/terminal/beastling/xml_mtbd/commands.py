import click

from pathlib import Path
from nanopath.beastling import MultiTypeBirthDeath


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
    "--prefix", "-p", required=False, type=str, default='bdss',
    help="Prefix for sample logs from BEAST [bdss]"
)
@click.option(
    "--trait", "-t", required=False, type=str, default='trait',
    help="Trait column name to use for demes in model [trait]"
)
def xml_mtbd(alignment, data, yaml, clock, mcmc, length, hot, prefix, trait):

    """ Pre-configured Multi-Type Birth-Death XML """

    mtbd = MultiTypeBirthDeath(
        alignment=alignment,
        data=data,
        clock_model=clock,
        chain_type=mcmc,
        chain_length=length,
        chain_number=hot+1,
        prefix=prefix,
        trait=trait
    )

    mtbd.print_configuration()
    mtbd.check_configuration()

    config = mtbd.read_config(file=yaml)

    # Set model configuration
    model_config = config.get('model')
    mtbd.set_model_config(model_config=model_config)

    # Set model prior configuration
    model_priors = config.get('priors').get('model')
    mtbd.set_model_priors(prior_config=model_priors, distribution=True)

    # Set clock prior configuration
    clock_priors = config.get('priors').get('clock')
    mtbd.set_clock(prior_config=clock_priors)

    # Check deme settings
    mtbd.check_deme_config()

    mtbd.construct_template(
        xml=Path(f'{prefix}.xml')
    )
