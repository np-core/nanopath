import click

from pathlib import Path
from nanopath.beastling import BirthDeathSkylineContemporary


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
    "--intervals", "-i", is_flag=True,
    help="Use sampling proportion intervals from configuration YAML [false]"
)
@click.option(
    "--prefix", "-p", required=False, type=str, default='bdss',
    help="Prefix for sample logs from BEAST [bdss]"
)
def xml_bdsc(alignment, data, yaml, clock, mcmc, length, hot, intervals, prefix):

    """ Pre-configured Birth-Death Skyline Contemporary XML """

    bdsc = BirthDeathSkylineContemporary(
        alignment=alignment,
        data=data,
        clock_model=clock,
        chain_type=mcmc,
        chain_length=length,
        chain_number=hot+1,
        prefix=prefix
    )

    bdsc.print_configuration()
    bdsc.check_configuration()

    config = bdsc.read_config(file=yaml)

    # Set model prior configuration
    model_priors = config.get('priors').get('model')
    bdsc.set_model_priors(prior_config=model_priors, distribution=True)

    # Set clock prior configuration
    clock_priors = config.get('priors').get('clock')
    bdsc.set_clock(prior_config=clock_priors)

    if intervals:
        # Set slice configurations and overwrite associated priors
        slice_config = config.get('priors').get('intervals')
        bdsc.set_slices(slice_config=slice_config)

    bdsc.construct_template(
        xml=Path(f'{prefix}.xml')
    )
