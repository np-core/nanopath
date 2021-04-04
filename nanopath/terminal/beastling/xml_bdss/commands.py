import click

from pathlib import Path
from nanopath.beastling import BirthDeathSkylineSerial
from nanopath.utils import set_nested_item

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
    "--dimensions", "-dm", type=str, default=None,
    help="Use the same template configuration but produce multiple XML "
         "with dimensional slices of Re over a range(format: 1..10) [none]"
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
def xml_bdss(
    alignment, data, yaml, clock, mcmc, length, hot, outdir, intervals, prefix, sample_prior, dimensions, model_prior
):

    """ Pre-configured Birth-Death Skyline Serial XML """

    yaml_files = {prefix: yaml}
    outdir.mkdir(parents=True, exist_ok=True)

    for prefix, y in yaml_files.items():

        bdss = BirthDeathSkylineSerial(
            alignment=alignment,
            data=data,
            clock_model=clock,
            chain_type=mcmc,
            chain_length=length,
            chain_number=hot+1,
            prefix=prefix,
            sample_prior=sample_prior
        )

        bdss.print_configuration()
        bdss.check_configuration()

        config = bdss.read_config(file=y)

        # Set model prior configuration
        model_priors = config.get('priors').get('model')

        # Modify the model prior configs if settings are passed

        print(model_priors)
        
        for mp in model_prior:
            mp_nest = mp.split(":")
            mp_path, mp_val = mp_nest[:-1], mp_nest[-1]
            print(mp_path, mp_val, mp_nest)
            model_priors = set_nested_item(model_priors, mp_path, mp_val)

        print(model_priors)

        bdss.set_model_priors(prior_config=model_priors, distribution=True)

        # Set clock prior configuration
        clock_priors = config.get('priors').get('clock')
        bdss.set_clock(prior_config=clock_priors)

        if intervals:
            # Set slice configurations and overwrite associated priors
            slice_config = config.get('priors').get('intervals')
            bdss.set_slices(slice_config=slice_config)

        if dimensions:
            df, dt = dimensions.split("..")
            drange = range(int(df), int(dt)+1)

            for d in drange:
                bdss.construct_template(
                    xml=outdir / f'{prefix}_d{d}.xml', r_dimension=d
                )
        else:
            bdss.construct_template(
                xml=outdir / f'{prefix}.xml'
            )