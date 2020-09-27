import click
import yaml

from pathlib import Path

# example
# data/models/templates/bdss/Exploration
# for ST in ST*.fasta; do st=$(basename ${ST} .fasta);
#   np beastling config -o bur_high -b ${st}.yaml --priors prior_configs/bur_high_prior.yaml
# done


@click.command()
@click.option(
    "--outdir", "-o", required=True, type=Path, default=None,
    help="Output directory for config files",
)
@click.option(
    "--model", "-m", required=False, type=str, default="bdss",
    help="Model to select configuration template: bdss, bdsc, mtbd, cosky",
)
@click.option(
    "--template", "-t", required=False, type=Path, default=Path("bdss.yaml"),
    help="Output of template configuration file for modification (YAML)",
)
@click.option(
    "--base", "-b", required=False, type=str, default=None,
    help="List of case configuration file, comma-delimited (YAML)",
)
@click.option(
    "--priors", "-p", required=False, type=Path, default=None,
    help="Prior configuration file to add to base; entries overwrite base if present (YAML)",
)
@click.option(
    "--dim", "-d", required=False, type=str, default="5,6,7,8,9,10",
    help="Comma delimited string of dimensions along which to modify the config files",
)
@click.option(
    "--dim_prior", "-dp", required=False, type=str, default="reproductive_number",
    help="Name of the prior to set the dimension string files to [reproductive_number]",
)
def config(model, template, outdir, base, priors, dim, dim_prior):

    """ Get a Beastling configuration template YAML """

    dimensions = [int(d) for d in dim.split(",")]
    base_configs = [Path(b) for b in base.split(",")]

    outdir.mkdir(parents=True, exist_ok=True)

    for b in base_configs:
        print(f"Creating config for base config: {b}")
        # For each base file generate the config
        with b.open("r") as base_config:
            bconf = yaml.load(base_config, Loader=yaml.FullLoader)

        print(f"Using prior configuration file: {priors}")
        with priors.open("r") as prior_config:
            pconf = yaml.load(prior_config, Loader=yaml.FullLoader)

        for prior_conf, settings in pconf.items():  # model, clock, intervals
            if prior_conf in bconf['priors'].keys():
                print(f"Found prior configuration: {prior_conf}")
                for prior in bconf['priors'][prior_conf].keys():
                    if prior in settings.keys():
                        print(f"Replacing prior {prior} with user settings")
                        bconf['priors'][prior_conf][prior] = settings[prior]

        bconf_dim = {}
        for d in dimensions:
            print(f"Creating configuration file variant for dimension {d} of model prior: {dim_prior}")
            dconf = bconf.copy()
            try:
                dconf['priors']['model'][dim_prior]['dimension'] = int(d)
            except KeyError:
                print(
                      f'Could not find dimension entry for '
                      f'model prior {dim_prior} at priors/model/{dim_prior}/dimension'
                )
                continue
            except TypeError:
                print('Something went wrong with the dimension integer list input')
                exit(1)
            finally:
                bconf_dim[d] = dconf
                fname = f"{b.stem}-{priors.stem}-{dim_prior}-d{d}.yaml"
                with (outdir / fname).open("w") as dconf_file:
                    yaml.dump(dconf, dconf_file)

                print(
                    f"Wrote configuration file variant for dimension {d} of"
                    f" model prior: {dim_prior} to {outdir / fname}"
                )

