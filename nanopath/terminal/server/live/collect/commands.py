import click

from nanopath.pipelines import NanoPathLive
from pathlib import Path


@click.command()
@click.option(
    '--path',
    '-p',
    type=Path,
    help='Path to pipeline output directory containing results'
)
def collect(path):
    """ Collect output from sepsis pipeline into server data """

    npl = NanoPathLive(path=path)

    online_data = npl.update_stats_view()
    telemetry_data, _ = npl.update_telemetry_view()
    run_view_data, read_lengths, read_qualities = npl.update_run_view(
        length_bin_min=0,
        length_bin_max=80000,
        length_bin_size=2000,
        quality_bin_min=0.0,
        quality_bin_max=21.0,
        quality_bin_size=1.0
    )

