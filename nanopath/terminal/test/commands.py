import click

from pathlib import Path
from nanopath.pipelines import NanoPathLive


@click.command()
@click.option(
    '--directory',
    '-d',
    type=Path,
    help=''
)
@click.option(
    '--outdir',
    '-o',
    type=Path,
    help=''
)
def test(directory, outdir):

    """ Task for testing command line functions """

    np = NanoPathLive(path=directory)

    online_data = np.update_online_view()
    telemetry_data, _ = np.update_telemetry_view()
    run_view_data, read_lengths, read_qualities = np.update_run_view(
        length_bin_min=0,
        length_bin_max=80000,
        length_bin_size=2000,
        quality_bin_min=0.0,
        quality_bin_max=21.0,
        quality_bin_size=1.0
    )
    print(online_data)
    print(telemetry_data)
    print(run_view_data)
    print(read_lengths)
    print(read_qualities)