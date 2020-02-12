import click
from pathlib import Path
from nanopath.pipelines import Metagenome


@click.command()
@click.option(
    '--report', '-r', type=str, help='Path or file glob to tax report files'
)
@click.option(
    '--prefix', '-p', type=str, help='Output prefix for plot file.'
)
@click.option(
    '--top', '-t', type=int, default=10,
    help='Show top taxonomic levels in plots [10]'
)
@click.option(
    '--color', '-c', type=str, default='Greens',
    help='Color palette for central donut plot.'
)
@click.option(
    '--title', '-t', type=str, default=None,
    help='Row titles for center plot, comma separated string.'
)
@click.option(
    '--sub', '-s', is_flag=True,
    help='Add subplot titles for each column.'
)
def process_kraken(report, prefix, top, color, title, sub):

    """ Process results and generate server data in the metagenome pipeline  """

    mg = Metagenome()
    mg.set_meta(
        uuid=None,
        user=None,
        name=None,
        submitted=None,
        host_removal=''
    )
    mg.process_kraken(report=report)
    mg.get_server_data(
        fout=Path(f'{prefix}.json')
    )
    mg.plot_kraken_summary(
        plot_file=Path(f'{prefix}.png'),
        palette=color,
        top_minor=top,
        title=title,
        subtitles=sub
    )
