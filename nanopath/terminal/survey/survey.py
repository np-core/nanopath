import click

from .filter import filter
from .collect import collect
from .download import download
from .query import query
from .find import find

VERSION = '0.1'


@click.group()
@click.version_option(version=VERSION)
def survey():

    """ Survey: genome surveillance and processing from ENA """

    pass


survey.add_command(query)
survey.add_command(collect)
survey.add_command(download)
survey.add_command(filter)
survey.add_command(find)