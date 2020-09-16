import click

from .list import list
from .create import create
from .store import store
from .delete import delete
from .inspect import inspect
from .retrieve import retrieve

VERSION = "1"


@click.group()
@click.version_option(version=VERSION)
def dataset():

    """ Hybrid genome assembly and pangenome pipeline in Nextflow """

    pass


dataset.add_command(inspect)
dataset.add_command(list)
dataset.add_command(delete)
dataset.add_command(store)
dataset.add_command(create)
dataset.add_command(retrieve)