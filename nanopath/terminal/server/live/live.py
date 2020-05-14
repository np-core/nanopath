
import click

from .collect import collect

VERSION = "0.1"


@click.group()
@click.version_option(version=VERSION)
def live():
    """ Live dashboard command line tasks for pipeline support """
    pass

import click

from .collect import collect

VERSION = "0.1"


@click.group()
@click.version_option(version=VERSION)
def live():
    """ Live dashboard command line tasks for pipeline support """
    pass

live.add_command(collect)