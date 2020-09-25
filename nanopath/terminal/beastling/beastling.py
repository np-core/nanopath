import click

from .config import config
from .xml_bdsc import xml_bdsc
from .xml_bdss import xml_bdss
from .xml_cosky import xml_cosky
from .xml_mtbd import xml_mtbd

VERSION = '1'


@click.group()
@click.version_option(version=VERSION)
def beastling():

    """ Beastling: pre-configured transmission models for bacteria in BEAST2 """

    pass


beastling.add_command(config)
beastling.add_command(xml_mtbd)
beastling.add_command(xml_bdss)
beastling.add_command(xml_bdsc)
beastling.add_command(xml_cosky)
