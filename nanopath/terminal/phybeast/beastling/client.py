import click

VERSION = '0.1'


from .xml_bdsc import xml_bdsc
from .xml_bdss import xml_bdss
from .xml_cosky import xml_cosky
from .xml_mtbd import xml_mtbd


@click.group()
@click.version_option(version=VERSION)
def beastling():

    """ Beastling: pre-configured transmission models for bacteria in BEAST """

    pass


beastling.add_command(xml_mtbd)
beastling.add_command(xml_bdss)
beastling.add_command(xml_bdsc)
beastling.add_command(xml_cosky)
