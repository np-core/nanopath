import click
from pathlib import Path
from nanopath.netflow import NanoPathClient, ClientError


@click.command()
@click.option(
    '--workdir', '-w',
    type=Path,
    required=True,
    help='Path to working directory'
)
@click.option(
    '--server_config', '-s',
    type=Path,
    required=True,
    help='Path to server and pipeline configuration file'
)
@click.option(
    '--fastq', '-fq',
    type=Path,
    help='Path to search for basecalled nanopore reads (.fq / .fastq)'
)
@click.option(
    '--illumina',
    is_flag=True,
    help='When activated, input reads are Illumina PE (.fq.gz / .fastq.gz)'
)
@click.option(
    '--recursive', '-r',
    is_flag=True,
    help='Search path for read files recursively'
)
@click.option(
    '--keep_open',
    is_flag=True,
    help='When activated, pipeline screens will be open and client will run indefinitely'
)
@click.option(
    '--online',
    is_flag=True,
    help='When activated, launches watchdogs and online monitoring of run'
)
@click.option(
    '--configure',
    is_flag=True,
    help='When activated, configure servers and conduct pre-launch checks only'
)
@click.option(
    '--nuke',
    is_flag=True,
    help='When activated, remove all working directories including '
         'on remote after exiting client'
)
@click.option(
    '--resume',
    is_flag=True,
    help='When activated, add the -resume flag to pipelines, '
         'requires existing working directory.'
)
def run(
    server_config, configure, online, fastq, workdir, nuke,
    illumina, keep_open, resume, recursive
):

    """ Initiate server connections and run pipelines  """

    if fastq is None:
        raise ClientError(
            'Either one or both inputs must be specified (--fast5 | --fastq)'
        )

    if resume and not workdir.exists():
        raise ClientError(
            f'If resuming operations in pipelines, '
            f'you need the existing working directory: {workdir}'
        )

    if nuke:
        click.confirm(
            'All working directories will be removed after exiting the program'
            'Do you want to continue?', abort=True
        )

    if keep_open:
        click.confirm(
            'Client will run indefinitely and requires user input to exit.'
            'Do you want to continue?', abort=True
        )

    cli = NanoPathClient(server_config=server_config, workdir=workdir)

    if configure:
        cli.configure_clients()
        exit(0)

    if online:
        pass
    else:
        cli.configure_clients()
        cli.launch_pipelines(
            fastq=fastq, illumina=illumina, recursive=recursive,
            nuke=nuke, keep_open=keep_open, resume=resume
        )
