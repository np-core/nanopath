import click
import pandas
from nanopath.pipelines import PathogenPipeline
from pathlib import Path


@click.command()
@click.option(
    '--path',
    '-p',
    type=Path,
    help='Path to pipeline output directory containing read files for assembly'
)
@click.option(
    '--outdir',
    '-o',
    type=Path,
    default=Path('rave'),
    help='Path to output directory [$PWD/rave]'
)
@click.option(
    '--file_glob',
    '-f',
    type=str,
    default="*.fastq",
    help='File glob to search for input reads files [*.fastq]'
)
@click.option(
    '--threads',
    '-t',
    type=int,
    default=None,
    help='Parallel processes to launch Raven assemblies of '
         'grouped files for total of: parallel*raven_cpu load[4]'
)
@click.option(
    '--regex',
    '-rg',
    type=str,
    default=r"barcode\d\d",
    help=r'Regex to group files for assembly by [barcode\d\d]'
)
@click.option(
    '--raven_cpu',
    '-rc',
    type=int,
    default=2,
    help='Number of processors to launch Raven with [2]'
)
@click.option(
    '--raven_args',
    '-ra',
    type=str,
    default="",
    help='Arguments to be passed to Raven assembler [none]'
)
def raven_mag(path, outdir, file_glob, regex, raven_cpu, raven_args, threads):

    """ Group files by regex and execute Raven on each group (for pipeline) """

    regex = fr"{regex}"

    if raven_args is None:
        raven_args = ""

    pp = PathogenPipeline(workdir=outdir)

    pp.logger.info('Welcome to fast metagenome assembly with Raven')

    # get regexed file groups
    files = [p for p in path.glob(file_glob) if p.is_file()]
    pp.logger.info(
        f'Detected {len(files)} read files for processing at: {path}'
    )
    grouped_files = pp.group_files_by(files=files, regex=regex)
    pp.logger.info(
        f'Detected {len(grouped_files)} regex ({regex}) groups:'
    )

    for key in grouped_files.keys():
        pp.logger.info(f'Regex group: {key}')

    if threads is None:
        for data in [
            list(data) + [raven_cpu, raven_args]
            for data in grouped_files.items()
        ]:
            pp._assemble_raven(group_data=data)
    else:
        pp.assemble(
            grouped_files=grouped_files,
            assembler="raven",
            assembler_cpu=raven_cpu,
            threads=threads,
            raven_args=raven_args
        )

    # pp.clean()  # don't remove workdir



