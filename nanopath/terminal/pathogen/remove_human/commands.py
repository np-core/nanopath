import click
import logging

from pathlib import Path
from nanopath.utils import smart_open, PoreLogger

try:
    import mappy as mp
except ModuleNotFoundError:
    pass  # windows


@click.command()
@click.option(
    '--fq',
    type=Path,
    help='Reads Oxford Nanopore or PacBio'
)
@click.option(
    '--fq1',
    type=Path,
    help='Forward mate pair PE Illumina [none]'
)
@click.option(
    '--fq2',
    type=Path,
    help='Reverse mate pair PE Illumina [none]'
)
@click.option(
    '--ref',
    type=Path,
    help='Minimap2 indexed reference genome index to map against [none]'
)
@click.option(
    '--out',
    type=str,
    default="-",
    help='Write reads to out [-]'
)
@click.option(
    '--mapq',
    type=int,
    default=30,
    help='Minimum mapping quality [30]'
)
@click.option(
    '--preset',
    type=str,
    default="map-ont",
    help='Mapping configuration preset for minimap2 [map-ont]'
)
@click.option(
    '--human_out',
    type=Path,
    help='File for human reads to write (slower) [none]'
)
@click.option(
    '--threads',
    type=int,
    default=2,
    help='Parallel processes for index parsing  [none]'
)
def remove_human(fq, fq1, fq2, ref, out, mapq, preset, human_out, threads):

    """ Filter a read file (single, paired) by mapping to human reference index """

    logger = PoreLogger(level=logging.INFO, name="HostRemoval").logger

    remove_by_alignment(
        fq=fq, ref=ref, out=out,
        mapq=mapq, preset=preset,
        human_out=human_out, threads=threads, logger=logger
    )


def remove_by_alignment(
    fq, ref, out, mapq, preset, human_out, threads, logger
):

    fout = smart_open(filename=out, mode="w")

    if human_out:
        hout = smart_open(filename=human_out, mode="w")
    else:
        hout = None

    logger.info(f"Starting to map reads against: {ref}")

    logger.info(f"Initiating aligner: {ref}")
    aligner = mp.Aligner(str(ref), preset=preset, n_threads=threads)

    logger.info(f"Opening file handle: {fq}")
    if fq:
        reads = mp.fastx_read(str(fq))
    else:
        reads = None  # PE

    ref_maps = 0
    total_reads = 0

    logger.info(f"Filtering mapped reads [Q >= {mapq}]")

    human = []
    not_human = []
    for name, seq, qual in reads:
        mapped = aligner.map(seq)
        for aln in mapped:
            if aln.mapq >= mapq:
                ref_maps += 1
                if name not in human:
                    human.append(name)
                if hout is not None:
                    hout.write(
                        str(f"@{name}\n{seq}\n+\n{qual}\n")
                    )
                continue

        if name not in human:
            fout.write(
                str(f"@{name}\n{seq}\n+\n{qual}\n")
            )
            if name not in not_human:
                not_human.append(name)

        total_reads += 1

    fout.close()

    if hout is not None:
        hout.close()

    logger.info(f"Computed {ref_maps} mappings against reference: {ref}")
    logger.info(f"Recovered  {len(not_human)} / {total_reads} reads from {fq}")