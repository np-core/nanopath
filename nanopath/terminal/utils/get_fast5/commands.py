import click
import pandas

from pyfastaq import sequences
from pathlib import Path
from nanopath.utils import run_cmd


@click.command()
@click.option(
    '--fastq',
    '-fq',
    type=Path,
    help='Path to sequence read files to get Fast5 [non-recursive]'
)
@click.option(
    '--fast5',
    '-f5',
    type=Path,
    help='Path to Fast5 files [recursive]'
)
@click.option(
    '--outdir',
    '-o',
    type=Path,
    default=Path.cwd() / 'fast5',
    help='Path to output directory with --fastq subdirectories containing Fast5'
)
@click.option(
    '--batch_size',
    '-b',
    type=int,
    default=4000,
    help='Path to output directory with --fastq subdirectories containing Fast5'
)
@click.option(
    '--extension',
    '-e',
    type=str,
    default=".fastq",
    help='Path to file with one name per line of isolates to include'
)
@click.option(
    '--subset',
    '-s',
    type=Path,
    default=None,
    help='Path to file with one file basename per line of isolate to include'
)
@click.option(
    '--threads',
    '-t',
    type=int,
    default=12,
    help='Threads to use for fetching Fast5'
)
def get_fast5(fastq, fast5, subset, batch_size, outdir, extension, threads):

    """ Get Fast5 reads from basecalled Fastq """

    outdir.mkdir(parents=True, exist_ok=True)
    fq_files = fastq.glob(f"*{extension}")

    if subset:
        df = pandas.read_csv(subset, sep='\t')
        names = df.iloc[:, 0].tolist()
    else:
        names = []

    for fq in fq_files:

        name = fq.stem

        # Subset checks
        if names:
            if name in names:
                pass
            else:
                print(f"{name} not in subset, ignoring")
                continue

        # Create files:

        read_ids = [seq.id for seq in sequences.file_reader(str(fq))]
        (outdir / name).mkdir(parents=True, exist_ok=True)

        read_id_list = outdir / name / f"{name}.txt"
        with open(read_id_list, 'w') as outfile:
            for read_id in read_ids:
                outfile.write(read_id + '\n')

        # Run ONT Fast5 API:

        print(f"Fetching Fast5 for: {fq}")

        run_cmd(
            f"fast5_subset --input {fast5} --save_path {outdir / name} --read_id_list {read_id_list} "
            f"--batch_size {batch_size} --recursive --filename_base {name} --threads {threads}"
        )



