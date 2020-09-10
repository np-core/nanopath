from pathlib import Path
import click
import json
import pandas
import os

from nanopath.surveillance import BioSampler


@click.command()
@click.option(
    '--file', '-f', type=Path, required=False,
    help='Data table to extract accessions from'
)
@click.option(
    '--outfile', '-o', type=Path, required=False, default="[queries.json]",
    help='Output response from queries by query key JSON format [queries.json]'
)
@click.option(
    '--query_column',  type=str, required=False, default="sample",
    help='Column in data containing accessions to query'
)
@click.option(
    '--json_file',  type=Path, required=False, default=None,
    help='JSON response file from queries to extract values for'
)
@click.option(
    '--link',  type=Path, required=False, default=None,
    help='Path to symlink sequence read files from [none]'
)
@click.option(
    '--link_column', type=str, default='accession',
    help='Column name to link files from link folder [none]'
)
@click.option(
    '--link_outdir',  type=Path, required=False, default=Path.cwd(),
    help='Path to symlink sequence read files to [none]'
)
@click.option(
    '--prefix',  type=str, required=False, default="complete",
    help='Prefix for output data [complete]'
)
def query(file, outfile, query_column, link_column, json_file, link, link_outdir, prefix):

    """ Query and process metaddata from the BioSample DB """

    bs = BioSampler()

    if json_file:
        with json_file.open() as infile:
            query_responses = json.load(infile)

        collection_dates = []
        total = 0
        locations = []
        combined = []
        for query, response in query_responses.items():
            characteristics = response['characteristics']
            try:
                collection_date = characteristics['collection date']
                collection_dates.append(collection_date)
                # print(f'Found collection date for {query}')
            except KeyError:
                # print(f'Could not find collection date for {query}')
                collection_date = None

            try:
                location1 = characteristics['geographic location (country and/or sea)']
            except KeyError:
                location1 = None
                pass

            try:
                location2 = characteristics['geo loc name']
            except KeyError:
                location2 = None
                pass

            if location1 or location2:
                if location1 is not None:
                    locations.append(location1)
                else:
                    locations.append(location2)

            if (location1 or location2) and collection_date:
                combo = {}
                if location1 is not None:
                    combo["location"] = location1[0]['text']
                else:
                    combo["location"] = location2[0]['text']

                combo['date'] = collection_date[0]['text']
                combo['query'] = query
                combined.append(combo)

            total += 1

            # print(query)

        print(f'Found {len(collection_dates)} dates out of {total} queries')
        print(f'Found {len(locations)} locations out of {total} queries')
        print(f'Found {len(combined)} combined out of {total} queries')

        complete = pandas.DataFrame(combined)

        if file:
            df = pandas.read_csv(file, sep='\t')

            df = df.merge(
                complete,
                left_on=query_column,
                right_on='query',
                how='inner'
            )

            if link:
                link_outdir.mkdir(exist_ok=True, parents=True)

                files = link.glob("*")
                link_names = tuple(
                    df[link_column].tolist()
                )
                detected = 0
                for file in files:
                    if file.is_file() and file.name.startswith(link_names):
                        try:
                            print(f'Linking file: {file}')
                            os.symlink(
                                str(file.absolute()),
                                str((link_outdir / file.name).absolute())
                            )
                        except FileExistsError:
                            print(f'Detected file: {file}')
                        detected += 1

                print(
                    f'Linked {detected//2}/{len(df)} metadata-complete files'
                )

            df.to_csv(f'{prefix}.tsv', sep='\t', index=None)

    else:
        bs.process_list(
            file=file,
            sep='\t',
            query_column=query_column,
            outfile=outfile
        )

    #     df = pandas.read_csv(file, sep=delimiter, header=0)
    #     df.columns = [c.lower() for c in df.columns]
    #     if 'accession' in df:
    #         accession = df['accession'].tolist()
    #         project, species, query = (None,) * 3
    #     elif 'project' in df:
    #         project = df['project']
    #         accession, species, query = (None,) * 3
    #     else:
    #         print('Could not find columns: accession or project in CSV.')
    #
    #
    # if not no_check_dir:
    #     Path(outdir).mkdir(exist_ok=True, parents=True)
    #
    # survey = Survey(outdir=outdir)

    # results = survey.query_ena(
    #     sample=accession,
    #     study=project,
    #     species=species,
    #     term=query,
    #     scheme=scheme,
    #     submitted_fastq=submitted,
    #     allow_merged_fastq=allow_merged,
    #     allow_missing=allow_missing
    # )
    #
    # if query_only:
    #     rdf = pandas.concat(
    #         [data['results'] for i, data in results.items()]
    #     )
    #     rdf.to_csv(f'{outdir}/{prefix}.tsv', sep='\t', index=False)
    #     exit(0)