import click
import pandas
import networkx as nx

from tqdm import tqdm
from numpy import array
from numpy import fill_diagonal, count_nonzero, where
from pathlib import Path


import matplotlib.pyplot as plt

@click.command()
@click.option(
    '--alignment',
    type=Path,
    help='Path to SNP alignment (.fasta) to compute outbreak distances from'
)
@click.option(
    '--matrix',
    type=Path,
    help='Path to tab-delimited pairwise SNP difference matrix with rownames [none]'
)
@click.option(
    '--threshold',
    type=int,
    default=4,
    help='SNP threshold to apply for connecting recent potential transmission events [4]'
)
@click.option(
    '--prefix',
    type=str,
    default='test',
    help='Output prefix for path to output files [test]'
)
@click.option(
    '--substitution_rate',
    type=float,
    help='Substitution rate if used in combination with genome size computes SNP threshold [none]'
)
@click.option(
    '--genome_size',
    type=int,
    help='Genome size from which SNPs are estimated to compute threshold with substitution rate [none]'
)
@click.option(
    '--data',
    type=Path,
    help='Data frame with header and group assignments of samples in alignment'
         ' to compute threshold statistics for [none]'
)
@click.option(
    '--column',
    type=str,
    default="longitudinal_color",
    help='Comma-separated column names to compute group statistics for [none]'
)
def snp_distance(
    alignment,
    threshold,
    matrix,
    substitution_rate,
    genome_size,
    prefix,
    data,
    column
):

    """ Compute SNP outbreak distance based on substitution thresholds  """

    if substitution_rate and genome_size:
        threshold = int(substitution_rate*genome_size)

    print(f'SNP threshold set to: {threshold}')

    if matrix:
        snp_dist = pandas.read_csv(matrix, sep='\t', index_col=0, header=None)
        samples = snp_dist.index.tolist()
    else:
        genotypes, samples = read_fasta_alignment(alignment=alignment)

        print(
            f'Computing pairwise difference between sites in alignment'
        )
        snp_dist = pairwise_snp_difference(genotypes=genotypes)
        snp_dist = pandas.DataFrame(snp_dist, index=samples)

        snp_dist.to_csv(
            f'{prefix}_snp_dist.tsv', sep='\t', index=True, header=None
        )

    snp_dist.columns = samples
    snp_dist_matrix = snp_dist.values

    snp_dist_matrix[snp_dist_matrix <= threshold] = 1
    snp_dist_matrix[snp_dist_matrix > threshold] = 0
    fill_diagonal(snp_dist_matrix, 0)

    adjacency_matrix = snp_dist_matrix.copy()

    print(
        f'Detected {count_nonzero(adjacency_matrix) // 2} events'
    )

    graph = nx.from_numpy_array(adjacency_matrix)

    graph.remove_nodes_from(
        list(nx.isolates(graph))
    )
    fig, ax = plt.subplots(
        nrows=1, ncols=1, figsize=(
            1 * 7, 1 * 4.5
        )
    )
    pos = nx.spring_layout(graph)

    if data:
        data = pandas.read_csv(data, sep='\t', header=0)

        # Remove nodes in longitudinal PNG samples of same patient
        nodes = [node for node in graph.nodes]

        connections = []
        removed_longitudinal = 0
        for node in nodes:
            name = samples[node]
            is_longitudinal = get_value_from_data(data, name, 'longitudinal')

            if is_longitudinal == 'True':
                graph.remove_node(node)
                removed_longitudinal += 1

        print(f'Removed {removed_longitudinal} longitudinal samples from graph')

        node_names = [samples[node] for node in graph.nodes]

        node_colors = [
            get_value_from_data(data, name, column) for name in node_names
        ]

        graph_data = data[data["Laboratory_ID"].isin(node_names)]
        #
        # for node in graph.nodes:
        #     edges = graph.edges(node)
        #     print(edges)

        state_counts_total = data['state'].value_counts()
        for name, gd in graph_data.groupby("state"):
            print(name, len(gd), state_counts_total[name])
            print(round(len(gd) / state_counts_total[name], 4))

    else:
        node_colors = "gray"

    print(
        f'Detected {len(graph.nodes)} out of {len(adjacency_matrix)} '
        f'samples in putative recent transmission.'
    )

    nx.draw(
        graph,
        pos,
        node_size=8,
        node_color=node_colors,
        width=1,
        with_labels=False,
        ax=ax
    )

    fig.subplots_adjust(hspace=0.8)

    plt.tight_layout()
    fig.savefig(f"{prefix}_graph.pdf")


def get_edges_for_nodes(nodes: list, graph: nx.Graph):

    pass


def get_value_from_data(data, name, column):

    return data.at[
        where(
            data['Laboratory_ID'] == name
        )[0][0], column
    ]


def pairwise_snp_difference(genotypes: list) -> array:

    matrix = []
    for i, g1 in tqdm(
        enumerate(genotypes), total=len(genotypes)
    ):
        row = []
        for j, g2 in enumerate(genotypes):
            row.append(
                sum(a1 != a2 for a1, a2 in zip(g1, g2))
            )
        matrix.append(row)

    return array(matrix)


def read_fasta_alignment(alignment: Path) -> (list, list):

    samples = []
    genotypes = []
    with alignment.open('r') as aln:
        for line in aln:
            if line.startswith('>'):
                name = line.strip().replace('>', '')
            else:
                genotypes.append(
                    line.strip()
                )
                samples.append(name)
    if len(samples) != len(genotypes):
        raise ValueError

    return genotypes, samples

