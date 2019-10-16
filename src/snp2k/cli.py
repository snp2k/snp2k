# -*- coding: utf-8 -*-

"""Command line interface for SNP2K."""

import json
import logging
import os
from typing import List

import click

from pybel.struct import count_functions, count_namespaces
from snp2k import ChromosomeMapper, get_graph

__all__ = [
    'main',
]


@click.command()
@click.option('-f', '--force', is_flag=True, help='Force reconstruction of all graphs')
@click.option('-x', '--force-global', is_flag=True, help='Only force reconstruction of global graph')
@click.option('-p', '--packages', help='Bio2BEL package prefixes')
@click.option('-d', '--directory', default=os.getcwd(), type=click.Path(dir_okay=True, exists=True, file_okay=False))
@click.option('-v', '--debug', is_flag=True, help='Set logging level to INFO')
def main(force: bool, force_global: bool, packages: List[str], directory: str, debug: bool):
    """Summarize the graph."""
    if debug:
        logging.basicConfig(level=logging.INFO)

    graph = get_graph(
        force=force,
        force_global=force_global,
        names=packages,
    )

    click.echo(graph.summary_str())
    click.echo(str(count_functions(graph)))
    click.echo(str(count_namespaces(graph)))

    mapper = ChromosomeMapper()
    mapper.update_chromosome_map(graph)

    with open(os.path.join(directory, 'chromosomes_single.json'), 'w') as file:
        json.dump(mapper.get_chromosome_count(), file, indent=2)

    cross_chromosome_count = mapper.get_cross_chromosome_count()
    with open(os.path.join(directory, 'chromosomes_paired.json'), 'w') as file:
        json.dump(
            [
                {
                    'source': s,
                    'target': t,
                    'count': count,
                }
                for (s, t), count in cross_chromosome_count.items()
            ],
            file,
            indent=2,
        )

    triple_df = mapper.get_triple_df()
    triple_df.to_csv(os.path.join(directory, 'triple_map.tsv'), sep='\t', index=False)


if __name__ == '__main__':
    main()
