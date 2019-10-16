# -*- coding: utf-8 -*-

"""The chromsome mapper identifies which edges in a BEL graph are incident to which chromosomes."""

import itertools as itt
import logging
import re
from collections import defaultdict
from typing import Any, Dict, List, Iterable, Mapping, Optional, Tuple

import pandas as pd
from tqdm import tqdm

import bio2bel_hgnc
import pybel
import pybel.dsl
import pybel.struct
from pybel.canonicalize import edge_to_bel
from bio2bel_entrez.parser import get_human_refseq_slim_df
from bio2bel_hgnc.models import HumanGene
from .constants import CHROMOSOMES

__all__ = [
    'ChromosomeMapper',
]

CHROMOSOME_SPLIT_RE = re.compile(r'p|q')

Edge = Tuple[pybel.BaseEntity, pybel.BaseEntity, str, Mapping[str, Any]]

logger = logging.getLogger(__name__)


class ChromosomeMapper:
    """The chromosome mapper identifies which edges in a BEL graph are incident to which chromosomes."""

    chromosome_to_edge_keys: Dict[str, List[Edge]]
    cross_chromosome_to_edge_keys: Dict[Tuple[str, str], List[Edge]]

    def __init__(self, hgnc_manager: Optional[bio2bel_hgnc.Manager] = None) -> None:
        if hgnc_manager is None:
            logger.info('getting Bio2BEL HGNC manager')
            hgnc_manager = bio2bel_hgnc.Manager()

        if not hgnc_manager.is_populated():
            logger.info('populating HGNC')
            hgnc_manager.populate()

        logger.info('generating hgnc symbol to chromosome mapping')
        self.hgnc_symbol_to_chromosome = {
            str(symbol): CHROMOSOME_SPLIT_RE.split(location)[0]
            for symbol, location in hgnc_manager.session.query(HumanGene.symbol, HumanGene.location)
            if location is not None

        }

        logger.info('generating hgnc id to chromosome mapping')
        self.hgnc_id_to_chromosome = {
            str(hgnc_id): CHROMOSOME_SPLIT_RE.split(location)[0]
            for hgnc_id, location in hgnc_manager.session.query(HumanGene.identifier, HumanGene.location)
            if location is not None
        }

        self.hgnc_id_to_positions = {}
        self.hgnc_symbol_to_positions = {}

        logger.info('generating hgnc symbol to chromosome mapping')
        self.entrez_id_to_hgnc_id = hgnc_manager.build_entrez_id_to_hgnc_id_mapping()

        logger.info('getting human refseq data')
        self.human_refseq_df = get_human_refseq_slim_df()

        logger.info('generating maps with refseq data')
        for entrez_id, symbol, start, end in self.human_refseq_df.values:
            hgnc_id = self.entrez_id_to_hgnc_id.get(str(entrez_id))
            if hgnc_id is None:
                logger.debug(f'Could not find ncbigene:{entrez_id} in HGNC. May be withdrawn')
                continue
            self.hgnc_id_to_positions[hgnc_id] = start, end
            self.hgnc_symbol_to_positions[symbol] = start, end

        # These will get populated as graphs are added with update_chromosome_map()
        self.chromosome_to_edge_keys = defaultdict(list)
        self.cross_chromosome_to_edge_keys = defaultdict(list)

    def get_node_chromosomes(self, node: pybel.dsl.BaseEntity) -> Iterable[str]:
        if isinstance(node, pybel.dsl.BaseAbundance) and node.namespace == 'hgnc':
            if node.identifier:
                chromosome = self.hgnc_id_to_chromosome.get(node.identifier)
            else:
                chromosome = self.hgnc_symbol_to_chromosome.get(node.name)

                # TODO remove this after fixing name/identifiers in nodes
            if chromosome is None:
                chromosome = self.hgnc_id_to_chromosome.get(node.name)

            yield chromosome

        elif isinstance(node, pybel.dsl.ComplexAbundance):
            for member in node.members:
                yield from self.get_node_chromosomes(member)

    def update_chromosome_map(self, graph: pybel.BELGraph, use_tqdm: bool = True):
        logger.info(f'updating chromosome map with {graph}')

        it = graph.edges(keys=True, data=True)
        if use_tqdm:
            it = tqdm(it, desc='Mapping edges')
        for u, v, key, d in it:
            all_relevant_chromosomes = list(itt.chain(
                self.get_node_chromosomes(u),
                self.get_node_chromosomes(v),
            ))

            for chromosome in all_relevant_chromosomes:
                self.chromosome_to_edge_keys[chromosome].append((u, v, key, d))

            for c1, c2 in itt.combinations(all_relevant_chromosomes, 2):
                self.cross_chromosome_to_edge_keys[c1, c2].append((u, v, key, d))

    def get_chromosome_count(self) -> Mapping[str, int]:
        return {
            chromosome: len(edges)
            for chromosome, edges in self.chromosome_to_edge_keys.items()
            if chromosome in CHROMOSOMES
        }

    def get_cross_chromosome_count(self) -> Mapping[Tuple[str, str], int]:
        return {
            (c1, c2): len(edges)
            for (c1, c2), edges in self.cross_chromosome_to_edge_keys.items()
            if c1 in CHROMOSOMES and c2 in CHROMOSOMES
        }

    def get_cross_chromosome_count_df(self) -> pd.DataFrame:
        cross_chromosome_count = self.get_cross_chromosome_count()

        return pd.DataFrame(
            [
                [
                    cross_chromosome_count.get((c1, c2), 0)
                    for c2 in CHROMOSOMES
                ]
                for c1 in CHROMOSOMES
            ],
            index=CHROMOSOMES,
            columns=CHROMOSOMES,
        )

    def get_triple_df(self) -> pd.DataFrame:
        """Get a dataframe of all BEL triples with their objects and positions.

        The resulting dataframe has the following columns:

        - BEL statement
        - Subject HGNC ID
        - Subject HGNC gene symbol
        - Subject chromosome
        - Subject locus
        - Object HGNC ID
        - Object HGNC gene symbol
        - Object chromosome
        - Object locus

        If a BEL statement contains more than two HGNC nodes, then it will give
        one row for each possible combination.
        """
        rows = self._get_triples_iter()
        return pd.DataFrame(rows, columns=[
            'BEL',
            'Subject HGNC ID',
            'Subject HGNC Symbol',
            'Subject chromosome',
            'Subject locus',
            'Object HGNC ID',
            'Object HGNC Symbol',
            'Object chromosome',
            'Object locus',
        ])

    def _get_triples_iter(self):
        """Get the information from a triple:
        - the BEL statement, and
        - for each node in the triple: identifier, name, chromosome it is on, and locus on that chromosome
        """
        for (u_chromosome, v_chromosome), edges in self.cross_chromosome_to_edge_keys.items():
            for u, v, k, d in edges:
                bel = edge_to_bel(u, v, d)
                it = itt.product(self._iter_id_name_loci(u), self._iter_id_name_loci(v))
                for (u_identifier, u_name, u_locus), (v_identifier, v_name, v_locus) in it:
                    yield (
                        bel,
                        u_identifier, u_name, u_chromosome, u_locus,
                        v_identifier, v_name, v_chromosome, v_locus,
                    )

    def _iter_id_name_loci(self, node: pybel.BaseEntity):
        """Iterate over all possible pairs from this node."""
        if isinstance(node, pybel.dsl.CentralDogma):
            # print(node.namespace)
            if node.namespace.lower() == "hgnc":
                identifier = node.identifier
                name = node.name
                locus = self.hgnc_id_to_positions[identifier]
                # try:
                #     locus = self.hgnc_id_to_positions[identifier]
                # except KeyError:
                #     print(node)
                yield identifier, name, locus

        elif isinstance(node, pybel.dsl.ListAbundance):
            for member in node.members:
                print(member.namespace)
                yield from self._iter_id_name_loci(member)

        else:
            logger.debug(f"Unhandled node: {node}")

