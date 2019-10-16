# -*- coding: utf-8 -*-

"""Acquisition of relevant Bio2BEL resources for SNP2K."""

import importlib
import logging
import os
from types import ModuleType
from typing import Any, List, Mapping, Optional, Tuple, Type, Union

import networkx as nx

import bio2bel_hgnc
import pybel
from bio2bel.manager.bel_manager import BELManagerMixin
from pybel import BELGraph, from_pickle, to_pickle
from pybel.struct import count_functions, count_namespaces, enrich_protein_and_rna_origins
from snp2k.constants import RESOURCES

__all__ = [
    'get_graph',
]

logger = logging.getLogger(__name__)

BELModule = Union[str, ModuleType, BELManagerMixin, Type[BELManagerMixin]]

NamesList = List[Tuple[BELModule, Mapping[str, Any]]]

DEFAULT_NAMES: NamesList = [
    ('neurommsig', {}),
    ('kegg', {}),
    ('reactome', {}),
    ('adeptus', {}),
    ('hippie', dict(namespace='hgnc')),
    ('gwascatalog', {}),
    ('wikipathways', {}),
    # ('phewascatalog', {}),
]

try:
    import compath_resources
except ImportError:
    pass
else:
    DEFAULT_NAMES.insert(0, (compath_resources.Manager, {}))

CACHE_NAME = 'snp2k.bel.pickle'


def get_graph(
    force: bool = False,
    force_global: bool = False,
    names: Optional[NamesList] = None,
    resources_directory: Optional[str] = None,
) -> BELGraph:
    """Get all resources in a combine BELGraph.

    :param force: Should cached files be overwritten?
    :param force_global: Should the global cache file be overwritten?
    :param names: The name of the bio2bel packages to use and arguments
    :param resources_directory: A non-default place to store the resources
    """
    pickle_path = os.path.join(resources_directory or RESOURCES, CACHE_NAME)
    if not force_global and os.path.exists(pickle_path):
        logger.info(f'Getting cached full graph')
        return from_pickle(pickle_path)

    if names is None:
        names = DEFAULT_NAMES

    logger.info('Generating graphs')
    graphs = []
    for name, to_bel_kwargs in names:
        _graph = get_graph_by_manager(name, force=force, to_bel_kwargs=to_bel_kwargs)
        logger.info(_graph.summary_str())
        graphs.append(_graph)

    logger.info('Merging graphs')
    graph = pybel.union(graphs)
    graph.name = f'Graph from: {", ".join(graph.name for graph in graphs)}'
    graph.version = '0.0.1'
    logger.info('Finished merging graphs')

    logger.info('Preparing HGNC mappings')
    hgnc_manager = bio2bel_hgnc.Manager()
    hgnc_symbol_to_id = hgnc_manager.build_hgnc_symbol_id_mapping()
    entrez_id_to_hgnc_symbol = hgnc_manager.build_entrez_id_to_hgnc_symbol_mapping()

    logger.info('Generating namespace mapping for nodes')
    mapping = {}
    for node in graph:
        namespace = node.get('namespace')
        if namespace is None:
            continue
        elif namespace.lower() in {'ncbigene', 'egid'} and node.identifier in entrez_id_to_hgnc_symbol:
            name = entrez_id_to_hgnc_symbol[node.identifier]
            identifier = hgnc_symbol_to_id[name]
            mapping[node] = node.__class__(
                namespace='hgnc',
                name=name,
                identifier=identifier,
            )

    logger.info('Relabeling nodes')
    nx.relabel_nodes(graph, mapping, copy=False)

    logger.info('Enriching central dogma')
    enrich_protein_and_rna_origins(graph)

    logger.info('Exporting snp2k pickle')
    to_pickle(graph, pickle_path)
    return graph


def get_graph_by_manager(
    module: Union[str, ModuleType, BELManagerMixin, Type[BELManagerMixin]],
    force: bool = False,
    to_bel_kwargs: Optional[Mapping[str, Any]] = None,
) -> BELGraph:
    """Get a graph for a manager."""
    if isinstance(module, str):  # get the cache or import that module
        _pickle_path = os.path.join(RESOURCES, f'{module}.bel.pickle')
        if os.path.exists(_pickle_path) and not force:
            logger.info(f'Getting {module} from pickle at {_pickle_path}')
            return from_pickle(_pickle_path)

        module_name = f'bio2bel_{module}'
        _module = importlib.import_module(module_name)
        manager = _module.Manager()
    elif isinstance(module, BELManagerMixin):
        manager = module
    elif isinstance(module, ModuleType):
        manager = module.Manager()
    elif isinstance(module, type):
        if not issubclass(module, BELManagerMixin):
            raise TypeError(f'{module} is not a subclass of BELManagerMixin')
        manager = module()
    else:
        raise TypeError(f'{module} has invalid type: {type(module)}')

    pickle_path = os.path.join(RESOURCES, f'{manager.module_name}.bel.pickle')
    if os.path.exists(pickle_path) and not force:
        logger.info(f'Getting {manager.module_name} from pickle at {pickle_path}')
        return from_pickle(pickle_path)

    if not manager.is_populated():
        logger.info(f'Populating manager for {manager.module_name}')
        manager.populate()

    graph = manager.to_bel(**(to_bel_kwargs or {}))
    logger.info(graph.summary_str())
    logger.info(str(count_namespaces(graph)))
    logger.info(str(count_functions(graph)))

    logger.info(f'Writing pickle for {pickle_path}')
    to_pickle(graph, pickle_path)
    return graph
