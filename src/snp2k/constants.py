# -*- coding: utf-8 -*-

"""Constants for SNP2K."""

import os

__all__ = [
    'HERE',
    'RESOURCES',
    'CHROMOSOMES',
    'CHROMOSOME_TO_INDEX',
]

HERE = os.path.abspath(os.path.dirname(__file__))
RESOURCES = os.path.join(HERE, 'resources')

CHROMOSOMES = [str(i) for i in range(1, 23)] + ['X', 'Y', 'mitochondria']
CHROMOSOME_TO_INDEX = {
    chromosome: i
    for i, chromosome in enumerate(CHROMOSOMES)
}
