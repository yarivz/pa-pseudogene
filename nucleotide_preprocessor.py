import gzip
import logging
import multiprocessing
import os
import shutil

from Bio import SeqIO

from logging_config import worker_configurer
from constants import DATA_DIR, STRAINS_DIR, NUMBER_OF_PROCESSES, FASTA_FILE_TYPE, PROTEIN_FILE_PATTERN, \
    CDS_FROM_GENOMIC_PATTERN, STRAIN_INDEX_FILE, COMBINED_STRAIN_PROTEINS_PREFIX, WORKER_PROTEIN_FILE_PREFIX, \
    COMBINED_PROTEINS_FILE_PATH


