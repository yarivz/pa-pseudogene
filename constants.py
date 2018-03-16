import os
import re

DATA_DIR = os.getcwd() + os.sep + "data"
STRAINS_DIR = DATA_DIR + os.sep + "strains"
NUMBER_OF_PROCESSES = os.cpu_count()
FASTA_FILE_TYPE = "fasta"
PROTEIN_FILE_PATTERN = "protein.faa"
CDS_FROM_GENOMIC_PATTERN = "cds_from_genomic.fna"
GENOMIC_PATTERN = "genomic.fna"
STRAIN_INDEX_FILE = "strain_index"
FEATURE_TABLE_PATTERN = "feature_table.txt"

COMBINED_STRAIN_PROTEINS_PREFIX = "combined_strain_proteins"
WORKER_PROTEIN_FILE_PREFIX = COMBINED_STRAIN_PROTEINS_PREFIX + "_worker_"
COMBINED_PROTEINS_FILE_PATH = os.path.join(DATA_DIR, COMBINED_STRAIN_PROTEINS_PREFIX + "_all.fasta")
COMBINED_STRAIN_CDS_PREFIX = "combined_strain_reps_pseudogenes_cds"
WORKER_CDS_FILE_PREFIX = COMBINED_STRAIN_CDS_PREFIX + "_worker_"
COMBINED_CDS_FILE_PATH = os.path.join(DATA_DIR, COMBINED_STRAIN_CDS_PREFIX + "_all.fasta")

CD_HIT_CLUSTER_REPS_OUTPUT_FILE = os.path.join(DATA_DIR, 'protein_clusters.txt')
CD_HIT_CLUSTERS_OUTPUT_FILE = CD_HIT_CLUSTER_REPS_OUTPUT_FILE + ".clstr"
GENOMIC_STATS_PKL = os.path.join(DATA_DIR, "genomic_stats.pkl")
PROTEIN_STATS_PKL = os.path.join(DATA_DIR, "protein_stats.pkl")
CONTIGS_PKL = os.path.join(DATA_DIR, "contigs.pkl")
PSEUDOGENES_PKL = os.path.join(DATA_DIR, "pseudogenes.pkl")
TOTAL_CLUSTERS_PKL = os.path.join(DATA_DIR, "total_clusters.pkl")
CORE_CLUSTERS_PKL = os.path.join(DATA_DIR, "core_clusters.pkl")
SINGLETON_CLUSTERS_PKL = os.path.join(DATA_DIR, "singleton_clusters.pkl")
CLUSTER_STRAIN_PATTERN = re.compile("[0-9a,> \t]+\[(\d+)\]\[(\d+)\]")

