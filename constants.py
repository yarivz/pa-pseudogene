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
CLUSTER_STRAIN_PATTERN = re.compile("[0-9a-z,> \t]+\[(\d+)\]\[(\d+)\]")
CLUSTER_PSEUDOGENE_PATTERN = re.compile(CLUSTER_STRAIN_PATTERN.pattern + "\[p")
CLUSTER_1ST_STAGE_REPRESENTATIVE_PATTERN = re.compile(CLUSTER_STRAIN_PATTERN.pattern + "\[cluster_(\d+)\]")
CLUSTER_2ND_STAGE_SEQ_LEN_PATTERN = re.compile("(\d+)nt,")

COMBINED_STRAIN_PROTEINS_PREFIX = "combined_strain_proteins"
WORKER_PROTEIN_FILE_PREFIX = COMBINED_STRAIN_PROTEINS_PREFIX + "_worker_"
COMBINED_PROTEINS_FILE_PATH = os.path.join(DATA_DIR, COMBINED_STRAIN_PROTEINS_PREFIX + "_all.fasta")
COMBINED_STRAIN_CDS_PREFIX = "combined_strain_reps_pseudogenes_cds"
WORKER_CDS_FILE_PREFIX = COMBINED_STRAIN_CDS_PREFIX + "_worker_"
COMBINED_CDS_FILE_PATH = os.path.join(DATA_DIR, COMBINED_STRAIN_CDS_PREFIX + "_all.fasta")

CD_HIT_CLUSTER_REPS_OUTPUT_FILE = os.path.join(DATA_DIR, 'protein_clusters.txt')
CD_HIT_CLUSTERS_OUTPUT_FILE = CD_HIT_CLUSTER_REPS_OUTPUT_FILE + ".clstr"
CD_HIT_EST_CLUSTER_REPS_OUTPUT_FILE = os.path.join(DATA_DIR, 'cds_clusters.txt')
CD_HIT_EST_CLUSTERS_OUTPUT_FILE = CD_HIT_EST_CLUSTER_REPS_OUTPUT_FILE + ".clstr"

GENOMIC_STATS_PKL = os.path.join(DATA_DIR, "genomic_stats.pkl")
PROTEIN_STATS_PKL = os.path.join(DATA_DIR, "protein_stats.pkl")
CONTIGS_PKL = os.path.join(DATA_DIR, "contigs.pkl")
PSEUDOGENES_PKL = os.path.join(DATA_DIR, "pseudogenes.pkl")
TOTAL_CLUSTERS_PKL = os.path.join(DATA_DIR, "total_clusters.pkl")
CORE_CLUSTERS_PKL = os.path.join(DATA_DIR, "core_clusters.pkl")
SINGLETON_CLUSTERS_PKL = os.path.join(DATA_DIR, "singleton_clusters.pkl")
FIRST_STAGE_STATS_PKL = os.path.join(DATA_DIR, "1st_stage_stats.pkl")
SECOND_STAGE_STRAIN_STATS_PKL = os.path.join(DATA_DIR, "2nd_stage_strain_stats.pkl")
SECOND_STAGE_CLUSTER_STATS_PKL = os.path.join(DATA_DIR, "2nd_stage_cluster_stats.pkl")
SECOND_STAGE_AGGREGATED_CLUSTER_STATS_PKL = os.path.join(DATA_DIR, "2nd_stage_aggregated_cluster_stats.pkl")

FIRST_STAGE_STATS_CSV = os.path.join(DATA_DIR, "1st_stage_stats.csv")
SECOND_STAGE_STATS_CSV = os.path.join(DATA_DIR, "2nd_stage_stats.csv")
