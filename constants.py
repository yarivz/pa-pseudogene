import os

DATA_DIR = os.getcwd() + os.sep + "data"
STRAINS_DIR = DATA_DIR + os.sep + "strains"
NUMBER_OF_PROCESSES = os.cpu_count()
FASTA_FILE_TYPE = "fasta"
PROTEIN_FILE_PATTERN = "protein.faa"
CDS_FROM_GENOMIC_PATTERN = "cds_from_genomic.fna"
STRAIN_INDEX_FILE = "strain_index"
FEATURE_TABLE_PATTERN = "feature_table.txt"