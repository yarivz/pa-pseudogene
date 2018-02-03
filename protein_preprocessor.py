import gzip
import logging
import multiprocessing
import os
import shutil

from Bio import SeqIO

from logging_config import worker_configurer, log_queue
from constants import DATA_DIR, STRAINS_DIR, NUMBER_OF_PROCESSES, FASTA_FILE_TYPE, PROTEIN_FILE_PATTERN, \
    CDS_FROM_GENOMIC_PATTERN, STRAIN_INDEX_FILE, COMBINED_STRAIN_PROTEINS_PREFIX, WORKER_PROTEIN_FILE_PREFIX, \
    COMBINED_PROTEINS_FILE_PATH


def create_all_strains_file_with_indices():
    """
    Preprocess all strains proteins in parallel and combine all worker output files into single fasta file for clustering
    """
    logger = logging.getLogger()
    logger.info("Indexing proteins by their strain index & protein index in strain gene")
    for file in os.listdir(DATA_DIR):
        if COMBINED_STRAIN_PROTEINS_PREFIX in file:
            os.remove(os.path.join(DATA_DIR, file))
    job_queue = multiprocessing.Queue()
    prepare_preprocessing_jobs(job_queue)
    workers = [multiprocessing.Process(target=preprocess_strain_proteins, args=(i, job_queue, worker_configurer))
               for i in range(NUMBER_OF_PROCESSES)]
    for w in workers:
        w.start()
    job_queue.put(None)
    for w in workers:
        w.join()
    worker_combined_protein_files = [f for f in os.listdir(DATA_DIR) if WORKER_PROTEIN_FILE_PREFIX in f]
    with open(COMBINED_PROTEINS_FILE_PATH, 'w') as dstd:
        for worker_file in worker_combined_protein_files:
            with open(os.path.join(DATA_DIR, worker_file), 'r') as srcd:
                shutil.copyfileobj(srcd, dstd)


def prepare_preprocessing_jobs(job_queue):
    """Put all downloaded strain dirs in job queue for workers"""
    downloaded_strains = os.listdir(STRAINS_DIR)
    for strain_dir in downloaded_strains:
        job_queue.put(strain_dir)


def preprocess_strain_proteins(worker_id, job_queue, configurer):
    """
    Preprocess all proteins for downloaded strains by indexing according to strain index and protein index within
    strain genome, using the cds_from_genome file
    """
    configurer(log_queue)
    logger = logging.getLogger(__name__ + "_worker_" + str(worker_id))
    worker_combined_proteins_file_path = os.path.join(DATA_DIR, WORKER_PROTEIN_FILE_PREFIX + str(worker_id))
    while True:
        strain_dir = job_queue.get()
        if strain_dir is None:
            job_queue.put(None)
            break
        strain_dir_files = os.listdir(os.path.join(STRAINS_DIR, strain_dir))
        protein_file_name = [f for f in strain_dir_files if PROTEIN_FILE_PATTERN in f][0]
        cds_file_name = [f for f in strain_dir_files if CDS_FROM_GENOMIC_PATTERN in f][0]
        if not protein_file_name or not cds_file_name:
            logger.warning(
                "Could not find protein file or cds_from_genomic file for strain %s, skipping" % strain_dir)
            continue
        protein_file = cds_file = strain_index_file = worker_combined_proteins_file = None
        try:
            strain_index_file = open(os.path.join(STRAINS_DIR, strain_dir, STRAIN_INDEX_FILE))
            strain_index = '[' + strain_index_file.readline() + ']'
            if protein_file_name.endswith('gz'):
                protein_file = gzip.open(os.path.join(STRAINS_DIR, strain_dir, protein_file_name), 'rt')
            else:
                protein_file = open(os.path.join(STRAINS_DIR, strain_dir, protein_file_name))
            if cds_file_name.endswith('gz'):
                cds_file = gzip.open(os.path.join(STRAINS_DIR, strain_dir, cds_file_name), 'rt')
            else:
                cds_file = open(os.path.join(STRAINS_DIR, strain_dir, cds_file_name))

            strain_cds_seq_dict = SeqIO.to_dict(SeqIO.parse(cds_file, FASTA_FILE_TYPE))
            strain_cds_seq_headers = list(strain_cds_seq_dict.keys())
            strain_protein_seq_iter = SeqIO.parse(protein_file, FASTA_FILE_TYPE)
            for strain_protein_seq in strain_protein_seq_iter:
                protein_id = strain_protein_seq.id
                cds_protein_header = [h for h in strain_cds_seq_headers if protein_id in h][0]
                protein_index_in_gene = '[' + cds_protein_header[cds_protein_header.rfind('_') + 1:] + ']'
                strain_protein_seq.description = strain_index + protein_index_in_gene + strain_protein_seq.description
                strain_protein_seq.id = ""
                worker_combined_proteins_file = open(worker_combined_proteins_file_path, 'a+')
                SeqIO.write(strain_protein_seq, worker_combined_proteins_file, FASTA_FILE_TYPE)
            logger.info(
                "Strain %s proteins were indexed and written to file" % strain_dir[strain_dir.rfind(']') + 1:])
        finally:
            if protein_file is not None:
                protein_file.close()
            if cds_file is not None:
                cds_file.close()
            if worker_combined_proteins_file is not None:
                worker_combined_proteins_file.close()
            if strain_index_file is not None:
                strain_index_file.close()
