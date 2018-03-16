import gzip
import logging
import multiprocessing
import os
import shutil

from Bio import SeqIO

from logging_config import worker_configurer
from constants import DATA_DIR, STRAINS_DIR, NUMBER_OF_PROCESSES, FASTA_FILE_TYPE, CDS_FROM_GENOMIC_PATTERN, \
    STRAIN_INDEX_FILE, CLUSTER_STRAIN_PATTERN, COMBINED_STRAIN_CDS_PREFIX, WORKER_CDS_FILE_PREFIX, \
    COMBINED_CDS_FILE_PATH, CD_HIT_CLUSTERS_OUTPUT_FILE


class Representative:
    def __init__(self, position, cluster_index):
        self.position = position
        self.cluster_index = cluster_index


class StrainData:
    def __init__(self, index, reps, strain_dir):
        self.index = index
        self.representatives = reps
        self.dir = strain_dir


def create_representatives_and_pseudogenes_file(log_queue):
    """
    Preprocess all strains representative and pseudogene nucleutide sequences in parallel and combine all worker output
    files into single fasta file for clustering
    """
    logger = logging.getLogger(__name__)
    logger.info("Preprocessing cds sequences for cluster representative proteins and pseudogenes")
    representatives = get_clusters_representatives(CD_HIT_CLUSTERS_OUTPUT_FILE)
    logger.info("representatives: %s" % representatives)
    for file in os.listdir(DATA_DIR):
        if COMBINED_STRAIN_CDS_PREFIX in file:
            os.remove(os.path.join(DATA_DIR, file))
    job_queue = multiprocessing.Queue()
    prepare_preprocessing_jobs(job_queue, representatives)
    workers = [multiprocessing.Process(target=preprocess_strain_cds, args=(i, job_queue, worker_configurer, log_queue))
               for i in range(NUMBER_OF_PROCESSES)]
    for w in workers:
        w.start()
    job_queue.put(None)
    for w in workers:
        w.join()
    worker_combined_cds_files = [f for f in os.listdir(DATA_DIR) if WORKER_CDS_FILE_PREFIX in f]
    with open(COMBINED_CDS_FILE_PATH, 'w') as dstd:
        for worker_file in worker_combined_cds_files:
            with open(os.path.join(DATA_DIR, worker_file), 'r') as srcd:
                shutil.copyfileobj(srcd, dstd)


def get_clusters_representatives(clusters_file):
    representatives = {}
    with open(clusters_file, 'r') as clusters_db:
        print("getting reps")
        for line in clusters_db:
            if line.startswith(">Cluster"):
                cluster_index = int(line.split()[-1])
                print("got cluster " + str(cluster_index))
            else:
                if line.endswith("*"):
                    match = CLUSTER_STRAIN_PATTERN.match(line)
                    strain_index = int(match.group(1))
                    print("adding rep " + match.group(0))
                    strain_reps = representatives[strain_index] if strain_index in representatives.keys() else []
                    strain_reps.append(Representative(int(match.group(2)), cluster_index))
                    representatives[strain_index] = strain_reps
    return representatives


def prepare_preprocessing_jobs(job_queue, representatives):
    """Put all downloaded strain data in job queue for workers"""
    downloaded_strains = os.listdir(STRAINS_DIR)
    for strain_dir in downloaded_strains:
        strain_index = get_strain_index(strain_dir)
        strain_reps = representatives[strain_index] if strain_index in representatives.keys() else []
        strain_data = StrainData(strain_index, strain_reps, strain_dir)
        job_queue.put(strain_data)


def get_strain_index(strain_dir):
    with open(os.path.join(STRAINS_DIR, strain_dir, STRAIN_INDEX_FILE)) as f:
        return int(f.readline())


def preprocess_strain_cds(worker_id, job_queue, configurer, log_queue):
    """
    Preprocess all cds for downloaded strains by indexing according to strain index and protein index within
    strain genome, using the cds_from_genome file
    """
    configurer(log_queue)
    logger = logging.getLogger(__name__ + "_worker_" + str(worker_id))
    worker_combined_cds_file_path = os.path.join(DATA_DIR, WORKER_CDS_FILE_PREFIX + str(worker_id))
    while True:
        strain_data = job_queue.get()
        if strain_data is None:
            job_queue.put(None)
            break
        strain_dir = strain_data.dir
        strain_dir_files = os.listdir(os.path.join(STRAINS_DIR, strain_dir))
        cds_file_name = [f for f in strain_dir_files if CDS_FROM_GENOMIC_PATTERN in f][0]
        cds_file = worker_combined_cds_file = None
        try:
            if cds_file_name.endswith('gz'):
                cds_file = gzip.open(os.path.join(STRAINS_DIR, strain_dir, cds_file_name), 'rt')
            else:
                cds_file = open(os.path.join(STRAINS_DIR, strain_dir, cds_file_name))
            strain_cds_seq_iter = SeqIO.parse(cds_file, FASTA_FILE_TYPE)
            for strain_cds_seq in strain_cds_seq_iter:
                seq_position_in_genome = int(strain_cds_seq.id[strain_cds_seq.id.rfind("_") + 1:])
                rep = [r for r in strain_data.representatives if seq_position_in_genome == r.position]
                pseudo = "pseudo=true" in strain_cds_seq.description
                if rep or pseudo:
                    strain_cds_seq.description = "[" + str(strain_data.index) + "]" + "[" + str(seq_position_in_genome) + "]"\
                                                 + "{info}".format(info="[cluster_" + str(rep[0].cluster_index) + "]" if rep else "[pseudo]")\
                                                 + strain_cds_seq.description
                    worker_combined_cds_file = open(worker_combined_cds_file_path, 'a+')
                    SeqIO.write(strain_cds_seq, worker_combined_cds_file, FASTA_FILE_TYPE)
            logger.info(
                "Strain %s reps and pseudogenes were indexed and written to file" % strain_dir[strain_dir.rfind(']') + 1:])
        finally:
            if cds_file is not None:
                cds_file.close()
            if worker_combined_cds_file is not None:
                worker_combined_cds_file.close()
