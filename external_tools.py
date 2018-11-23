import logging
import multiprocessing
import os
from subprocess import run

from constants import CD_HIT_CLUSTER_REPS_OUTPUT_FILE, CLUSTERS_NT_SEQS_DIR, CLUSTERS_ALIGNMENTS_DIR, STRAINS_DIR, \
    NUMBER_OF_PROCESSES
from logging_config import worker_configurer


def perform_clustering_on_proteins(aggregated_proteins_file_path):
    """Run the CD-HIT program to perform clustering on the strains"""
    logger = logging.getLogger()
    logger.info("Running CD-HIT on combined proteins file to create clustering")
    cd_hit_args = " ".join(["cd-hit", "-i", aggregated_proteins_file_path, "-o", CD_HIT_CLUSTER_REPS_OUTPUT_FILE, "-c 0.70",
                   "-n 5", "-M 16000", "-g 1", "-p 1"])
    cd_hit_return_code = run(cd_hit_args, shell=True).returncode
    logger.info("Finished running CD-HIT with return code %d" % cd_hit_return_code)
    return cd_hit_return_code


def perform_clustering_on_cds(input_file, output_file):
    """Run the CD-HIT-EST program to perform clustering on the strains representatives and pseudogenes"""
    logger = logging.getLogger()
    logger.info("Running CD-HIT-EST on combined representative and pseudogene cds file to create clustering")
    cd_hit_est_args = " ".join(["cd-hit-est", "-i", input_file, "-o", output_file, "-c 0.8",
                   "-n 5", "-M 16000", "-g 1", "-p 1", "-d 30"])
    cd_hit_est_return_code = run(cd_hit_est_args, shell=True).returncode
    logger.info("Finished running CD-HIT with return code %d" % cd_hit_est_return_code)
    return cd_hit_est_return_code


def perform_alignment_on_core_clusters(log_queue):
    """Run MAFFT & Gblocks tools on fasta files of protein nucleotide seqs for each core cluster"""
    logger = logging.getLogger(__name__)
    logger.info("Running MAFFT & Gblocks on core clusters for alignment")
    if not os.path.exists(CLUSTERS_NT_SEQS_DIR):
        logger.error("No clusters dir found, exiting")
        exit(1)
    if not os.path.exists(CLUSTERS_ALIGNMENTS_DIR):
        os.makedirs(CLUSTERS_ALIGNMENTS_DIR)

    job_queue = multiprocessing.Queue()
    prepare_alignment_jobs(job_queue)
    workers = [
        multiprocessing.Process(target=perform_alignment_and_pruning, args=(i, job_queue, worker_configurer, log_queue))
        for i in range(NUMBER_OF_PROCESSES)]
    for w in workers:
        w.start()
    job_queue.put(None)
    for w in workers:
        w.join()
    logger.info("Finished running MAFFT for all clusters")


def prepare_alignment_jobs(job_queue):
    """Put all downloaded strain dirs in job queue for workers"""
    core_clusters = os.listdir(CLUSTERS_NT_SEQS_DIR)
    for cluster_file in core_clusters:
        job_queue.put(cluster_file)


def perform_alignment_and_pruning(worker_id, job_queue, configurer, log_queue):
    """
    Perform MAFFT alignment and Gblocks pruning for a core cluster fasta file
    """
    configurer(log_queue)
    logger = logging.getLogger(__name__ + "_worker_" + str(worker_id))
    while True:
        cluster_file = job_queue.get()
        if cluster_file is None:
            job_queue.put(None)
            break
        logger.info("Running MAFFT for %s" % cluster_file)
        cluster_alignment_filename = cluster_file + "_alignment"
        cluster_alignment_file = open(os.path.join(CLUSTERS_ALIGNMENTS_DIR, cluster_alignment_filename), 'w')
        mafft_args = " ".join(["mafft", "--auto", os.path.join(CLUSTERS_NT_SEQS_DIR, cluster_file)])
        mafft_return_code = run(mafft_args, shell=True, stdout=cluster_alignment_file).returncode
        cluster_alignment_file.close()
        logger.info("Finished running MAFFT for %s with return code %d" % (cluster_file, mafft_return_code))
        logger.info("Running GBlocks for %s" % cluster_file)
        gblocks_args = " ".join(["Gblocks", os.path.join(CLUSTERS_ALIGNMENTS_DIR, cluster_alignment_filename), "-t=d", "-b5=a"])
        gblocks_return_code = run(gblocks_args, shell=True).returncode
        logger.info(
            "Finished running Gblocks for alignment %s with return code %d" % (cluster_alignment_filename, gblocks_return_code))


# def perform_pruning_on_alignments():
#     """Run Gblocks tool on MAFFT alignments of protein nucleotide seqs for each core cluster"""
#     logger = logging.getLogger()
#     logger.info("Running Gblocks on core clusters alignments")
#     if not os.path.exists(CLUSTERS_ALIGNMENTS_DIR):
#         logger.error("No clusters alignments dir found, exiting")
#         exit(1)
#
#     for alignment_file in os.listdir(CLUSTERS_ALIGNMENTS_DIR):
#         gblocks_args = " ".join(["Gblocks", os.path.join(CLUSTERS_ALIGNMENTS_DIR, alignment_file), "-t=d", "-b5=a"])
#         gblocks_return_code = run(gblocks_args, shell=True).returncode
#         logger.info("Finished running Gblocks for alignment %s with return code %d" % (alignment_file, gblocks_return_code))
#     logger.info("Finished running Gblocks for all alignments")

   # for cluster_file in os.listdir(CLUSTERS_NT_SEQS_DIR):
   #      cluster_alignment_filename = cluster_file + "_alignment"
   #      cluster_alignment_file = open(os.path.join(CLUSTERS_ALIGNMENTS_DIR, cluster_alignment_filename), 'w')
   #      mafft_args = " ".join(["mafft", "--auto", os.path.join(CLUSTERS_NT_SEQS_DIR, cluster_file)])
   #      mafft_return_code = run(mafft_args, shell=True, stdout=cluster_alignment_file, stderr=cluster_alignment_file).returncode
   #      logger.info("Finished running MAFFT for %s with return code %d" % (cluster_file, mafft_return_code))
   #      cluster_alignment_file.close()