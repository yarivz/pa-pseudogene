import logging
import os
from subprocess import run

from constants import CD_HIT_CLUSTER_REPS_OUTPUT_FILE, CLUSTERS_NT_SEQS_DIR, CLUSTERS_ALIGNMENTS_DIR


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


def perform_alignment_on_core_clusters():
    """Run MAFFT tool on fasta files of protein nucleotide seqs for each core cluster"""
    logger = logging.getLogger()
    logger.info("Running MAFFT on core clusters for alignment")
    if not os.path.exists(CLUSTERS_NT_SEQS_DIR):
        logger.error("No clusters dir found, exiting")
        exit(1)
    if not os.path.exists(CLUSTERS_ALIGNMENTS_DIR):
        os.makedirs(CLUSTERS_ALIGNMENTS_DIR)

    for cluster_file in os.listdir(CLUSTERS_NT_SEQS_DIR):
        cluster_alignment_file = cluster_file + "_alignment"
        mafft_args = " ".join(["mafft", "--auto", os.path.join(CLUSTERS_NT_SEQS_DIR, cluster_file)])
        mafft_return_code = run(mafft_args, shell=True,
                                stdout=open(os.path.join(CLUSTERS_ALIGNMENTS_DIR, cluster_alignment_file)),
                                stderr=open(os.path.join(CLUSTERS_ALIGNMENTS_DIR, cluster_alignment_file))).returncode
        logger.info("Finished running MAFFT for %s with return code %d" % (cluster_file, mafft_return_code))
    logger.info("Finished running MAFFT for all clusters")


def perform_pruning_on_alignments():
    """Run Gblocks tool on MAFFT alignments of protein nucleotide seqs for each core cluster"""
    logger = logging.getLogger()
    logger.info("Running Gblocks on core clusters alignments")
    if not os.path.exists(CLUSTERS_ALIGNMENTS_DIR):
        logger.error("No clusters alignments dir found, exiting")
        exit(1)

    for alignment_file in os.listdir(CLUSTERS_ALIGNMENTS_DIR):
        gblocks_args = " ".join(["Gblocks", os.path.join(CLUSTERS_ALIGNMENTS_DIR, alignment_file), "-t=d", "-b5=a"])
        gblocks_return_code = run(gblocks_args, shell=True).returncode
        logger.info("Finished running Gblocks for alignment %s with return code %d" % (alignment_file, gblocks_return_code))
    logger.info("Finished running Gblocks for all alignments")
