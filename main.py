import argparse
import multiprocessing
import sys
import logging
import os
from subprocess import run
import pandas

from data_visualization import create_1st_stage_charts
from nucleotide_preprocessor import create_representatives_and_pseudogenes_file
from constants import STRAINS_DIR, COMBINED_PROTEINS_FILE_PATH, CD_HIT_CLUSTER_REPS_OUTPUT_FILE, \
    CD_HIT_CLUSTERS_OUTPUT_FILE, CD_HIT_EST_CLUSTER_CDS_OUTPUT_FILE, COMBINED_CDS_FILE_PATH, \
    FIRST_STAGE_STATS_PKL
from data_analysis import get_1st_stage_stats_per_strain
from ftp_handler import download_strain_files
from logging_config import listener_process, listener_configurer, worker_configurer
from protein_preprocessor import create_all_strains_file_with_indices


def main():
    parser = init_args_parser()
    args = parser.parse_args()
    if not len(sys.argv) > 1:
        parser.print_help()
        exit(0)
    log_queue = multiprocessing.Queue(-1)
    listener = multiprocessing.Process(target=listener_process,
                                       args=(log_queue, listener_configurer))
    listener.start()
    worker_configurer(log_queue)
    logger = logging.getLogger()

    if not os.path.exists(STRAINS_DIR):
        os.makedirs(STRAINS_DIR)

    try:
        stats_df = None
        logger.info("Starting work")
        if args.download:
            download_strain_files(STRAINS_DIR, log_queue, sample_size=args.sample_size)
        if args.preprocess_proteins:
            if os.listdir(STRAINS_DIR):
                create_all_strains_file_with_indices(log_queue)
            else:
                logger.error("Cannot preprocess strain proteins without downloaded strains")
        if args.preprocess_cds:
            if os.listdir(STRAINS_DIR) and os.path.exists(CD_HIT_CLUSTERS_OUTPUT_FILE):
                create_representatives_and_pseudogenes_file(log_queue)
        if args.cluster_proteins:
            if os.path.exists(COMBINED_PROTEINS_FILE_PATH):
                perform_clustering_on_proteins(COMBINED_PROTEINS_FILE_PATH)
            else:
                logger.error("Cannot run clustering without pre-processed proteins file")
        if args.cluster_cds:
            if os.path.exists(COMBINED_CDS_FILE_PATH):
                perform_clustering_on_cds(COMBINED_CDS_FILE_PATH)
            else:
                logger.error("Cannot run clustering without pre-processed cds file")
        if args.stats:
            logger.info("Gathering genomic and 1st stage clusters statistics per strain")
            if os.path.exists(CD_HIT_CLUSTERS_OUTPUT_FILE):
                stats_df = get_1st_stage_stats_per_strain()
                stats_df.to_pickle(FIRST_STAGE_STATS_PKL)
            else:
                logger.error("Cannot perform analysis without clusters file")
        if args.graph:
            logger.info("Plotting charts from statistics")
            if stats_df is None:
                logger.info("retrieving 1st stage stats from pkl file")
                stats_df = pandas.read_pickle(FIRST_STAGE_STATS_PKL)
            create_1st_stage_charts(stats_df)

        logger.info("Finished work, exiting")
    finally:
        log_queue.put_nowait(None)
        listener.join()


def init_args_parser():
    parser = argparse.ArgumentParser(description='Data processing pipeline for pseudogene search '
                                                 'in Pseudomonas Areguinosa strains')
    parser.add_argument('-d', '--download', action="store_true", help='Download all valid PA strains from the refseq ftp for analysis')
    parser.add_argument('--sample', type=int, dest='sample_size', default=None,
                        help='Specify a sample size to limit the amount of strains downloaded')
    parser.add_argument('-p', '--preprocess_proteins', action="store_true", help='Preprocess downloaded PA strains proteins')
    parser.add_argument('-c', '--cluster_proteins', action="store_true", help='Run CD-HIT clustering on preprocessed PA strains proteins')
    parser.add_argument('-s', '--stats', action="store_true", help='Get stats from CD-HIT clustering output')
    parser.add_argument('-g', '--graph', action="store_true", help='Plot graphs from strain stats')
    parser.add_argument('-r', '--preprocess_cds', action="store_true",
                        help='Preprocess clustered PA strains proteins representative and pseudogene cds')
    parser.add_argument('-x', '--cluster_cds', action="store_true",
                        help='Run CD-HIT clustering on preprocessed PA strains cds of representatives and pseudogenes')
    return parser


def perform_clustering_on_proteins(aggregated_proteins_file_path):
    """Run the CD-HIT program to perform clustering on the strains"""
    logger = logging.getLogger()
    logger.info("Running CD-HIT on combined proteins file to create clustering")
    cd_hit_args = " ".join(["cd-hit", "-i", aggregated_proteins_file_path, "-o", CD_HIT_CLUSTER_REPS_OUTPUT_FILE, "-c 0.70",
                   "-n 5", "-M 16000", "-g 1", "-p 1"])
    cd_hit_return_code = run(cd_hit_args, shell=True).returncode
    logger.info("Finished running CD-HIT with return code %d" % cd_hit_return_code)
    return cd_hit_return_code


def perform_clustering_on_cds(aggregated_cds_file_path):
    """Run the CD-HIT-EST program to perform clustering on the strains representatives and pseudogenes"""
    logger = logging.getLogger()
    logger.info("Running CD-HIT-EST on combined representative and pseudogene cds file to create clustering")
    cd_hit_est_args = " ".join(["cd-hit-est", "-i", aggregated_cds_file_path, "-o", CD_HIT_EST_CLUSTER_CDS_OUTPUT_FILE, "-c 0.8",
                   "-n 5", "-M 16000", "-g 1", "-p 1"])
    cd_hit_est_return_code = run(cd_hit_est_args, shell=True).returncode
    logger.info("Finished running CD-HIT with return code %d" % cd_hit_est_return_code)
    return cd_hit_est_return_code


if __name__ == '__main__':
    sys.exit(main())
