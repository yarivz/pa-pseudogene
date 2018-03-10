import argparse
import multiprocessing
import sys
import logging
import os
from subprocess import run

import pandas

from constants import STRAINS_DIR, COMBINED_PROTEINS_FILE_PATH, CD_HIT_CLUSTER_REPS_OUTPUT_FILE, \
    CD_HIT_CLUSTERS_OUTPUT_FILE
from data_analysis import get_strains_stats, get_genomic_stats_per_strain, create_strains_clusters_map
from ftp_handler import download_strain_files
from logging_config import listener_process, listener_configurer, worker_configurer
from protein_preprocessor import create_all_strains_file_with_indices
import matplotlib.pyplot as plt


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
        combined_proteins_file_path = None
        contigs = pseudogenes = total_clusters = core_clusters = singleton_clusters = None
        logger.info("Starting work")
        if args.download:
            download_strain_files(STRAINS_DIR, log_queue, sample_size=args.sample_size)
        if args.preprocess:
            if os.listdir(STRAINS_DIR):
                combined_proteins_file_path = create_all_strains_file_with_indices(log_queue)
            else:
                logger.error("Cannot preprocess strain proteins without downloaded strains")
        if args.cluster:
            if combined_proteins_file_path is not None:
                perform_clustering_on_strains(combined_proteins_file_path)
            elif os.path.exists(COMBINED_PROTEINS_FILE_PATH):
                perform_clustering_on_strains(COMBINED_PROTEINS_FILE_PATH)
            else:
                logger.error("Cannot run clustering without pre-processed proteins file")
        if args.stats:
            logger.info("Gathering statistics")
            genomic_stats, contigs, pseudogenes = get_genomic_stats_per_strain()
            for stat in genomic_stats:
                logger.info("Strain index: %s" % stat)
                logger.info(genomic_stats[stat])
            if os.path.exists(CD_HIT_CLUSTERS_OUTPUT_FILE):
                cluster_stats, total_clusters, core_clusters, singleton_clusters = get_strains_stats(CD_HIT_CLUSTERS_OUTPUT_FILE)
                for stat in cluster_stats:
                    logger.info("Strain index: %s" % stat)
                    logger.info(cluster_stats[stat])
            else:
                logger.error("Cannot perform analysis without clusters file")
        if args.graph:
            logger.info("Plotting charts from statistics")
            plt.use('Agg')
            strains_map, total_strains_count = create_strains_clusters_map(CD_HIT_CLUSTERS_OUTPUT_FILE)
            x_strains = y_clusters = []
            for index, strain in strains_map.items():
                for c in strain.containing_clusters.keys():
                    x_strains.append(index)
                    y_clusters.append(c)
            plt.scatter(x_strains, y_clusters).savefig('clusters_by_strain_scatterplot.png')
            if total_clusters:
                plt.hist(total_clusters).savefig('total_clusters_by_strain_index.png')
            if core_clusters:
                plt.hist(core_clusters).savefig('core_clusters_by_strain_index.png')
            if singleton_clusters:
                plt.hist(singleton_clusters).savefig('singleton_clusters_by_strain_index.png')
            if contigs:
                plt.hist(contigs).savefig('contigs_by_strain_index.png')
            if pseudogenes:
                plt.hist(pseudogenes).savefig('pseudogenes_by_strain_index.png')
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
    parser.add_argument('-p', '--preprocess', action="store_true", help='Preprocess downloaded PA strains proteins')
    parser.add_argument('-c', '--cluster', action="store_true", help='Run CD-HIT clustering on preprocessed PA strains proteins')
    parser.add_argument('-s', '--stats', action="store_true", help='Get stats from CD-HIT clustering output')
    parser.add_argument('-g', '--graph', action="store_true", help='Plot graphs from strain stats')
    return parser


def perform_clustering_on_strains(aggregated_proteins_file_path):
    """Run the CD-HIT program to perform clustering on the strains"""
    logger = logging.getLogger()
    logger.info("Running CD-HIT on combined proteins file to create clustering")
    cd_hit_args = " ".join(["cd-hit", "-i", aggregated_proteins_file_path, "-o", CD_HIT_CLUSTER_REPS_OUTPUT_FILE, "-c 0.70",
                   "-n 5", "-M 16000", "-g 1", "-p 1"])
    cd_hit_return_code = run(cd_hit_args, shell=True).returncode
    logger.info("Finished running CD-HIT with return code %d" % cd_hit_return_code)
    return cd_hit_return_code


if __name__ == '__main__':
    sys.exit(main())
