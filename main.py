import argparse
import multiprocessing
import sys
import logging
import os
from subprocess import run
import pickle
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from constants import STRAINS_DIR, COMBINED_PROTEINS_FILE_PATH, CD_HIT_CLUSTER_REPS_OUTPUT_FILE, \
    CD_HIT_CLUSTERS_OUTPUT_FILE, GENOMIC_STATS_PKL, PROTEIN_STATS_PKL, TOTAL_CLUSTERS_PKL, CORE_CLUSTERS_PKL, \
    SINGLETON_CLUSTERS_PKL, CONTIGS_PKL, PSEUDOGENES_PKL
from data_analysis import get_strains_stats, get_genomic_stats_per_strain, create_strains_clusters_map
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
            with open(GENOMIC_STATS_PKL, 'wb') as f:
                pickle.dump(genomic_stats, f)
            with open(CONTIGS_PKL, 'wb') as f:
                pickle.dump(contigs, f)
            with open(PSEUDOGENES_PKL, 'wb') as f:
                pickle.dump(pseudogenes, f)
            if os.path.exists(CD_HIT_CLUSTERS_OUTPUT_FILE):
                cluster_stats, total_clusters, core_clusters, singleton_clusters = get_strains_stats(CD_HIT_CLUSTERS_OUTPUT_FILE)
                with open(PROTEIN_STATS_PKL, 'wb') as f:
                    pickle.dump(cluster_stats, f)
                with open(TOTAL_CLUSTERS_PKL, 'wb') as f:
                    pickle.dump(total_clusters, f)
                with open(CORE_CLUSTERS_PKL, 'wb') as f:
                    pickle.dump(core_clusters, f)
                with open(SINGLETON_CLUSTERS_PKL, 'wb') as f:
                    pickle.dump(singleton_clusters, f)
            else:
                logger.error("Cannot perform analysis without clusters file")
        if args.graph:
            logger.info("Plotting charts from statistics")
            if not total_clusters:
                with open(TOTAL_CLUSTERS_PKL, 'rb') as f:
                    total_clusters = pickle.load(f)
            if not core_clusters:
                with open(CORE_CLUSTERS_PKL, 'rb') as f:
                    core_clusters = pickle.load(f)
            if not singleton_clusters:
                with open(SINGLETON_CLUSTERS_PKL, 'rb') as f:
                    singleton_clusters = pickle.load(f)
            if not contigs:
                with open(CONTIGS_PKL, 'rb') as f:
                    contigs = pickle.load(f)
            if not pseudogenes:
                with open(PSEUDOGENES_PKL, 'rb') as f:
                    pseudogenes = pickle.load(f)
            # plt.figure(figsize=(70, 70))
            ax = plt.subplot()  # Defines ax variable by creating an empty plot
            # Set the tick labels font
            for label in (ax.get_xticklabels() + ax.get_yticklabels()):
                label.set_fontsize(1)
            strains_map, total_strains_count = create_strains_clusters_map(CD_HIT_CLUSTERS_OUTPUT_FILE)
            x_strains = []
            y_clusters = []
            logger.info("Plotting strains to clusters scatter chart")
            for index, strain in strains_map.items():
                for c in strain.containing_clusters.keys():
                    x_strains.append(index)
                    y_clusters.append(c)
            plt.scatter(x_strains, y_clusters)
            plt.xlabel("strains (indices)")
            plt.ylabel("clusters (indices)")
            plt.title("strains to clusters heatmap")
            plt.savefig('clusters_by_strain_scatterplot.pdf', format="pdf")
            plt.close()
            if total_clusters:
                logger.info("Plotting strains to clusters histogram")
                plt.hist(total_clusters, color='green')
                plt.ylabel("strains #")
                plt.xlabel("clusters #")
                plt.title("strains to clusters histogram")
                plt.xticks(list(range(4000, 7300, 100)))
                plt.savefig('total_clusters_by_strain_index.pdf', format="pdf")
                plt.close()
            if core_clusters:
                logger.info("Plotting strains to core clusters histogram")
                plt.hist(core_clusters, color='green')
                plt.ylabel("strains #")
                plt.xlabel("clusters #")
                plt.title("strains to core clusters histogram")
                plt.xticks(list(range(4500, 5300, 50)))
                plt.savefig('core_clusters_by_strain_index.pdf', format="pdf")
                plt.close()
            if singleton_clusters:
                logger.info("Plotting strains to singleton clusters histogram")
                plt.hist(singleton_clusters, color='green')
                plt.ylabel("strains #")
                plt.xlabel("clusters #")
                plt.title("strains to singleton clusters histogram")
                plt.xticks(list(range(0, 50, 1)) + list(range(50, 250, 20)))
                plt.savefig('singleton_clusters_by_strain_index.pdf', format="pdf")
                plt.close()
            if contigs:
                logger.info("Plotting strains to contigs histogram")
                plt.hist(contigs, color='green')
                plt.ylabel("strains #")
                plt.xlabel("contigs #")
                plt.title("strains to contigs histogram")
                plt.xticks(list(range(0, 10, 1)) + list(range(20, 250, 10)) + list(range(250, 2000, 100)))
                plt.savefig('contigs_by_strain_index.pdf', format="pdf")
                plt.close()
            if pseudogenes:
                logger.info("Plotting strains to pseudogenes histogram")
                plt.hist(pseudogenes, color='green')
                plt.ylabel("strains #")
                plt.xlabel("pseudogenes #")
                plt.title("strains to pseudogenes histogram")
                plt.xticks(list(range(0, 10, 1)) + list(range(20, 250, 10)) + list(range(250, 2000, 100)))
                plt.savefig('pseudogenes_by_strain_index.pdf', format="pdf")
                plt.close()
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
