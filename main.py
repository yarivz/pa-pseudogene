import argparse
import multiprocessing
import sys
import logging
import os
from subprocess import run
import pandas

from data_visualization import create_1st_stage_charts, create_2nd_stage_charts
from nucleotide_preprocessor import create_representatives_and_pseudogenes_file
from constants import STRAINS_DIR, COMBINED_PROTEINS_FILE_PATH, CD_HIT_CLUSTER_REPS_OUTPUT_FILE, \
    CD_HIT_CLUSTERS_OUTPUT_FILE, CD_HIT_EST_CLUSTER_REPS_OUTPUT_FILE, COMBINED_CDS_FILE_PATH, \
    FIRST_STAGE_STATS_PKL, SECOND_STAGE_STRAIN_STATS_PKL, SECOND_STAGE_CLUSTER_STATS_PKL, FIRST_STAGE_STATS_CSV, \
    CD_HIT_EST_CLUSTERS_OUTPUT_FILE, SECOND_STAGE_AGGREGATED_CLUSTER_STATS_PKL, SECOND_STAGE_STATS_CSV
from data_analysis import get_1st_stage_stats_per_strain, get_2nd_stage_stats_per_strain, \
    get_2nd_stage_stats_per_cluster, filter_2nd_stage_clusters_with_multiple_proteins, \
    split_2nd_stage_combined_fasta_to_reps_pseudogenes, get_pseudogenes_without_blast_hits_fasta
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
        stats_df = strains_df = clusters_df = None
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
            input_file = args.input if args.input else COMBINED_CDS_FILE_PATH
            output_file = args.output if args.output else CD_HIT_EST_CLUSTER_REPS_OUTPUT_FILE
            if os.path.exists(input_file):
                perform_clustering_on_cds(input_file, output_file)
            else:
                logger.error("Cannot run clustering without pre-processed cds file")
        if args.protein_stats:
            logger.info("Gathering genomic and 1st stage clusters statistics per strain")
            if os.path.exists(CD_HIT_CLUSTERS_OUTPUT_FILE):
                stats_df = get_1st_stage_stats_per_strain()
                stats_df.to_pickle(FIRST_STAGE_STATS_PKL)
            else:
                logger.error("Cannot perform analysis without clusters file")
        if args.nucleotide_stats:
            logger.info("Gathering 2nd stage clusters statistics per strain")
            if stats_df is None:
                logger.info("retrieving 1st stage stats from pkl file")
                stats_df = pandas.read_pickle(FIRST_STAGE_STATS_PKL)
            strains_df, clusters_df = get_2nd_stage_stats_per_strain(stats_df)
            strains_df.to_pickle(SECOND_STAGE_STRAIN_STATS_PKL)
            clusters_df.to_pickle(SECOND_STAGE_CLUSTER_STATS_PKL)
        if args.graph_1st_stage:
            logger.info("Plotting charts from 1st stage statistics")
            if stats_df is None:
                logger.info("retrieving 1st stage stats from pkl file")
                stats_df = pandas.read_pickle(FIRST_STAGE_STATS_PKL)
            create_1st_stage_charts(stats_df)
        if args.graph_2nd_stage:
            logger.info("Plotting charts from 2nd stage statistics")
            if strains_df is None:
                logger.info("retrieving 1st stage stats from pkl file")
                strains_df = pandas.read_pickle(SECOND_STAGE_STRAIN_STATS_PKL)
            if clusters_df is None:
                logger.info("retrieving 1st stage stats from pkl file")
                clusters_df = pandas.read_pickle(SECOND_STAGE_CLUSTER_STATS_PKL)
            create_2nd_stage_charts(strains_df, clusters_df)
        if args.get_1st_stage_stats_csv:
            logger.info("Gathering genomic and protein clusters statistics per strain as CSV file")
            if os.path.exists(CD_HIT_CLUSTERS_OUTPUT_FILE):
                stats_df = get_1st_stage_stats_per_strain()
                stats_df.to_pickle(FIRST_STAGE_STATS_PKL)
            else:
                logger.error("Cannot perform analysis without clusters file")
            stats_df.to_csv(FIRST_STAGE_STATS_CSV)
        if args.get_2nd_stage_stats_csv:
            logger.info("Gathering 2nd stage clusters statistics per cluster as CSV file")
            if os.path.exists(CD_HIT_EST_CLUSTERS_OUTPUT_FILE):
                cluster_stats = get_2nd_stage_stats_per_cluster()
                cluster_stats.to_pickle(SECOND_STAGE_AGGREGATED_CLUSTER_STATS_PKL)
            else:
                logger.error("Cannot perform analysis without clusters file")
            cluster_stats.to_csv(SECOND_STAGE_STATS_CSV)
        if args.filter_clusters:
            filter_2nd_stage_clusters_with_multiple_proteins()
        if args.split_2nd_stage_fasta:
            split_2nd_stage_combined_fasta_to_reps_pseudogenes()
        if args.get_pseudogenes_no_hits_fasta:
            get_pseudogenes_without_blast_hits_fasta()

        logger.info("Finished work, exiting")
    finally:
        log_queue.put_nowait(None)
        listener.join()


def init_args_parser():
    parser = argparse.ArgumentParser(description='Data processing pipeline for pseudogene search '
                                                 'in Pseudomonas Areguinosa strains')
    parser.add_argument('-dl', '--download', action="store_true", help='Download all valid PA strains from the refseq ftp for analysis')
    parser.add_argument('--sample', type=int, dest='sample_size', default=None,
                        help='Specify a sample size to limit the amount of strains downloaded')
    parser.add_argument('-p', '--preprocess_proteins', action="store_true", help='Preprocess downloaded PA strains proteins')
    parser.add_argument('-c', '--cluster_proteins', action="store_true", help='Run CD-HIT clustering on preprocessed PA strains proteins')
    parser.add_argument('-s1', '--protein_stats', action="store_true", help='Get stats from CD-HIT clustering output')
    parser.add_argument('-s1csv', '--get_1st_stage_stats_csv', action="store_true", help='Get stage 1 stats in csv')
    parser.add_argument('-s2', '--nucleotide_stats', action="store_true", help='Get stats from CD-HIT-EST clustering output')
    parser.add_argument('-s2csv', '--get_2nd_stage_stats_csv', action="store_true", help='Get stage 2 nucleotide clusters stats in csv')
    parser.add_argument('-g1', '--graph_1st_stage', action="store_true", help='Plot graphs from 1st stage strain stats')
    parser.add_argument('-g2', '--graph_2nd_stage', action="store_true", help='Plot graphs from 2nd stage strain stats')
    parser.add_argument('-r', '--preprocess_cds', action="store_true",
                        help='Preprocess clustered PA strains proteins representative and pseudogene cds')
    parser.add_argument('-x', '--cluster_cds', action="store_true",
                        help='Run CD-HIT clustering on preprocessed PA strains cds of representatives and pseudogenes')
    parser.add_argument('-f', '--filter_clusters', action="store_true",
                        help='Filter 2nd stage clusters with multiple proteins')
    parser.add_argument('-sp', '--split_2nd_stage_fasta', action="store_true",
                        help='Split 2nd stage combined fasta to representatives and pseudogenes files')
    parser.add_argument('-pnh', '--get_pseudogenes_no_hits_fasta', action="store_true",
                        help='Get pseudogenes without blast hits fasta')
    parser.add_argument('-in', '--input', nargs=1, help='Get input file')
    parser.add_argument('-out', '--output', nargs=1, help='Get output file')
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


def perform_clustering_on_cds(input_file, output_file):
    """Run the CD-HIT-EST program to perform clustering on the strains representatives and pseudogenes"""
    logger = logging.getLogger()
    logger.info("Running CD-HIT-EST on combined representative and pseudogene cds file to create clustering")
    cd_hit_est_args = " ".join(["cd-hit-est", "-i", input_file, "-o", output_file, "-c 0.8",
                   "-n 5", "-M 16000", "-g 1", "-p 1", "-d 30"])
    cd_hit_est_return_code = run(cd_hit_est_args, shell=True).returncode
    logger.info("Finished running CD-HIT with return code %d" % cd_hit_est_return_code)
    return cd_hit_est_return_code


if __name__ == '__main__':
    sys.exit(main())
