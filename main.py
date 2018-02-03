import argparse
import multiprocessing
import sys
import logging
import os
from subprocess import call
from constants import STRAINS_DIR, DATA_DIR
from ftp_handler import download_strain_files
from logging_config import listener_process, listener_configurer, worker_configurer, log_queue
from protein_preprocessor import create_all_strains_file_with_indices


def main():
    parser = init_args_parser()
    args = parser.parse_args()
    if not len(sys.argv) > 1:
        parser.print_help()
        exit(0)

    listener = multiprocessing.Process(target=listener_process,
                                       args=(log_queue, listener_configurer))
    listener.start()
    worker_configurer(log_queue)
    logger = logging.getLogger()

    if not os.path.exists(STRAINS_DIR):
        os.makedirs(STRAINS_DIR)

    try:
        logger.info("Starting work")
        if args.download:
            download_strain_files(STRAINS_DIR, sample_size=args.sample_size)
        if args.cluster:
            aggregated_proteins_file_path = create_all_strains_file_with_indices()
            perform_clustering_on_strains(aggregated_proteins_file_path)
        logger.info("Finished work, exiting")
    finally:
        log_queue.put_nowait(None)
        listener.join()


def init_args_parser():
    parser = argparse.ArgumentParser(description='Data processing pipeline for pseudogene search '
                                                 'in Pseudomonas Areguinosa strains')
    parser.add_argument('-d', '--download', action="store_true",
                        help='Download all valid PA strains from the refseq ftp for analysis')
    parser.add_argument('-s', '--sample', type=int, dest='sample_size', default=None,
                        help='Specify a sample size to limit the amount of strains downloaded')
    parser.add_argument('-c', '--cluster', action="store_true",
                        help='Run CD-HIT clustering on downloaded PA strains proteins')
    return parser


def perform_clustering_on_strains(aggregated_proteins_file_path):
    """Run the CD-HIT program to perform clustering on the strains"""
    logger = logging.getLogger()
    logger.info("Running CD-HIT on combined proteins file to create clustering")
    clustering_output_file = os.path.join(DATA_DIR, 'protein_clusters.txt')
    cd_hit_args = ['cd-hit', '-i', aggregated_proteins_file_path, '-o', clustering_output_file, '-c 0.70',
                   '-n 5', '-M 16000', '-g 1', '-p 1']
    call(cd_hit_args)


if __name__ == '__main__':
    sys.exit(main())
