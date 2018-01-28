import argparse
import sys
import logging
import threading
import os
from glob import iglob

from multiprocessing import Queue
from ftp_handler import FtpHandler
from Bio import SeqIO

TARGET_DIR = os.getcwd() + os.sep + "data" + os.sep + "strains"


def logger_thread(log_queue):
    while True:
        record = log_queue.get()
        if record is None:
            break
        logger = logging.getLogger(record.name)
        logger.handle(record)


def main():
    parser = argparse.ArgumentParser(description='Data processing pipeline for pseudogene search in Pseudomonas Areguinosa strains')
    parser.add_argument('-d', '--download', action="store_true",
                        help='Download all valid PA strains from the refseq ftp for analysis')
    parser.add_argument('-s', '--sample', type=int, dest='sample_size', default=None,
                        help='Specify a sample size to limit the amount of strains downloaded')
    args = parser.parse_args()

    if not os.path.exists(TARGET_DIR):
        os.makedirs(TARGET_DIR)

    log_queue = Queue()
    logger = logging.getLogger(__name__)
    logging.basicConfig(level=logging.DEBUG)
    handler = logging.FileHandler('pseudogene.log', 'w')
    handler.setLevel(logging.INFO)
    logger.addHandler(handler)
    lp = threading.Thread(target=logger_thread, args=(log_queue,))
    lp.start()

    if args.download:
        ftp_handler = FtpHandler()
        logger.info("Downloading Strain files from NCBI FTP (refseq)")
        ftp_handler.download_strain_files(TARGET_DIR, log_queue, sample_size=args.sample_size)

    add_strain_indices()
    perform_clustering_on_strains()

    logger.info("Finished work, exiting")
    log_queue.put(None)
    lp.join()


def add_strain_indices():
    """Add indices to each strain sequence header for easier parsing"""
    for root, dirs, files in os.walk(TARGET_DIR):
        for dir in dirs:
            protein_file = [f for f in os.listdir(dir) if os.path.isfile(f) and f.endswith('protein.faa')]
            cds_file = [f for f in os.listdir(dir) if os.path.isfile(f) and f.endswith('cds_from_genomic.fna')]
            fasta_sequences = SeqIO.parse(open(protein_file), 'fasta')
            for sequence in fasta_sequences:
                print(sequence)


def perform_clustering_on_strains():
    """Run the CD-HIT program to perform clustering on the strains"""
    pass


if __name__ == '__main__':
    sys.exit(main())
