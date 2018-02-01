import argparse
import gzip
import sys
import logging
import threading
import os

from multiprocessing import Queue
from Bio import SeqIO
from ftp_handler import FtpHandler

DATA_DIR = os.getcwd() + os.sep + "data"
STRAINS_DIR = DATA_DIR + os.sep + "strains"


def logger_thread(log_queue):
    while True:
        record = log_queue.get()
        if record is None:
            break
        logger = logging.getLogger(record.name)
        logger.handle(record)

# Add with statements to automatically close resources
def main():
    parser = argparse.ArgumentParser(description='Data processing pipeline for pseudogene search in Pseudomonas Areguinosa strains')
    parser.add_argument('-d', '--download', action="store_true",
                        help='Download all valid PA strains from the refseq ftp for analysis')
    parser.add_argument('-s', '--sample', type=int, dest='sample_size', default=None,
                        help='Specify a sample size to limit the amount of strains downloaded')
    args = parser.parse_args()

    if not os.path.exists(STRAINS_DIR):
        os.makedirs(STRAINS_DIR)

    log_queue = Queue()
    logger = logging.getLogger(__name__)
    logging.basicConfig(level=logging.DEBUG)
    handler = logging.FileHandler('pseudogene.log', 'w')
    handler.setLevel(logging.INFO)
    logger.addHandler(handler)
    logger.addHandler(logging.StreamHandler())
    lp = threading.Thread(target=logger_thread, args=(log_queue,))
    lp.start()

    if args.download:
        ftp_handler = FtpHandler()
        logger.info("Downloading Strain files from NCBI FTP (refseq)")
        ftp_handler.download_strain_files(STRAINS_DIR, log_queue, sample_size=args.sample_size)

    all_strains_file = create_all_strains_file_with_indices(logger)
    perform_clustering_on_strains(all_strains_file, logger)

    logger.info("Finished work, exiting")
    log_queue.put(None)
    lp.join()


def create_all_strains_file_with_indices(logger):
    """Add indices to each strain sequence header for easier parsing"""
    logger.info("Indexing proteins by strain index + location in gene index")
    output_file_path = os.path.join(DATA_DIR, 'all_strains_proteins.fasta')
    if os.path.exists(output_file_path):
        os.remove(output_file_path)
    for root, dirs, files in os.walk(STRAINS_DIR):
        for dir in dirs:
            dir_files = os.listdir(os.path.join(root, dir))
            protein_file_name = [f for f in dir_files if 'protein.faa' in f][0]
            cds_file_name = [f for f in dir_files if 'cds_from_genomic.fna' in f][0]
            if not protein_file_name or not cds_file_name:
                logger.warn("Could not find protein file or cds_from_genomic file for strain {}, skipping", dir)
                continue
            strain_index = dir[:dir.find(']') + 1]
            if protein_file_name.endswith('gz'):
                protein_file = gzip.open(os.path.join(root, dir, protein_file_name), 'rt')
            else:
                protein_file = open(os.path.join(root, dir, protein_file_name))
            if cds_file_name.endswith('gz'):
                cds_file = gzip.open(os.path.join(root, dir, cds_file_name), 'rt')
            else:
                cds_file = open(os.path.join(root, dir, cds_file_name))

            strain_cds_seq_dict = SeqIO.to_dict(SeqIO.parse(cds_file, 'fasta'))
            strain_cds_seq_headers = list(strain_cds_seq_dict.keys())
            strain_protein_seq_iter = SeqIO.parse(protein_file, 'fasta')
            for strain_protein_seq in strain_protein_seq_iter:
                protein_id = strain_protein_seq.id
                cds_protein_header = [h for h in strain_cds_seq_headers if protein_id in h][0]
                protein_index_in_gene = '[' + cds_protein_header[cds_protein_header.rfind('_') + 1:] + ']'
                strain_protein_seq.description = strain_index + protein_index_in_gene + strain_protein_seq.description
                strain_protein_seq.id = ""
                SeqIO.write(strain_protein_seq, open(output_file_path, 'a+'), 'fasta')
    logger.info("Finished adding all strains protein sequences to combined file")
    return output_file_path


def perform_clustering_on_strains(all_strains_file, logger):
    """Run the CD-HIT program to perform clustering on the strains"""
    logger.info("Running CD-HIT on combined proteins file to create clustering")
    clusters_file = os.path.join(DATA_DIR, 'protein_clusters.txt')
    cd_hit_args = ['cd-hit', '-i', all_strains_file, '-o', clusters_file, '-c 0.70', '-n 5', '-M 16000', '-g 1', '-p 1']
    os.subprocess.call(cd_hit_args)


if __name__ == '__main__':
    sys.exit(main())
