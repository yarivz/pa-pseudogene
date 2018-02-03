import argparse
import gzip
import multiprocessing
import sys
import logging
import os

from multiprocessing import Queue
from Bio import SeqIO
from ftp_handler import download_strain_files
from logging_config import listener_process, listener_configurer, worker_configurer

DATA_DIR = os.getcwd() + os.sep + "data"
STRAINS_DIR = DATA_DIR + os.sep + "strains"


def main():
    parser = argparse.ArgumentParser(description='Data processing pipeline for pseudogene search '
                                                 'in Pseudomonas Areguinosa strains')
    parser.add_argument('-d', '--download', action="store_true",
                        help='Download all valid PA strains from the refseq ftp for analysis')
    parser.add_argument('-s', '--sample', type=int, dest='sample_size', default=None,
                        help='Specify a sample size to limit the amount of strains downloaded')
    parser.add_argument('-c', '--cluster', action="store_true",
                        help='Run CD-HIT clustering on downloaded PA strains proteins')
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
        logger.info("Starting work")
        if args.download:
            download_strain_files(STRAINS_DIR, log_queue, sample_size=args.sample_size)
        if args.cluster:
            aggregated_proteins_file_path = create_all_strains_file_with_indices()
            perform_clustering_on_strains(aggregated_proteins_file_path)
        logger.info("Finished work, exiting")
    finally:
        log_queue.put_nowait(None)
        listener.join()


def create_all_strains_file_with_indices():
    """
    Add indices to each strain sequence header for easier parsing
    and combine all strains protein sequences into single fasta file for clustering
    """
    logger = logging.getLogger()
    logger.info("Indexing proteins by their strain index & protein index in strain gene")
    aggregated_proteins_file_path = os.path.join(DATA_DIR, 'all_strains_proteins.fasta')
    if os.path.exists(aggregated_proteins_file_path):
        os.remove(aggregated_proteins_file_path)
    for root, strain_dirs, files in os.walk(STRAINS_DIR):
        for strain_dir in strain_dirs:
            strain_dir_files = os.listdir(os.path.join(root, strain_dir))
            protein_file_name = [f for f in strain_dir_files if 'protein.faa' in f][0]
            cds_file_name = [f for f in strain_dir_files if 'cds_from_genomic.fna' in f][0]
            if not protein_file_name or not cds_file_name:
                logger.warning("Could not find protein file or cds_from_genomic file for strain %s, skipping" % strain_dir)
                continue
            protein_file = cds_file = aggregated_proteins_file = strain_index_file = None
            try:
                strain_index_file = open(os.path.join(root, strain_dir, 'strain_index'))
                strain_index = '[' + strain_index_file.readline() + ']'
                if protein_file_name.endswith('gz'):
                    protein_file = gzip.open(os.path.join(root, strain_dir, protein_file_name), 'rt')
                else:
                    protein_file = open(os.path.join(root, strain_dir, protein_file_name))
                if cds_file_name.endswith('gz'):
                    cds_file = gzip.open(os.path.join(root, strain_dir, cds_file_name), 'rt')
                else:
                    cds_file = open(os.path.join(root, strain_dir, cds_file_name))

                strain_cds_seq_dict = SeqIO.to_dict(SeqIO.parse(cds_file, 'fasta'))
                strain_cds_seq_headers = list(strain_cds_seq_dict.keys())
                strain_protein_seq_iter = SeqIO.parse(protein_file, 'fasta')
                for strain_protein_seq in strain_protein_seq_iter:
                    protein_id = strain_protein_seq.id
                    cds_protein_header = [h for h in strain_cds_seq_headers if protein_id in h][0]
                    protein_index_in_gene = '[' + cds_protein_header[cds_protein_header.rfind('_') + 1:] + ']'
                    strain_protein_seq.description = strain_index + protein_index_in_gene + strain_protein_seq.description
                    strain_protein_seq.id = ""
                    aggregated_proteins_file = open(aggregated_proteins_file_path, 'a+')
                    SeqIO.write(strain_protein_seq, aggregated_proteins_file, 'fasta')
                logger.info("Strain %s proteins were indexed and written to file" % strain_dir[strain_dir.rfind(']') + 1:])
            finally:
                if protein_file is not None:
                    protein_file.close()
                if cds_file is not None:
                    cds_file.close()
                if aggregated_proteins_file is not None:
                    aggregated_proteins_file.close()
                if strain_index_file is not None:
                    strain_index_file.close()
    logger.info("Finished adding all strains protein sequences to combined file")
    return aggregated_proteins_file_path


def perform_clustering_on_strains(aggregated_proteins_file_path):
    """Run the CD-HIT program to perform clustering on the strains"""
    logger = logging.getLogger()
    logger.info("Running CD-HIT on combined proteins file to create clustering")
    clusterring_output_file = os.path.join(DATA_DIR, 'protein_clusters.txt')
    cd_hit_args = ['cd-hit', '-i', aggregated_proteins_file_path, '-o', clusterring_output_file, '-c 0.70',
                   '-n 5', '-M 16000', '-g 1', '-p 1']


if __name__ == '__main__':
    sys.exit(main())
