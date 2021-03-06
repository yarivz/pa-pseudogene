import logging
import multiprocessing
import os
from subprocess import run

from Bio import SeqIO, AlignIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from constants import CD_HIT_CLUSTER_REPS_OUTPUT_FILE, CLUSTERS_NT_SEQS_DIR, CLUSTERS_ALIGNMENTS_DIR, \
    NUMBER_OF_PROCESSES, FASTA_FILE_TYPE, ALIGNMENTS_FOR_TREE_DIR, DATA_DIR, ALIGNMENT_STRAIN_PATTERN, STRAINS_COUNT
from data_analysis import build_strain_names_map
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
        alignment_stdout = open("alignment_stdout.log", "w")
        alignment_stderr = open("alignment_stderr.log", "w")
        cluster_alignment_filename = cluster_file + "_alignment"
        if not os.path.exists(os.path.join(CLUSTERS_ALIGNMENTS_DIR, cluster_alignment_filename)):
            cluster_alignment_file = open(os.path.join(CLUSTERS_ALIGNMENTS_DIR, cluster_alignment_filename), 'w')
            mafft_args = " ".join(["mafft", "--auto", os.path.join(CLUSTERS_NT_SEQS_DIR, cluster_file)])
            mafft_return_code = run(mafft_args, shell=True, stdout=cluster_alignment_file, stderr=alignment_stderr).returncode
            logger.info("Finished running MAFFT for %s with return code %d" % (cluster_file, mafft_return_code))
            cluster_alignment_file.close()

        logger.info("Running GBlocks for %s" % cluster_file)
        gblocks_args = " ".join(["Gblocks", os.path.join(CLUSTERS_ALIGNMENTS_DIR, cluster_alignment_filename), "-t=d", "-b5=a", "-p=n"])
        gblocks_return_code = run(gblocks_args, shell=True, stdout=alignment_stdout, stderr=alignment_stderr).returncode
        logger.info(
            "Finished running Gblocks for alignment %s with return code %d" % (cluster_alignment_filename, gblocks_return_code))


def prepare_alignments_for_tree(log_queue):
    """Edit each alignment to remove invariant positions, pad missing strain seqs & concatenate all alignments"""
    logger = logging.getLogger(__name__)
    logger.info("Preparing core clusters alignments for tree")
    if not os.path.exists(CLUSTERS_ALIGNMENTS_DIR):
        logger.error("No alignments dir found, exiting")
        exit(1)
    if not os.path.exists(ALIGNMENTS_FOR_TREE_DIR):
        os.makedirs(ALIGNMENTS_FOR_TREE_DIR)

    job_queue = multiprocessing.Queue()
    prepare_alignment_editing_jobs(job_queue)
    workers = [
        multiprocessing.Process(target=perform_alignment_editing, args=(i, job_queue, worker_configurer, log_queue))
        for i in range(NUMBER_OF_PROCESSES)]
    for w in workers:
        w.start()
    job_queue.put(None)
    for w in workers:
        w.join()
    logger.info("Finished editing all alignments, concatenating")
    edited_alignment_files = os.listdir(ALIGNMENTS_FOR_TREE_DIR)
    concatenated_alignment = None
    concatenated_alignment_file = os.path.join(DATA_DIR, "all_alignments")
    for edited_alignment_file in edited_alignment_files:
        logger.info("Concatenating alignment %s" % edited_alignment_file)
        with open(os.path.join(ALIGNMENTS_FOR_TREE_DIR, edited_alignment_file), "r") as f:
            edited_alignment = AlignIO.read(f, FASTA_FILE_TYPE)
            if not concatenated_alignment:
                concatenated_alignment = edited_alignment[:, :]
            else:
                concatenated_alignment += edited_alignment[:, :]
    AlignIO.write(concatenated_alignment, open(concatenated_alignment_file, "w"), FASTA_FILE_TYPE)
    logger.info("Finished concatenating all alignments, written to %s" % concatenated_alignment_file)


def prepare_alignment_editing_jobs(job_queue):
    """Put all downloaded strain dirs in job queue for workers"""
    alignments = os.listdir(CLUSTERS_ALIGNMENTS_DIR)
    for alignment_file in alignments:
        if alignment_file.endswith("-gb"):
            job_queue.put(alignment_file)


def perform_alignment_editing(worker_id, job_queue, configurer, log_queue):
    """
    Perform alignment editing
    """
    configurer(log_queue)
    logger = logging.getLogger(__name__ + "_worker_" + str(worker_id))
    while True:
        alignment_file = job_queue.get()
        if alignment_file is None:
            job_queue.put(None)
            break
        logger.info("Editing alignment %s" % alignment_file)
        alignment = AlignIO.read(open(os.path.join(CLUSTERS_ALIGNMENTS_DIR, alignment_file), "r"), FASTA_FILE_TYPE)
        edited_alignment = None
        for col_idx in range(alignment.get_alignment_length()):
            col = alignment[:, col_idx:col_idx + 1]
            col_str = alignment[:, col_idx]
            if not all(c == col_str[0] for c in col_str):
                if not edited_alignment:
                    edited_alignment = col
                else:
                    edited_alignment += col
        alignment_seq_len = edited_alignment.get_alignment_length()
        logger.info("alignment_seq_len = %d" % alignment_seq_len)
        strain_idx = 0
        while strain_idx < STRAINS_COUNT:
            logger.info("in while - strain_idx = %d" % strain_idx)
            if len(edited_alignment) > strain_idx:
                seq = edited_alignment[strain_idx]
                seq_strain_idx = int(ALIGNMENT_STRAIN_PATTERN.match(seq.id).group(1))
                logger.info("checking if strain idx %d < seq_strain_idx %d" % (strain_idx, seq_strain_idx))
                if strain_idx < seq_strain_idx:
                    for i in range(seq_strain_idx - strain_idx):
                        logger.info("adding padded seq at idx %d" % (strain_idx + i))
                        edited_alignment._records.insert(strain_idx + i, SeqRecord(Seq(alignment_seq_len * '-'), id="[%d] padding" % (strain_idx + i)))
                    strain_idx += (seq_strain_idx - strain_idx + 1)
                    continue
                strain_idx += 1
            else:
                logger.info("adding padded seq at end of alignment list")
                edited_alignment.append(SeqRecord(Seq(alignment_seq_len * '-'), id="[%d] padding" % strain_idx))
                strain_idx += 1
        alignment_file_edited = os.path.join(ALIGNMENTS_FOR_TREE_DIR, alignment_file)
        logger.info("Finished padding alignment - writing to file %s" % alignment_file_edited)
        AlignIO.write(edited_alignment, open(alignment_file_edited, "w"), FASTA_FILE_TYPE)


def format_concatenated_alignment():
    logger = logging.getLogger(__name__)
    strain_names_map = build_strain_names_map()
    tree_alignment = AlignIO.read(open(os.path.join(DATA_DIR, "all_alignments"), "r"), FASTA_FILE_TYPE)
    tree_alignment_filtered = AlignIO.MultipleSeqAlignment([])
    for id, strain in zip(range(STRAINS_COUNT), tree_alignment):
        if all(c == '-' for c in strain.seq):
            logger.info("skipping filtered strain %d" % id)
        else:
            logger.info("adding id to strain %d" % id)
            strain.id = "[" + str(id) + "]" + strain_names_map[id]
            strain.description = ''
            tree_alignment_filtered.append(strain)
    AlignIO.write(tree_alignment_filtered, open(os.path.join(DATA_DIR, "filtered_tree_alignment"), "w"), FASTA_FILE_TYPE)




