import multiprocessing
import os
import logging
from ftplib import FTP, error_temp
from io import StringIO
from time import sleep

from logging_config import worker_configurer

NCBI_FTP_SITE = "ftp.ncbi.nlm.nih.gov"
PA_LATEST_REFSEQ_URL = "/genomes/refseq/bacteria/Pseudomonas_aeruginosa/latest_assembly_versions"

NUMBER_OF_PROCESSES = os.cpu_count()
ftp_handle = FTP(NCBI_FTP_SITE)
#TODO add syncing for already downloaded strains


def download_valid_strains(worker_id, job_queue, log_queue, configurer, download_dir, strains_downloaded_counter):
    """
    Filter out any strain that does not have a features_table / cds_from_genomic
    as well as all suppressed strains from the strain list
    """
    configurer(log_queue)
    logger = logging.getLogger(__name__ + "_worker_" + str(worker_id))
    ftp_con = FTP(NCBI_FTP_SITE)
    ftp_con.login()
    ftp_con.cwd(PA_LATEST_REFSEQ_URL)
    num_of_strains_downloaded = 0
    while True:
        strain_dir = job_queue.get()
        if strain_dir is None:
            job_queue.put(None)
            break
        try:
            strain_dir_files_list = ftp_con.nlst(strain_dir)
            feature_table = [file for file in strain_dir_files_list if FEATURE_TABLE_PATTERN in file][0]
            cds_from_genomic = [file for file in strain_dir_files_list if CDS_FROM_GENOMIC_PATTERN in file][0]
            genomic_sequences = [file for file in strain_dir_files_list if strain_dir + "_genomic.fna" in file][0]
            protein_sequences = [file for file in strain_dir_files_list if PROTEIN_FILE_PATTERN in file][0]
            if protein_sequences:
                if feature_table or cds_from_genomic:
                    status_file = [file for file in strain_dir_files_list if "assembly_status" in file]
                    download_strain = True
                    if status_file:
                        line_reader = StringIO()
                        ftp_con.retrlines('RETR ' + status_file[0], line_reader.write)
                        if "suppressed" in line_reader.getvalue():
                            logger.info("skipping suppressed strain %s" % strain_dir)
                            download_strain = False
                            line_reader.close()

                    if download_strain:
                        strain_download_dir = download_dir + os.sep + strain_dir + os.sep
                        if not os.path.exists(strain_download_dir):
                            os.mkdir(strain_download_dir)
                        if feature_table:
                            ftp_con.retrbinary("RETR " + feature_table,
                                               open(strain_download_dir + feature_table[len(strain_dir) + 1:], 'wb').write)
                        if cds_from_genomic:
                            ftp_con.retrbinary("RETR " + cds_from_genomic,
                                               open(strain_download_dir + cds_from_genomic[len(strain_dir) + 1:],
                                                    'wb').write)
                        if genomic_sequences:
                            ftp_con.retrbinary("RETR " + genomic_sequences,
                                               open(strain_download_dir + genomic_sequences[len(strain_dir) + 1:],
                                                    'wb').write)
                        if protein_sequences:
                            ftp_con.retrbinary("RETR " + protein_sequences,
                                               open(strain_download_dir + protein_sequences[len(strain_dir) + 1:],
                                                    'wb').write)
                        with strains_downloaded_counter.get_lock():
                            strain_index = strains_downloaded_counter.value
                            with open(os.path.join(strain_download_dir, STRAIN_INDEX_FILE), 'w') as index_file:
                                index_file.write(str(strain_index))
                            strains_downloaded_counter.value += 1
                        num_of_strains_downloaded += 1

                        logger.info("Downloaded files for strain %s" % strain_dir)
                else:
                    logger.info("No feature_table or cds_from_genomic files found for strain %s" % strain_dir)
            else:
                logger.info("No protein sequences found for strain %s" % strain_dir)
        except error_temp:
            job_queue.put(strain_dir)
            sleep(2)
    ftp_con.quit()
    exit(0)


def download_strain_files(download_dir, log_queue, sample_size=None):
    logger = logging.getLogger()
    ftp_handle.login()
    job_queue = multiprocessing.Queue()
    get_strains_from_ftp(job_queue, sample_size)
    logger.info("Starting download of strain files")
    strains_downloaded_counter = multiprocessing.Value('L', 0)
    workers = [multiprocessing.Process(target=download_valid_strains, args=(i, job_queue, log_queue, worker_configurer,
                                                                            download_dir, strains_downloaded_counter))
               for i in range(NUMBER_OF_PROCESSES)]
    for w in workers:
        w.start()
    job_queue.put(None)
    for w in workers:
        w.join()
    with strains_downloaded_counter.get_lock():
        logger.info("Finished downloading strain files, %s strains downloaded" % strains_downloaded_counter.value)
    job_queue.close()


def get_strains_from_ftp(job_queue, sample_size):
    logger = logging.getLogger()
    logger.info("Listing latest PA strains for download")
    ftp_handle.cwd(PA_LATEST_REFSEQ_URL)
    strains_dir_listing = ftp_handle.nlst()
    ftp_handle.quit()
    if sample_size:
        strains_dir_listing = strains_dir_listing[:sample_size]
    for strain_dir in strains_dir_listing:
        job_queue.put(strain_dir)



