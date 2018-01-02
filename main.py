import os
import sys
from ftplib import FTP
from io import StringIO

NCBI_FTP_SITE = "ftp.ncbi.nlm.nih.gov"
PA_LATEST_REFSEQ_URL = "/genomes/refseq/bacteria/Pseudomonas_aeruginosa/latest_assembly_versions"
TARGET_DIR = os.getcwd() + os.sep + "data" + os.sep + "strains"


def main():
    """Main entry point for the script."""
    if not os.path.exists(TARGET_DIR):
        os.makedirs(TARGET_DIR)
    ftp = FTP(NCBI_FTP_SITE)
    ftp.login()
    strains_dir_listing = get_ftp_root_dir_listing(ftp)
    download_filtered_strains(ftp, strains_dir_listing)
    ftp.quit()


def get_ftp_root_dir_listing(ftp):
    ftp.cwd(PA_LATEST_REFSEQ_URL)
    return ftp.nlst()


def download_filtered_strains(ftp, strains_dir_listing):
    """
    Filter out any strain that does not have a features_table / cds_from_genomic from the strain list
    Filter out all suppressed strains from the strain list
    """

    for strain_dir in strains_dir_listing:
        strain_dir_files_list = ftp.nlst(strain_dir)
        feature_table = [file for file in strain_dir_files_list if strain_dir + "_feature_table.txt" in file][0]
        cds_from_genomic = [file for file in strain_dir_files_list if strain_dir + "_cds_from_genomic.fna" in file][0]
        genomic_sequences = [file for file in strain_dir_files_list if strain_dir + "_genomic.fna" in file][0]
        protein_sequences = [file for file in strain_dir_files_list if strain_dir + "_protein.faa" in file][0]
        if feature_table or cds_from_genomic:
            status_file = [file for file in strain_dir_files_list if "assembly_status" in file]
            if status_file:
                line_reader = StringIO()
                ftp.retrlines('RETR ' + status_file[0], line_reader.write)
                if "suppressed" in line_reader.getvalue():
                    print("skipping suppressed strain %s" % strain_dir)
                    line_reader.close()
                else:
                    strain_download_dir = TARGET_DIR + os.sep + strain_dir + os.sep
                    os.mkdir(strain_download_dir)
                    if feature_table:
                        ftp.retrbinary("RETR " + feature_table,
                                       open(strain_download_dir + feature_table[len(strain_dir) + 1:], 'wb').write)
                    if cds_from_genomic:
                        ftp.retrbinary("RETR " + cds_from_genomic,
                                       open(strain_download_dir + cds_from_genomic[len(strain_dir) + 1:], 'wb').write)
                    if genomic_sequences:
                        ftp.retrbinary("RETR " + genomic_sequences,
                                       open(strain_download_dir + genomic_sequences[len(strain_dir) + 1:], 'wb').write)
                    if protein_sequences:
                        ftp.retrbinary("RETR " + protein_sequences,
                                       open(strain_download_dir + protein_sequences[len(strain_dir) + 1:], 'wb').write)
                    print("Downloaded files for strain %s" % strain_dir)
        else:
            print("No feature_table or cds_from_genomic files found for strain %s" % strain_dir)


def add_strain_indices():
    """Add indices to each strain sequence header for easier parsing"""
    pass


def perform_clustering_on_strains():
    """Run the CD-HIT program to perform clustering on the strains"""
    pass


if __name__ == '__main__':
    sys.exit(main())
