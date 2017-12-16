import sys
from ftplib import FTP

# Add retries for network operations since ftp can be slow to respond
#
#
#


def main():
    """Main entry point for the script."""
    all_strains = get_ftp_files()


def get_ftp_files():
    ftp = FTP('ftp.ncbi.nlm.nih.gov')
    ftp.login()
    ftp.cwd("/genomes/refseq/bacteria/Pseudomonas_aeruginosa/latest_assembly_versions")
    files = ftp.dir()
    print(files)


def filter_suppressed_strains():
    """Filter out all suppressed strains from the strain list"""
    #Look for an assembly_status.txt file with suppressed in it
    pass


def filter_missing_data():
    """Filter out any strain that does not have a features_table / cds_from_genomic from the strain list"""
    pass


def add_strain_indices():
    """Add indices to each strain sequence header for easier parsing"""
    pass


def perform_clustering_on_strains():
    """Run the CD-HIT program to perform clustering on the strains"""
    pass


if __name__ == '__main__':
    sys.exit(main())
