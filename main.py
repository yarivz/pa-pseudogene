import sys
from ftplib import FTP


def main():
    """Main entry point for the script."""
    get_ftp_files()


def get_ftp_files():
    ftp = FTP('ftp.ncbi.nlm.nih.gov')
    ftp.login()
    ftp.cwd("/genomes/refseq/bacteria/Pseudomonas_aeruginosa/latest_assembly_versions")
    files = ftp.dir()
    print(files)


if __name__ == '__main__':
    sys.exit(main())

