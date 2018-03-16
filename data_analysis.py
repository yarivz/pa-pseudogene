import gzip
from collections import defaultdict

import os

from constants import STRAINS_DIR, CDS_FROM_GENOMIC_PATTERN, GENOMIC_PATTERN, STRAIN_INDEX_FILE, CLUSTER_STRAIN_PATTERN


class Cluster:
    def __init__(self, index):
        self.index = index
        self.member_strains = defaultdict(int)
        self.total_members = 0

    def add_strain(self, strain_index):
        self.member_strains[strain_index] += 1
        self.total_members += 1

    def get_cluster_strains_num(self):
        return len(self.member_strains)


class Strain:
    def __init__(self, index):
        self.index = index
        self.containing_clusters = {}

    def add_cluster(self, cluster):
        self.containing_clusters[cluster.index] = cluster

    def get_strain_core_clusters(self, total_strains_count):
        return [c for c in self.containing_clusters.values() if ((c.get_cluster_strains_num() / total_strains_count) >= 0.9)]

    def get_strain_singleton_clusters(self):
        return [c for c in self.containing_clusters.values() if c.get_cluster_strains_num() is 1]


def get_strains_stats(clusters_file):
    strains_map, total_strains_count = create_strains_clusters_map(clusters_file)
    return get_cluster_stats_per_strain(strains_map, total_strains_count)


def get_cluster_stats_per_strain(strains_map, total_strains_count):
    stats = {}
    total_clusters_to_strains = []
    core_clusters_to_strains = []
    singleton_clusters_to_strains = []
    for strain in strains_map.values():
        total_clusters = len(strain.containing_clusters)
        total_clusters_to_strains.append(total_clusters)
        core_clusters = len(strain.get_strain_core_clusters(total_strains_count))
        core_clusters_to_strains.append(core_clusters)
        singleton_clusters = len(strain.get_strain_singleton_clusters())
        singleton_clusters_to_strains.append(singleton_clusters)
        stats[strain.index] = "Total clusters: " + str(total_clusters) + ", Core:" + str(core_clusters) + ", Singletons: " + str(singleton_clusters)
    return stats, total_clusters_to_strains, core_clusters_to_strains, singleton_clusters_to_strains


def create_strains_clusters_map(clusters_file):
    strains_map = {}
    with open(clusters_file, 'r') as clusters_db:
        cur_cluster = None
        for line in clusters_db:
            if line.startswith(">Cluster"):
                cluster_index = int(line.split()[-1])
                cur_cluster = Cluster(cluster_index)
            else:
                match = CLUSTER_STRAIN_PATTERN.match(line)
                strain_index = int(match.group(1))
                cur_cluster.add_strain(strain_index)
                cur_strain = strains_map[strain_index] if strain_index in strains_map.keys() else Strain(strain_index)
                cur_strain.add_cluster(cur_cluster)
                strains_map[strain_index] = cur_strain
    total_strains_count = len(strains_map)
    return strains_map, total_strains_count


def get_strain_contigs(strain_genomic_file):
    contigs = 0
    for line in strain_genomic_file:
        if line.startswith(">") and "plasmid" not in line:
            contigs += 1
    return contigs


def get_strain_pseudogenes(strain_cds_file):
    pseudogenes = 0
    for line in strain_cds_file:
        if "pseudo=true" in line:
            pseudogenes += 1
    return pseudogenes


def get_genomic_stats_per_strain():
    genomic_stats = {}
    contigs_to_strains = []
    pseudogenes_to_strains = []
    for strain_dir in os.listdir(STRAINS_DIR):
        strain_dir_files = os.listdir(os.path.join(STRAINS_DIR, strain_dir))
        cds_file_name = [f for f in strain_dir_files if CDS_FROM_GENOMIC_PATTERN in f][0]
        genomic_file_name = [f for f in strain_dir_files if GENOMIC_PATTERN in f and CDS_FROM_GENOMIC_PATTERN not in f][0]
        cds_file = strain_index_file = genomic_file = None
        try:
            strain_index_file = open(os.path.join(STRAINS_DIR, strain_dir, STRAIN_INDEX_FILE))
            strain_index = int(strain_index_file.readline())
            if genomic_file_name.endswith('gz'):
                genomic_file = gzip.open(os.path.join(STRAINS_DIR, strain_dir, genomic_file_name), 'rt')
            else:
                genomic_file = open(os.path.join(STRAINS_DIR, strain_dir, genomic_file_name))
            if cds_file_name.endswith('gz'):
                cds_file = gzip.open(os.path.join(STRAINS_DIR, strain_dir, cds_file_name), 'rt')
            else:
                cds_file = open(os.path.join(STRAINS_DIR, strain_dir, cds_file_name))
            strain_contigs = get_strain_contigs(genomic_file)
            contigs_to_strains.append(strain_contigs)
            strain_pseudogenes = get_strain_pseudogenes(cds_file)
            pseudogenes_to_strains.append(strain_pseudogenes)
            genomic_stats[strain_index] = "Contigs: " + str(strain_contigs) + ", Pseudogenes: " + str(strain_pseudogenes)
        finally:
            if genomic_file is not None:
                genomic_file.close()
            if cds_file is not None:
                cds_file.close()
            if strain_index_file is not None:
                strain_index_file.close()
    return genomic_stats, contigs_to_strains, pseudogenes_to_strains
