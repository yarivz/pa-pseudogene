import gzip
import logging
from collections import defaultdict
import os
import pandas
from constants import STRAINS_DIR, CDS_FROM_GENOMIC_PATTERN, GENOMIC_PATTERN, STRAIN_INDEX_FILE, CLUSTER_STRAIN_PATTERN, \
    CD_HIT_CLUSTERS_OUTPUT_FILE, CLUSTER_1ST_STAGE_REPRESENTATIVE_PATTERN, \
    CD_HIT_EST_CLUSTER_REPS_OUTPUT_FILE, CD_HIT_EST_CLUSTERS_OUTPUT_FILE, CLUSTER_PSEUDOGENE_PATTERN

logger = logging.getLogger(__name__)


class Cluster:
    def __init__(self, index):
        self.index = index
        self.member_strains = defaultdict(int)

    def add_strain(self, strain_index):
        self.member_strains[strain_index] += 1

    def get_cluster_strains_num(self):
        return len(self.member_strains)


class NucleotideCluster(Cluster):
    def __init__(self, index):
        super().__init__(index)
        self.member_pseudogenes = defaultdict(int)
        self.member_1st_stage_reps = defaultdict(int)

    def add_nucleotide(self, strain_index, seq_index, is_pseudogene):
        if is_pseudogene:
            strain_pseudogenes = self.member_pseudogenes[
                strain_index] if strain_index in self.member_pseudogenes.keys() else []
            strain_pseudogenes.append(seq_index)
            self.member_pseudogenes[strain_index] = strain_pseudogenes
        else:
            strain_1st_stage_reps = self.member_1st_stage_reps[
                strain_index] if strain_index in self.member_1st_stage_reps.keys() else []
            strain_1st_stage_reps.append(seq_index)
            self.member_1st_stage_reps[strain_index] = strain_1st_stage_reps

    def has_reps(self):
        return len(self.member_1st_stage_reps) == 0


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

    def get_strain_pseudogenes_in_clusters_without_rep(self):
        clusters_without_rep = [c for c in self.containing_clusters.values() if not c.member_1st_stage_reps]
        strain_pseudogenes = []
        for c in clusters_without_rep:
            if c.member_pseudogenes.get(self.index):
                strain_pseudogenes.extend(c.member_pseudogenes.get(self.index))
        return strain_pseudogenes


def create_strains_clusters_map(clusters_file):
    strains_map = {}
    clusters_map = {}
    with open(clusters_file, 'r') as clusters_db:
        cur_cluster = None
        for line in clusters_db:
            if line.startswith(">Cluster"):
                cluster_index = int(line.split()[-1])
                cur_cluster = Cluster(cluster_index)
                clusters_map[cluster_index] = cur_cluster
            else:
                match = CLUSTER_STRAIN_PATTERN.match(line)
                strain_index = int(match.group(1))
                cur_cluster.add_strain(strain_index)
                cur_strain = strains_map[strain_index] if strain_index in strains_map.keys() else Strain(strain_index)
                cur_strain.add_cluster(cur_cluster)
                strains_map[strain_index] = cur_strain
    total_strains_count = len(strains_map)
    total_core_clusters = len([c for c in clusters_map.values() if c.get_cluster_strains_num() / total_strains_count >= 0.9])
    return strains_map, total_strains_count, total_core_clusters


def create_clusters_map(clusters_file):
    clusters_map = {}
    with open(clusters_file, 'r') as clusters_db:
        cur_cluster = None
        for line in clusters_db:
            if line.startswith(">Cluster"):
                cluster_index = int(line.split()[-1])
                cur_cluster = Cluster(cluster_index)
                clusters_map[cluster_index] = cur_cluster
            else:
                match = CLUSTER_STRAIN_PATTERN.match(line)
                strain_index = int(match.group(1))
                cur_cluster.add_strain(strain_index)
    return clusters_map


def create_1st_stage_representative_clusters_map(representatives_file):
    """
    :param representatives_file:
    :return: dict(strain index : dict(seq index : cluster index))
    """
    representative_map = {}
    with open(representatives_file, 'r') as reps_db:
        for line in reps_db:
            if line.startswith(">") and CLUSTER_1ST_STAGE_REPRESENTATIVE_PATTERN.search(line):
                match = CLUSTER_1ST_STAGE_REPRESENTATIVE_PATTERN.match(line)
                strain_index = int(match.group(1))
                seq_index = int(match.group(2))
                cluster_index = int(match.group(3))
                cur_strain = representative_map[strain_index] if strain_index in representative_map.keys() else {}
                cur_strain[seq_index] = cluster_index
                representative_map[strain_index] = cur_strain
    return representative_map


def create_nucleotide_clusters_map(clusters_file):
    strains_map = {}
    clusters_map = {}
    with open(clusters_file, 'r') as clusters_db:
        cur_cluster = None
        for line in clusters_db:
            print(line)
            if line.startswith(">Cluster"):
                cluster_index = int(line.split()[-1])
                cur_cluster = NucleotideCluster(cluster_index)
                clusters_map[cluster_index] = cur_cluster
            else:
                print(CLUSTER_PSEUDOGENE_PATTERN)
                print(CLUSTER_STRAIN_PATTERN)
                is_pseudogene = CLUSTER_PSEUDOGENE_PATTERN.search(line)
                match = CLUSTER_STRAIN_PATTERN.match(line)
                strain_index = int(match.group(1))
                seq_index = int(match.group(2))
                cur_cluster.add_strain(strain_index)
                cur_cluster.add_nucleotide(strain_index, seq_index, is_pseudogene)
                cur_strain = strains_map[strain_index] if strain_index in strains_map.keys() else Strain(strain_index)
                cur_strain.add_cluster(cur_cluster)
                strains_map[strain_index] = cur_strain
    return strains_map, clusters_map


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


def get_1st_stage_stats_per_strain():
    strains_map, total_strains_count, total_core_clusters = create_strains_clusters_map(CD_HIT_CLUSTERS_OUTPUT_FILE)
    df = pandas.DataFrame(index=range(total_strains_count), columns=('total_clusters', 'core_clusters', 'missing_core', 'singletons', 'contigs', 'pseudogenes'))
    for strain in strains_map.values():
        total_clusters = len(strain.containing_clusters)
        core_clusters = len(strain.get_strain_core_clusters(total_strains_count))
        missing_core = 100 - (core_clusters / total_core_clusters * 100)
        singletons = len(strain.get_strain_singleton_clusters())
        df.loc[strain.index] = [total_clusters, core_clusters, missing_core, singletons, 0, 0]
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
            strain_pseudogenes = get_strain_pseudogenes(cds_file)
            df.loc[strain_index]['contigs'] = strain_contigs
            df.loc[strain_index]['pseudogenes'] = strain_pseudogenes
        finally:
            if genomic_file is not None:
                genomic_file.close()
            if cds_file is not None:
                cds_file.close()
            if strain_index_file is not None:
                strain_index_file.close()
    return df


def get_2nd_stage_stats_per_strain(first_stage_data):
    logger.info("Creating 1st stage clusters map from CD-HIT output")
    first_stage_clusters_map = create_clusters_map(CD_HIT_CLUSTERS_OUTPUT_FILE)
    logger.info("Creating 1st stage reps to cluster id map from CD-HIT-EST output")
    first_stage_reps_to_cluster_id_map = create_1st_stage_representative_clusters_map(CD_HIT_EST_CLUSTER_REPS_OUTPUT_FILE)
    logger.info("Creating 2nd stage strains & clusters maps from CD-HIT-EST output")
    second_stage_strains_map, second_stage_clusters_map = create_nucleotide_clusters_map(CD_HIT_EST_CLUSTERS_OUTPUT_FILE)
    total_strains_count = len(second_stage_strains_map.keys())
    total_clusters_count = len(second_stage_clusters_map.keys())

    strains_df = pandas.DataFrame(index=range(total_strains_count), columns=('totatl_pseudogenes', 'pseudogenes_in_cluster_without_reps'))
    for strain in second_stage_strains_map.values():
        strains_df.loc[strain.index]['totatl_pseudogenes'] = len(first_stage_data.loc[strain.index]['pseudogenes'])
        strains_df.loc[strain.index]['pseudogenes_in_cluster_without_reps'] = 0
        for cluster in second_stage_clusters_map.values():
            if not cluster.has_reps:
                strains_df.loc[strain.index]['pseudogenes_in_cluster_without_reps'] += len(cluster.member_pseudogenes)

    clusters_df = pandas.DataFrame(index=range(total_clusters_count), columns=('total_strains', '1st_stage_reps', 'strains_in_rep_1st_stage_cluster'))
    for cluster in second_stage_clusters_map.values():
        clusters_df.loc[cluster.index]['total_strains'] = cluster.get_cluster_strains_num()
        clusters_df.loc[cluster.index]['1st_stage_reps'] = cluster.member_1st_stage_reps
        if len(cluster.member_1st_stage_reps) == 1:
            representative_strain_seq = cluster.member_1st_stage_reps.items
            representative_strain_index = representative_strain_seq[0]
            representative_seq_index = representative_strain_seq[1]
            representative_cluster_id = first_stage_reps_to_cluster_id_map[representative_strain_index][representative_seq_index]
            representative_cluster = first_stage_clusters_map[representative_cluster_id]
            clusters_df.loc[cluster.index]['strains_in_rep_1st_stage_cluster'] = representative_cluster.get_cluster_strains_num()
    return strains_df, clusters_df













