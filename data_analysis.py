import re
from collections import defaultdict

CLUSTER_STRAIN_PATTERN = re.compile("[0-9a,> \t]+\[(\d+)\]\[(\d+)\]")


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
    for strain in strains_map.values():
        total_clusters = len(strain.containing_clusters)
        core_clusters = len(strain.get_strain_core_clusters(total_strains_count))
        singleton_clusters = len(strain.get_strain_singleton_clusters())
        stats[strain.index] = {"total_clusters": total_clusters, "core": core_clusters, "singletons": singleton_clusters}
    return stats


def create_strains_clusters_map(clusters_file):
    strains_map = {}
    with open(clusters_file, 'r') as clusters_db:
        cur_cluster = None
        for line in clusters_db.readlines():
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
