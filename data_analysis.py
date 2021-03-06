import gzip
import logging
import sys
from collections import defaultdict
import os
import pandas
from Bio import SeqIO

from constants import STRAINS_DIR, CDS_FROM_GENOMIC_PATTERN, GENOMIC_PATTERN, STRAIN_INDEX_FILE, CLUSTER_STRAIN_PATTERN, \
    CD_HIT_CLUSTERS_OUTPUT_FILE, CD_HIT_EST_CLUSTERS_OUTPUT_FILE, CLUSTER_PSEUDOGENE_PATTERN, \
    CLUSTER_2ND_STAGE_SEQ_LEN_PATTERN, CD_HIT_EST_MULTIPLE_PROTEIN_CLUSTERS_OUTPUT_FILE, COMBINED_CDS_FILE_PATH, \
    FASTA_FILE_TYPE, COMBINED_STRAIN_REPS_CDS_PATH, COMBINED_STRAIN_PSEUDOGENES_PATH, BLAST_RESULTS_FILE, \
    BLAST_PSEUDOGENE_PATTERN, COMBINED_PSEUDOGENES_WITHOUT_BLAST_HIT_PATH, CLUSTERS_NT_SEQS_DIR, \
    PROTEIN_CORE_CLUSTERS_PKL, MLST_GENES, STRAINS_COUNT, DATA_DIR
from nucleotide_preprocessor import get_strain_index

logger = logging.getLogger(__name__)


class Cluster:
    def __init__(self, index):
        self.index = index
        self.member_strains = defaultdict(int)
        self.member_strains_seqs = {}

    def add_strain(self, strain_index):
        self.member_strains[strain_index] += 1

    def add_strain_seq(self, strain_index, seq_index):
        strain_seqs = self.member_strains_seqs[strain_index] if strain_index in self.member_strains_seqs.keys() else []
        strain_seqs.append(seq_index)
        self.member_strains_seqs[strain_index] = strain_seqs

    def get_cluster_strains_num(self):
        return len(self.member_strains.keys())

    def get_cluster_proteins_num(self):
        return sum(self.member_strains.values())


class NucleotideCluster(Cluster):
    def __init__(self, index):
        super().__init__(index)
        self.member_pseudogenes = {}
        self.member_protein_seqs = {}
        self.total_protein_len = 0
        self.total_pseudogene_len = 0
        self.representative = ""

    def add_nucleotide(self, strain_index, seq_index, seq_len, is_pseudogene, is_cluster_rep):
        if is_pseudogene:
            strain_pseudogenes = self.member_pseudogenes[
                strain_index] if strain_index in self.member_pseudogenes.keys() else []
            strain_pseudogenes.append(seq_index)
            self.member_pseudogenes[strain_index] = strain_pseudogenes
            self.total_pseudogene_len += seq_len
        else:
            strain_1st_stage_reps = self.member_protein_seqs[
                strain_index] if strain_index in self.member_protein_seqs.keys() else []
            strain_1st_stage_reps.append(seq_index)
            self.member_protein_seqs[strain_index] = strain_1st_stage_reps
            self.total_protein_len += seq_len
        if is_cluster_rep:
            self.representative = "Strain idx: %s, Seq idx: %s, Pseudogene: %r" % (str(strain_index), str(seq_index), bool(is_pseudogene))

    def has_reps(self):
        return len(self.member_protein_seqs) == 0

    def get_avg_protein_seq_len(self):
        return int(round(self.total_protein_len / len(self.member_protein_seqs))) if len(self.member_protein_seqs) > 0 else 0

    def get_avg_pseudogene_seq_len(self):
        return int(round(self.total_pseudogene_len / len(self.member_pseudogenes))) if len(self.member_pseudogenes) > 0 else 0


class Strain:
    def __init__(self, index):
        self.index = index
        self.containing_clusters = {}
        self.seq_clusters = {}

    def add_cluster(self, cluster):
        self.containing_clusters[cluster.index] = cluster

    def add_seq_cluster(self, seq_index, cluster_index):
        self.seq_clusters[seq_index] = cluster_index

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
    return strains_map, clusters_map, total_strains_count, total_core_clusters


def create_1st_stage_sequences_clusters_map(clusters_file):
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
                seq_index = int(match.group(2))
                cur_cluster.add_strain(strain_index)
                cur_cluster.add_strain_seq(strain_index, seq_index)
                cur_strain = strains_map[strain_index] if strain_index in strains_map.keys() else Strain(strain_index)
                cur_strain.add_seq_cluster(seq_index, cluster_index)
                strains_map[strain_index] = cur_strain
    return strains_map, clusters_map


def create_nucleotide_clusters_map(clusters_file):
    strains_map = {}
    clusters_map = {}
    with open(clusters_file, 'r') as clusters_db:
        cur_cluster = None
        for line in clusters_db:
            if line.startswith(">Cluster"):
                cluster_index = int(line.split()[-1])
                cur_cluster = NucleotideCluster(cluster_index)
                clusters_map[cluster_index] = cur_cluster
            else:
                is_pseudogene = CLUSTER_PSEUDOGENE_PATTERN.search(line)
                is_representative = '*' in line
                index_match = CLUSTER_STRAIN_PATTERN.match(line)
                len_match = CLUSTER_2ND_STAGE_SEQ_LEN_PATTERN.search(line)
                if not index_match:
                    raise ValueError("line in clusters file %s does not match strain index pattern %s" % (line, CLUSTER_STRAIN_PATTERN.pattern))
                if not len_match:
                    raise ValueError("line in clusters file %s does not match seq len pattern %s" % (line, CLUSTER_2ND_STAGE_SEQ_LEN_PATTERN.pattern))
                strain_index = int(index_match.group(1))
                seq_index = int(index_match.group(2))
                seq_len = int(len_match.group(1))
                cur_cluster.add_strain(strain_index)
                cur_cluster.add_nucleotide(strain_index, seq_index, seq_len, is_pseudogene, is_representative)
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
    genes = 0
    pseudogenes = 0
    for line in strain_cds_file:
        if line.startswith('>'):
            if "pseudo=true" in line:
                pseudogenes += 1
            else:
                genes += 1
    return genes, pseudogenes


def get_1st_stage_stats_per_strain():
    strains_map, _, total_strains_count, total_core_clusters = create_strains_clusters_map(CD_HIT_CLUSTERS_OUTPUT_FILE)
    df = pandas.DataFrame(index=range(total_strains_count), columns=('strain_name', 'total_clusters', 'core_clusters', 'missing_core', 'singletons', 'contigs', 'pseudogenes', 'genes'))
    for strain in strains_map.values():
        total_clusters = len(strain.containing_clusters)
        core_clusters = len(strain.get_strain_core_clusters(total_strains_count))
        missing_core = 100 - (core_clusters / total_core_clusters * 100)
        singletons = len(strain.get_strain_singleton_clusters())
        df.loc[strain.index] = ['', total_clusters, core_clusters, missing_core, singletons, 0, 0, 0]
    for strain_dir in os.listdir(STRAINS_DIR):
        strain_dir_files = os.listdir(os.path.join(STRAINS_DIR, strain_dir))
        cds_file_name = [f for f in strain_dir_files if CDS_FROM_GENOMIC_PATTERN in f][0]
        genomic_file_name = [f for f in strain_dir_files if GENOMIC_PATTERN in f and CDS_FROM_GENOMIC_PATTERN not in f][0]
        cds_file = strain_index_file = genomic_file = None
        try:
            strain_index_file = open(os.path.join(STRAINS_DIR, strain_dir, STRAIN_INDEX_FILE))
            strain_index = int(strain_index_file.readline())
            df.loc[strain_index]['strain_name'] = strain_dir
            if genomic_file_name.endswith('gz'):
                genomic_file = gzip.open(os.path.join(STRAINS_DIR, strain_dir, genomic_file_name), 'rt')
            else:
                genomic_file = open(os.path.join(STRAINS_DIR, strain_dir, genomic_file_name))
            if cds_file_name.endswith('gz'):
                cds_file = gzip.open(os.path.join(STRAINS_DIR, strain_dir, cds_file_name), 'rt')
            else:
                cds_file = open(os.path.join(STRAINS_DIR, strain_dir, cds_file_name))
            strain_contigs = get_strain_contigs(genomic_file)
            strain_genes, strain_pseudogenes = get_strain_pseudogenes(cds_file)
            df.loc[strain_index]['contigs'] = strain_contigs
            df.loc[strain_index]['genes'] = strain_genes
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
    first_stage_strain_seq_cluster_map, first_stage_clusters_map = create_1st_stage_sequences_clusters_map(CD_HIT_CLUSTERS_OUTPUT_FILE)
    logger.info("Creating 2nd stage strains & clusters maps from CD-HIT-EST output")
    second_stage_strains_map, second_stage_clusters_map = create_nucleotide_clusters_map(CD_HIT_EST_CLUSTERS_OUTPUT_FILE)
    total_strains_count = len(second_stage_strains_map.keys())
    total_clusters_count = len(second_stage_clusters_map.keys())

    strains_df = pandas.DataFrame(index=range(total_strains_count), columns=('total_pseudogenes', 'pseudogenes_in_clusters_without_reps'))
    for strain in second_stage_strains_map.values():
        strains_df.loc[strain.index]['total_pseudogenes'] = first_stage_data.loc[strain.index]['pseudogenes']
        strains_df.loc[strain.index]['pseudogenes_in_clusters_without_reps'] = 0
        for cluster in second_stage_clusters_map.values():
            if not cluster.has_reps() and strain.index in cluster.member_pseudogenes.keys():
                strains_df.loc[strain.index]['pseudogenes_in_clusters_without_reps'] += len(cluster.member_pseudogenes[strain.index])

    clusters_df = pandas.DataFrame(index=range(total_clusters_count), columns=('total_strains', '1st_stage_reps', 'strains_in_rep_1st_stage_cluster'))
    for cluster in second_stage_clusters_map.values():
        clusters_df.loc[cluster.index]['total_strains'] = cluster.get_cluster_strains_num()
        clusters_df.loc[cluster.index]['1st_stage_reps'] = len(cluster.member_protein_seqs)
        if len(cluster.member_protein_seqs) == 1:
            [(representative_strain_index, representative_seq_index)] = cluster.member_protein_seqs.items()
            representative_cluster_id = first_stage_strain_seq_cluster_map[representative_strain_index].seq_clusters[representative_seq_index[0]]
            representative_cluster = first_stage_clusters_map[representative_cluster_id]
            clusters_df.loc[cluster.index]['strains_in_rep_1st_stage_cluster'] = representative_cluster.get_cluster_strains_num()
    return strains_df, clusters_df


def get_1st_stage_strains_per_clusters_stats():
    logger.info("Creating 1st stage clusters map from CD-HIT output")
    _, first_stage_clusters_map, total_strains_count, _ = create_strains_clusters_map(CD_HIT_CLUSTERS_OUTPUT_FILE)
    strains_percentage_per_cluster = []
    for cluster in first_stage_clusters_map.values():
        strains_percentage_per_cluster.append((cluster.get_cluster_strains_num() / total_strains_count) * 100)
    return strains_percentage_per_cluster


def get_2nd_stage_stats_per_cluster():
    logger.info("Creating 2nd stage strains & clusters maps from CD-HIT-EST output")
    type1 = type2 = type3 = type4 = 0
    first_stage_strain_seq_cluster_map, first_stage_clusters_map = create_1st_stage_sequences_clusters_map(CD_HIT_CLUSTERS_OUTPUT_FILE)
    second_stage_strains_map, second_stage_clusters_map = create_nucleotide_clusters_map(CD_HIT_EST_CLUSTERS_OUTPUT_FILE)
    total_clusters_count = len(second_stage_clusters_map.keys())
    clusters_df = pandas.DataFrame(index=range(total_clusters_count), columns=('protein_seqs', 'strains_of_protein_seqs',
                                                                               'pseudogenes', 'strains_of_pseudogenes',
                                                                               'avg_protein_len', 'avg_pseudogene_len',
                                                                               'cluster_type', 'representative',
                                                                               'strains_in_rep_1st_stage_cluster',
                                                                               'proteins_in_rep_1st_stage_cluster'))
    for cluster in second_stage_clusters_map.values():
        clusters_df.loc[cluster.index]['protein_seqs'] = len([seq for strain_proteins in cluster.member_protein_seqs.values() for seq in strain_proteins])
        clusters_df.loc[cluster.index]['strains_of_protein_seqs'] = len(cluster.member_protein_seqs.keys())
        clusters_df.loc[cluster.index]['pseudogenes'] = len([seq for strain_pseudos in cluster.member_pseudogenes.values() for seq in strain_pseudos])
        clusters_df.loc[cluster.index]['strains_of_pseudogenes'] = len(cluster.member_pseudogenes.keys())
        clusters_df.loc[cluster.index]['avg_protein_len'] = cluster.get_avg_protein_seq_len()
        clusters_df.loc[cluster.index]['avg_pseudogene_len'] = cluster.get_avg_pseudogene_seq_len()
        if clusters_df.loc[cluster.index]['protein_seqs'] == 1:
            if clusters_df.loc[cluster.index]['pseudogenes'] > 0:
                cluster_type = 1
                type1 += 1
            else:
                cluster_type = 3
                type3 += 1
            [(representative_strain_index, representative_seq_index)] = cluster.member_protein_seqs.items()
            representative_cluster_id = first_stage_strain_seq_cluster_map[representative_strain_index].seq_clusters[representative_seq_index[0]]
            representative_cluster = first_stage_clusters_map[representative_cluster_id]
            clusters_df.loc[cluster.index]['strains_in_rep_1st_stage_cluster'] = representative_cluster.get_cluster_strains_num()
            clusters_df.loc[cluster.index]['proteins_in_rep_1st_stage_cluster'] = representative_cluster.get_cluster_proteins_num()
        elif clusters_df.loc[cluster.index]['protein_seqs'] > 1:
            cluster_type = 4
            type4 += 1
        else:
            cluster_type = 2
            type2 += 1
        clusters_df.loc[cluster.index]['cluster_type'] = cluster_type
        clusters_df.loc[cluster.index]['representative'] = cluster.representative
    logger.info("Cluster type amounts:")
    logger.info("Type 1 (1 protein + pseudogenes): %d" % type1)
    logger.info("Type 2 (0 proteins + pseudogenes): %d" % type2)
    logger.info("Type 3 (1 protein + 0 pseudogenes): %d" % type3)
    logger.info("Type 4 (multiple proteins +- pseudogenes): %d" % type4)

    return clusters_df


def filter_2nd_stage_clusters_with_multiple_proteins():
    with open(CD_HIT_EST_CLUSTERS_OUTPUT_FILE, 'r') as clusters_db:
        with open(CD_HIT_EST_MULTIPLE_PROTEIN_CLUSTERS_OUTPUT_FILE, 'w') as output:
            cluster = ''
            num_of_proteins = 0
            for line in clusters_db:
                if line.startswith(">Cluster"):
                    if num_of_proteins > 1:
                        output.write(cluster)
                    cluster = ''
                    num_of_proteins = 0
                    cluster += line
                else:
                    cluster += line
                    if not CLUSTER_PSEUDOGENE_PATTERN.search(line):
                        num_of_proteins += 1
            if num_of_proteins > 1:
                output.write(cluster)
            return


def split_2nd_stage_combined_fasta_to_reps_pseudogenes():
    with open(COMBINED_STRAIN_REPS_CDS_PATH, "w") as reps_cds_file:
        with open(COMBINED_STRAIN_PSEUDOGENES_PATH, "w") as pseudogenes_file:
            with open(COMBINED_CDS_FILE_PATH) as cds_file:
                strain_cds_seq_iter = SeqIO.parse(cds_file, FASTA_FILE_TYPE)
                for strain_cds_seq in strain_cds_seq_iter:
                    if "pseudo=true" in strain_cds_seq.description:
                        SeqIO.write(strain_cds_seq, pseudogenes_file, FASTA_FILE_TYPE)
                    else:
                        SeqIO.write(strain_cds_seq, reps_cds_file, FASTA_FILE_TYPE)


def get_pseudogenes_from_blast_results():
    pseudogenes = {}
    with open(BLAST_RESULTS_FILE) as blast_result:
        prefix = "bla"
        for line in blast_result:
            if not line.startswith(prefix):
                pseudogene_prefix = BLAST_PSEUDOGENE_PATTERN.match(line.lstrip())
                prefix = pseudogene_prefix.group()
                strain_idx = pseudogene_prefix.group(1)
                strain_seqs = pseudogenes[strain_idx] if strain_idx in pseudogenes.keys() else []
                seq_idx = pseudogene_prefix.group(2)
                if seq_idx not in strain_seqs:
                    strain_seqs.append(seq_idx)
                    pseudogenes[strain_idx] = strain_seqs
    print("found %s pseudogenes in blast results" % len(pseudogenes))
    return pseudogenes


def get_pseudogenes_without_blast_hits_fasta():
    pseudogenes_with_hits = get_pseudogenes_from_blast_results()
    all_pseudogenes_iter = SeqIO.parse(open(COMBINED_STRAIN_PSEUDOGENES_PATH), FASTA_FILE_TYPE)
    with open(COMBINED_PSEUDOGENES_WITHOUT_BLAST_HIT_PATH, "w") as pseudogenes_file:
        for seq in all_pseudogenes_iter:
            pseudogene_prefix = BLAST_PSEUDOGENE_PATTERN.match(seq.description.lstrip())
            strain_idx = pseudogene_prefix.group(1)
            seq_idx = pseudogene_prefix.group(2)
            strain_seqs = pseudogenes_with_hits[strain_idx] if strain_idx in pseudogenes_with_hits.keys() else None
            if strain_seqs is None or seq_idx not in strain_seqs:
                SeqIO.write(seq, pseudogenes_file, FASTA_FILE_TYPE)


def get_core_clusters():
    first_stage_strain_seq_cluster_map, first_stage_clusters_map = create_1st_stage_sequences_clusters_map(
        CD_HIT_CLUSTERS_OUTPUT_FILE)
    total_strains_count = len(first_stage_strain_seq_cluster_map)
    core_clusters_multiple_strain_seqs = {}
    clusters_to_remove = []
    core_clusters = {ind: cluster for (ind, cluster) in first_stage_clusters_map.items() if cluster.get_cluster_strains_num() / total_strains_count >= 0.9}
    for cluster in core_clusters.values():
        if not all(i == 1 for i in cluster.member_strains.values()):
            core_clusters_multiple_strain_seqs[cluster.index] = core_clusters[cluster.index]
    for index in core_clusters_multiple_strain_seqs.keys():
        core_clusters.pop(index)
    for cluster in core_clusters_multiple_strain_seqs.values():
        strains_to_remove = [s for s in cluster.member_strains.keys() if cluster.member_strains[s] > 1]
        for strain in strains_to_remove:
            cluster.member_strains.pop(strain)
        if cluster.get_cluster_strains_num() / total_strains_count < 0.9:
            clusters_to_remove.append(cluster.index)
    for cluster_index in clusters_to_remove:
        core_clusters_multiple_strain_seqs.pop(cluster_index)
    return core_clusters, core_clusters_multiple_strain_seqs


def export_protein_clusters_to_nucleotide_fasta_files():
    import pickle
    if os.path.exists(PROTEIN_CORE_CLUSTERS_PKL):
        logger.info("Loading protein core clusters from pickle")
        with open(PROTEIN_CORE_CLUSTERS_PKL, "rb") as f:
            core_clusters = pickle.load(f)
    else:
        logger.info("Generating protein core clusters")
        core_clusters, _ = get_core_clusters()
        with open(PROTEIN_CORE_CLUSTERS_PKL, "wb") as f:
            pickle.dump(core_clusters, f)
    if not os.path.exists(CLUSTERS_NT_SEQS_DIR):
        os.makedirs(CLUSTERS_NT_SEQS_DIR)

    for strain_index, strain_dir, cds_file in iterate_strains_cds():
        strain_cds = list(SeqIO.parse(cds_file, FASTA_FILE_TYPE))
        logger.info("Parsing strain %s index %d" % (str(strain_dir), strain_index))
        for cluster in core_clusters.values():
            if strain_index in cluster.member_strains_seqs.keys():
                strain_protein_seq_index = cluster.member_strains_seqs[strain_index][0]
                strain_protein_seq = strain_cds[strain_protein_seq_index - 1]
                strain_protein_seq.description = "[" + str(strain_index) + "][" + str(strain_protein_seq_index) + "]" + strain_protein_seq.description
                cluster.member_strains_seqs[strain_index] = strain_protein_seq

        if cds_file is not None:
            cds_file.close()

    for cluster in core_clusters.values():
        logger.info("Saving cluster %d to file" % cluster.index)
        cluster_protein_nts_file_path = os.path.join(CLUSTERS_NT_SEQS_DIR, "cluster_" + str(cluster.index))
        cluster_protein_nts_file = open(cluster_protein_nts_file_path, "a+")
        cluster_seqs = [seq for idx, seq in sorted(cluster.member_strains_seqs.items())]
        SeqIO.write(cluster_seqs, cluster_protein_nts_file, FASTA_FILE_TYPE)
        cluster_protein_nts_file.close()


def shorten_seq_names_in_clusters():
    for cluster_file in os.listdir(CLUSTERS_NT_SEQS_DIR):
        cluster_file_short_seq_names = cluster_file + "_short_names"
        logger.info("Shortening seq names for %s" % cluster_file)
        if not os.path.exists(os.path.join(CLUSTERS_NT_SEQS_DIR, cluster_file_short_seq_names)):
            with open(os.path.join(CLUSTERS_NT_SEQS_DIR, cluster_file), "r") as f1:
                cluster_cds = list(SeqIO.parse(f1, FASTA_FILE_TYPE))
                for cds in cluster_cds:
                        cds.id = cds.description.split(' ')[1]
                        cds.description = ''
                with open(os.path.join(CLUSTERS_NT_SEQS_DIR, cluster_file_short_seq_names), "w") as f2:
                    SeqIO.write(cluster_cds, f2, FASTA_FILE_TYPE)


def build_strain_names_map():
    strain_names_map = {}
    downloaded_strains = os.listdir(STRAINS_DIR)
    for strain_dir in downloaded_strains:
        strain_index = get_strain_index(strain_dir)
        strain_names_map[strain_index] = strain_dir
    return strain_names_map


def iterate_strains_cds():
    """Iterate over the CDS file for each downloaded strain and return the strain index and cds file handle"""
    downloaded_strains = os.listdir(STRAINS_DIR)
    for strain_dir in downloaded_strains:
        strain_index = get_strain_index(strain_dir)
        strain_dir_files = os.listdir(os.path.join(STRAINS_DIR, strain_dir))
        cds_file_name = [f for f in strain_dir_files if CDS_FROM_GENOMIC_PATTERN in f][0]
        if cds_file_name is None:
            raise RuntimeError("Failed to find a cds file for strain %s" % str(strain_dir))
        if cds_file_name.endswith('gz'):
            cds_file = gzip.open(os.path.join(STRAINS_DIR, strain_dir, cds_file_name), 'rt')
        else:
            cds_file = open(os.path.join(STRAINS_DIR, strain_dir, cds_file_name))
        yield strain_index, strain_dir, cds_file


def get_strains_mlst_genes():
    mlst_variants_records = {}
    for gene in MLST_GENES:
        mlst_variants_records[gene] = list(SeqIO.parse(open(gene + ".fas"), FASTA_FILE_TYPE))
    strains_mlst_vectors = pandas.DataFrame(index=range(STRAINS_COUNT), columns=MLST_GENES)
    for strain_index, strain_dir, cds_file in iterate_strains_cds():
        strain_cds = SeqIO.parse(cds_file, FASTA_FILE_TYPE)
        genes_found = 0
        for strain_gene in strain_cds:
            if genes_found == len(MLST_GENES):
                break
            for mlst_gene in MLST_GENES:
                if mlst_gene in strain_gene.description:
                    logger.info("found mlst gene " + mlst_gene + " in strain " + strain_gene.id + ", idx " + strain_index)
                    genes_found += 1
                    gene_variants = mlst_variants_records[mlst_gene]
                    for variant in gene_variants:
                        if str(variant.seq) in str(strain_gene.seq):
                            logger.info("found match: gene " + mlst_gene + ", variant " + variant.id + ", seq: " + variant.seq + ""
                                         "contained in strain gene" + strain_gene.id + ", seq: " + strain_gene.seq)
                            variant_id = variant.id[5:]
                            strains_mlst_vectors.loc[strain_index][mlst_gene] = variant_id
                            break
                    break
        cds_file.close()
    return strains_mlst_vectors

#TODO compare each strain vector to allelic profiles table - pandas.read_table(MLST_ALLELIC_PROFILE_PATH)
#TODO Build tree diagram and paint with mlst value per strain




