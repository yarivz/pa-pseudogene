import logging
import os

import matplotlib

from constants import FIRST_STAGE_GRAPHS_DIR, SECOND_STAGE_GRAPHS_DIR
from data_analysis import get_1st_stage_strains_per_clusters_stats

matplotlib.use('Agg')
import matplotlib.pyplot as plt

logger = logging.getLogger(__name__)
width = 0.35


def create_1st_stage_charts(stats_df):
    if not os.path.exists(FIRST_STAGE_GRAPHS_DIR):
        os.mkdir(FIRST_STAGE_GRAPHS_DIR)
    os.chdir(FIRST_STAGE_GRAPHS_DIR)

    logger.info("Plotting clusters per strain")
    set_labels_font_size()
    sorted_stats = stats_df.sort_values('total_clusters', ascending=False).reset_index(drop=True)
    plt.bar(sorted_stats.index.values, sorted_stats['total_clusters'])
    plt.xlabel("Strains #")
    plt.ylabel("Clusters #")
    plt.title("Clusters per strain")
    plt.savefig('clusters_per_strain.pdf', format="pdf")
    plt.close()

    logger.info("Plotting core clusters per strain")
    set_labels_font_size()
    sorted_stats = stats_df.sort_values('core_clusters', ascending=False).reset_index(drop=True)
    plt.bar(sorted_stats.index.values, sorted_stats['core_clusters'])
    plt.xlabel("Strains #")
    plt.ylabel("Core Clusters #")
    plt.title("Core Clusters per strain")
    plt.savefig('core_clusters_per_strain.pdf', format="pdf")
    plt.close()

    logger.info("Plotting missing core clusters per strain")
    set_labels_font_size()
    sorted_stats = stats_df.sort_values('missing_core', ascending=False).reset_index(drop=True)
    plt.bar(sorted_stats.index.values, sorted_stats['missing_core'])
    plt.xlabel("Strains #")
    plt.ylabel("Missing Core %")
    plt.title("Missing Core % per strain")
    plt.savefig('missing_core_clusters_per_strain.pdf', format="pdf")
    plt.close()

    logger.info("Plotting singleton clusters per strain")
    set_labels_font_size()
    sorted_stats = stats_df.sort_values('singletons', ascending=False).reset_index(drop=True)
    plt.bar(sorted_stats.index.values, sorted_stats['singletons'])
    plt.xlabel("Strains #")
    plt.ylabel("Singletons #")
    plt.title("Singletons per strain")
    plt.savefig('singleton_clusters_per_strain.pdf', format="pdf")
    plt.close()

    logger.info("Plotting pseudogenes per strain")
    set_labels_font_size()
    sorted_stats = stats_df.sort_values('pseudogenes', ascending=False).reset_index(drop=True)
    plt.bar(sorted_stats.index.values, sorted_stats['pseudogenes'])
    plt.xlabel("Strains #")
    plt.ylabel("Pseudogenes #")
    plt.title("Pseudogenes per strain")
    plt.savefig('pseudogenes_per_strain.pdf', format="pdf")
    plt.close()

    logger.info("Plotting contigs per strain")
    set_labels_font_size()
    sorted_stats = stats_df.sort_values('contigs', ascending=False).reset_index(drop=True)
    plt.bar(sorted_stats.index.values, sorted_stats['contigs'])
    plt.xlabel("Strains #")
    plt.ylabel("Contigs #")
    plt.title("Contigs per strain")
    plt.savefig('contigs_per_strain.pdf', format="pdf")
    plt.close()

    logger.info("Plotting contigs VS singletons per strain")
    set_labels_font_size()
    sorted_stats = stats_df.sort_values('contigs', ascending=False).reset_index(drop=True)
    chart1 = plt.bar(sorted_stats.index.values, sorted_stats['contigs'], width)
    chart2 = plt.bar(sorted_stats.index.values + width, sorted_stats['singletons'], width)
    plt.xlabel("Strains #")
    plt.ylabel("Contigs # / Singletons #")
    plt.title("Contigs VS Singletons per strain")
    plt.legend((chart1[0], chart2[0]), ("Contigs #", "Singletons #"))
    plt.savefig('contigs_vs_singletons_per_strain.pdf', format="pdf")
    plt.close()

    logger.info("Plotting contigs VS missing core % per strain")
    set_labels_font_size()
    sorted_stats = stats_df.sort_values('contigs', ascending=False).reset_index(drop=True)
    chart1 = plt.bar(sorted_stats.index.values, sorted_stats['contigs'], width)
    chart2 = plt.bar(sorted_stats.index.values + width, sorted_stats['missing_core'], width)
    plt.xlabel("Strains #")
    plt.ylabel("Contigs # / Missing Core %")
    plt.title("Contigs VS Missing Core % per strain")
    plt.legend((chart1[0], chart2[0]), ("Contigs #", "Missing Core %"))
    plt.savefig('contigs_vs_missing_core_per_strain.pdf', format="pdf")
    plt.close()

    logger.info("Plotting contigs VS pseudogenes per strain")
    set_labels_font_size()
    sorted_stats = stats_df.sort_values('contigs', ascending=False).reset_index(drop=True)
    chart1 = plt.bar(sorted_stats.index.values, sorted_stats['pseudogenes'], width)
    chart2 = plt.bar(sorted_stats.index.values + width, sorted_stats['contigs'], width)
    plt.xlabel("Strains #")
    plt.ylabel("Contigs # / Pseudogenes #")
    plt.title("Contigs VS Pseudogenes per strain")
    plt.legend((chart1[0], chart2[0]), ("Pseudogenes #", "Contigs #"))
    plt.savefig('contigs_vs_pseudogenes_per_strain.pdf', format="pdf")
    plt.close()

    logger.info("Plotting singletons VS missing core % per strain")
    set_labels_font_size()
    sorted_stats = stats_df.sort_values('singletons', ascending=False).reset_index(drop=True)
    chart1 = plt.bar(sorted_stats.index.values, sorted_stats['singletons'], width)
    chart2 = plt.bar(sorted_stats.index.values + width, sorted_stats['missing_core'], width)
    plt.xlabel("Strains #")
    plt.ylabel("Singletons # / Missing Core %")
    plt.title("Singletons VS Missing core % per strain")
    plt.legend((chart1[0], chart2[0]), ("Singletons #", "Missing Core %"))
    plt.savefig('singletons_vs_missing_core_per_strain.pdf', format="pdf")
    plt.close()

    logger.info("Plotting pseudogenes VS missing core % per strain")
    set_labels_font_size()
    sorted_stats = stats_df.sort_values('pseudogenes', ascending=False).reset_index(drop=True)
    chart1 = plt.bar(sorted_stats.index.values, sorted_stats['pseudogenes'], width)
    chart2 = plt.bar(sorted_stats.index.values + width, sorted_stats['missing_core'], width)
    plt.xlabel("Strains #")
    plt.ylabel("Pseudogenes # / Missing Core %")
    plt.title("strains to pseudogenes VS missing core bar chart")
    plt.legend((chart1[0], chart2[0]), ("Pseudogenes #", "Missing Core %"))
    plt.savefig('pseudogenes_vs_missing_core_per_strain.pdf', format="pdf")
    plt.close()

    logger.info("Plotting pseudogenes VS singletons per strain")
    set_labels_font_size()
    sorted_stats = stats_df.sort_values('pseudogenes', ascending=False).reset_index(drop=True)
    chart1 = plt.bar(sorted_stats.index.values, sorted_stats['pseudogenes'], width)
    chart2 = plt.bar(sorted_stats.index.values + width, sorted_stats['singletons'], width)
    plt.xlabel("Strains")
    plt.ylabel("Pseudogenes / Singletons")
    plt.title("strains to pseudogenes VS singletons bar chart")
    plt.legend((chart1[0], chart2[0]), ("Pseudogenes", "Singletons"))
    plt.savefig('pseudogenes_vs_singletons_per_strain.pdf', format="pdf")
    plt.close()

    logger.info("Plotting % of total strains per number of clusters")
    set_labels_font_size()
    data = get_1st_stage_strains_per_clusters_stats()
    plt.hist(data)
    plt.xlabel("Clusters #")
    plt.ylabel("% of total strains")
    plt.title("% of total strains per clusters histogram")
    plt.savefig('percentage_of_total_strains_per_clusters_hist.pdf', format="pdf")
    plt.close()


def create_2nd_stage_charts(strains_df, clusters_df):
    if not os.path.exists(SECOND_STAGE_GRAPHS_DIR):
        os.mkdir(SECOND_STAGE_GRAPHS_DIR)
    os.chdir(SECOND_STAGE_GRAPHS_DIR)

    logger.info("Plotting pseudogenes VS pseudogenes in clusters without reps per strain")
    set_labels_font_size()
    sorted_stats = strains_df.sort_values('total_pseudogenes', ascending=False).reset_index(drop=True)
    chart1 = plt.bar(sorted_stats.index.values, sorted_stats['total_pseudogenes'], width)
    chart2 = plt.bar(sorted_stats.index.values + width, sorted_stats['pseudogenes_in_clusters_without_reps'], width)
    plt.xlabel("Strains")
    plt.ylabel("Total Pseudogenes /\nPseudogenes in clusters without protein representatives")
    plt.title("Total strain pseudogenes VS Strain pseudogenes in clusters\nwithout protein representatives per strain")
    plt.legend((chart1[0], chart2[0]), ("Total Strain Pseudogenes", "Strain Pseudogenes in clusters\nwithout protein representatives"))
    plt.savefig('pseudogenes_vs_pseudogenes_in_repless_clusters_per_strain.pdf', format="pdf")
    plt.close()

    logger.info("Plotting strains per 2nd stage cluster")
    set_labels_font_size()
    sorted_stats = clusters_df.sort_values('total_strains', ascending=False).reset_index(drop=True)
    plt.bar(sorted_stats.index.values, sorted_stats['total_strains'], width)
    plt.xlabel("Clusters")
    plt.ylabel("Strains")
    plt.title("Total strains per 2nd stage cluster")
    plt.savefig('strains_per_2nd_stage_cluster.pdf', format="pdf")
    plt.close()

    logger.info("Plotting strains per 2nd stage clusters without protein sequences")
    set_labels_font_size()
    clusters_without_reps = clusters_df[clusters_df['1st_stage_reps'] == 0]
    sorted_stats = clusters_without_reps.sort_values('total_strains', ascending=False).reset_index(drop=True)
    plt.bar(sorted_stats.index.values, sorted_stats['total_strains'], width)
    plt.xlabel("Clusters")
    plt.ylabel("Strains")
    plt.title("Total strains per 2nd stage cluster without protein sequences")
    plt.savefig('strains_per_2nd_stage_cluster_without_protein_sequences.pdf', format="pdf")
    plt.close()

    logger.info("Plotting strains in protein sequence's 1st stage cluster VS pseudogenes in 2nd stage cluster")
    set_labels_font_size()
    clusters_with_reps = clusters_df[clusters_df['1st_stage_reps'] == 1]
    sorted_stats = clusters_with_reps.sort_values('strains_in_rep_1st_stage_cluster', ascending=False).reset_index(drop=True)
    chart1 = plt.bar(sorted_stats.index.values, sorted_stats['strains_in_rep_1st_stage_cluster'], width)
    chart2 = plt.bar(sorted_stats.index.values + width, sorted_stats['total_strains'] - 1, width)
    plt.xlabel("Clusters")
    plt.ylabel("Strains in protein rep 1st stage cluster /\nPseudogenes")
    plt.title("Strains in protein rep 1st stage cluster VS\nPseudogenes per 2nd stage cluster")
    plt.legend((chart1[0], chart2[0]), ("Strains in protein rep\n1st stage cluster", "Pseudogenes in 2nd\nstage cluster"))
    plt.savefig('protein_rep_1st_cluster_strains_vs_pseudogenes_per_2nd_stage_cluster.pdf', format="pdf")
    plt.close()


def set_labels_font_size():
    ax = plt.subplot()
    for label in (ax.get_xticklabels() + ax.get_yticklabels()):
        label.set_fontsize(3)
