import logging
import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt

logger = logging.getLogger(__name__)


def create_1st_stage_charts(stats_df):
    logger.info("Plotting clusters per strain")
    set_labels_font_size()
    sorted_stats = stats_df.sort_values('total_clusters', ascending=False).reset_index(drop=True)
    plt.bar(sorted_stats.index.values, sorted_stats['total_clusters'])
    plt.xlabel("strains #")
    plt.ylabel("clusters #")
    plt.title("strains to clusters bar chart")
    plt.savefig('total_clusters_by_strain.pdf', format="pdf")
    plt.close()

    logger.info("Plotting missing core clusters per strain")
    set_labels_font_size()
    sorted_stats = stats_df.sort_values('missing_core', ascending=False).reset_index(drop=True)
    plt.bar(sorted_stats.index.values, sorted_stats['missing_core'])
    plt.xlabel("strains #")
    plt.ylabel("% of missing core clusters")
    plt.title("strains to missing core clusters bar chart")
    plt.savefig('missing_core_clusters_by_strain.pdf', format="pdf")
    plt.close()

    logger.info("Plotting singleton clusters per strain")
    set_labels_font_size()
    sorted_stats = stats_df.sort_values('singletons', ascending=False).reset_index(drop=True)
    plt.bar(sorted_stats.index.values, sorted_stats['singletons'])
    plt.xlabel("strains #")
    plt.ylabel("singleton clusters #")
    plt.title("strains to singleton clusters bar chart")
    plt.savefig('singleton_clusters_by_strain.pdf', format="pdf")
    plt.close()

    logger.info("Plotting pseudogenes per strain")
    set_labels_font_size()
    sorted_stats = stats_df.sort_values('pseudogenes', ascending=False).reset_index(drop=True)
    plt.bar(sorted_stats.index.values, sorted_stats['pseudogenes'])
    plt.xlabel("strains #")
    plt.ylabel("pseudogenes #")
    plt.title("strains to pseudogenes bar chart")
    plt.savefig('pseudogenes_by_strain.pdf', format="pdf")
    plt.close()

    logger.info("Plotting contigs per strain")
    set_labels_font_size()
    sorted_stats = stats_df.sort_values('contigs', ascending=False).reset_index(drop=True)
    plt.bar(sorted_stats.index.values, sorted_stats['contigs'])
    plt.xlabel("strains #")
    plt.ylabel("contigs #")
    plt.title("strains to contigs bar chart")
    plt.savefig('contigs_by_strain.pdf', format="pdf")
    plt.close()

    logger.info("Plotting contigs VS singletons per strain")
    set_labels_font_size()
    sorted_stats = stats_df.sort_values('contigs', ascending=False).reset_index(drop=True)
    plt.bar(sorted_stats.index.values, sorted_stats['contigs'])
    plt.bar(sorted_stats.index.values, sorted_stats['singletons'])
    plt.xlabel("strains #")
    plt.ylabel("contigs #")
    plt.title("strains to contigs VS singletons bar chart")
    plt.savefig('contigs_vs_singletons_by_strain.pdf', format="pdf")
    plt.close()

    logger.info("Plotting contigs VS missing core % per strain")
    set_labels_font_size()
    sorted_stats = stats_df.sort_values('contigs', ascending=False).reset_index(drop=True)
    plt.bar(sorted_stats.index.values, sorted_stats['contigs'])
    plt.bar(sorted_stats.index.values, sorted_stats['missing_core'])
    plt.xlabel("strains #")
    plt.ylabel("contigs #")
    plt.title("strains to contigs VS missing core bar chart")
    plt.savefig('contigs_vs_missing_core_by_strain.pdf', format="pdf")
    plt.close()

    logger.info("Plotting contigs VS pseudogenes per strain")
    set_labels_font_size()
    sorted_stats = stats_df.sort_values('contigs', ascending=False).reset_index(drop=True)
    plt.bar(sorted_stats.index.values, sorted_stats['pseudogenes'])
    plt.bar(sorted_stats.index.values, sorted_stats['contigs'])
    plt.xlabel("strains #")
    plt.ylabel("contigs #")
    plt.title("strains to contigs VS pseudogenes bar chart")
    plt.savefig('contigs_vs_pseudogenes_by_strain.pdf', format="pdf")
    plt.close()

    logger.info("Plotting singletons VS missing core % per strain")
    set_labels_font_size()
    sorted_stats = stats_df.sort_values('singletons', ascending=False).reset_index(drop=True)
    plt.bar(sorted_stats.index.values, sorted_stats['singletons'])
    plt.bar(sorted_stats.index.values, sorted_stats['missing_core'])
    plt.xlabel("strains #")
    plt.ylabel("singletons #")
    plt.title("strains to singletons VS missing core bar chart")
    plt.savefig('singletons_vs_missing_core_by_strain.pdf', format="pdf")
    plt.close()

    logger.info("Plotting pseudogenes VS missing core % per strain")
    set_labels_font_size()
    sorted_stats = stats_df.sort_values('pseudogenes', ascending=False).reset_index(drop=True)
    plt.bar(sorted_stats.index.values, sorted_stats['pseudogenes'])
    plt.bar(sorted_stats.index.values, sorted_stats['missing_core'])
    plt.xlabel("strains #")
    plt.ylabel("pseudogenes #")
    plt.title("strains to pseudogenes VS missing core bar chart")
    plt.savefig('pseudogenes_vs_missing_core_by_strain.pdf', format="pdf")
    plt.close()

    logger.info("Plotting pseudogenes VS singletons per strain")
    set_labels_font_size()
    sorted_stats = stats_df.sort_values('pseudogenes', ascending=False).reset_index(drop=True)
    plt.bar(sorted_stats.index.values, sorted_stats['pseudogenes'])
    plt.bar(sorted_stats.index.values, sorted_stats['singletons'])
    plt.xlabel("strains #")
    plt.ylabel("pseudogenes #")
    plt.title("strains to pseudogenes VS singletons bar chart")
    plt.savefig('pseudogenes_vs_singletons_by_strain.pdf', format="pdf")
    plt.close()


def set_labels_font_size():
    ax = plt.subplot()
    for label in (ax.get_xticklabels() + ax.get_yticklabels()):
        label.set_fontsize(3)
