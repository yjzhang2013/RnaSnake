configfile: "config.yaml"
rule all:
    input:
        config["Assembly_Dir"] + "/corset.fasta",
	config["Assembly_Dir"] + "/corset.fasta.txt",
        config["Assembly_Dir"] + "/corset-clusters_trans_map.txt",
        config["Assembly_Dir"] + "/SuperDuper.fasta",
        config["Assembly_Dir"] + "/run_busco/short_summary_busco.txt",
        config["Annotation_Dir"] + "/my.emapper.annotations",
        config["Annotation_Dir"] + "/anno_stat.txt",
        directory(config["Annotation_Dir"] + "/org.My.eg.db"),
        config["Quantification_Dir"] + "/my.gene.counts.matrix",
        config["Quantification_Dir"] + "/my.gene.TMM.EXPR.matrix",
        config["ExprAnalysis_Dir"] + "/sample_cor/my.minRow10.sample_cor.dat",
        config["ExprAnalysis_Dir"] + "/sample_cor/my.minRow10.sample_cor_matrix.pdf",
        config["ExprAnalysis_Dir"] + "/pca/my.minRow10.PCA.prcomp.scores",
        config["ExprAnalysis_Dir"] + "/pca/my.minRow10.prcomp.principal_components.pdf",
        dynamic(config["ExprAnalysis_Dir"] + "/deg/my.gene.counts.matrix.{Case_vs_Control}.DESeq2.DE_results"),
        config["ExprAnalysis_Dir"] + "/deg/my.gene.counts.matrix.Case_vs_Control.DESeq2.DE_results.ekp_results.txt"


######################################################
# Includes
######################################################

include: "Includes/Assembly.snake"
include: "Includes/Annotation.snake"
include: "Includes/Quantification.snake"
include: "Includes/ExprAnalysis.snake"
include: "Includes/Enrichment.snake"
