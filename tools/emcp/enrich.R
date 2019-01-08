# parse parameter ---------------------------------------------------------
library(argparser, quietly=TRUE)
# Create a parser
p <- arg_parser("do GO enrichment")

# Add command line arguments
p <- add_argument(p, "--orgdb_dir", help="orgDB build from emapper ", type="character")
p <- add_argument(p, "--de_result", help="input de_result file, from run_DE_analysis.pl", type="character")
p <- add_argument(p, "--de_log2FoldChange", help="log2FoldChange cutoff", type="numeric", default = 1)
p <- add_argument(p, "--de_padj", help="adjust pvalue cutoff", type="numeric", default = 0.05)

# Parse the command line arguments
argv <- parse_args(p)

out_prefix <- argv$de_result

# load library ------------------------------------------------------------
library(tidyverse)
library(clusterProfiler)
library(DOSE)

is.installed <- function(mypkg) is.element(mypkg, installed.packages()[,1]) 
if (! is.installed("org.My.eg.db"))
	install.packages(paste(argv$orgdb_dir, "org.My.eg.db", sep = "/"), repos = NULL)

library(org.My.eg.db)


# load gene list or de_result ---------------------------------------------
de_result <- read.delim(argv$de_result)

de_result_filter <- mutate(de_result, GID = rownames(de_result)) %>%
  filter(abs(log2FoldChange) > argv$de_log2FoldChange, pvalue < argv$de_padj)

gene_list <- as.character(de_result_filter$GID)


# do GO enrich -------------------------------------------------------------
ego <- enrichGO(gene          = gene_list, #差异基因 vector
                keyType       = "GID",  #差异基因的 ID 类型，需要是 OrgDb 支持的
                OrgDb         = org.My.eg.db, #对应的OrgDb
                ont           = "CC", #GO 分类名称，CC BP MF 
                pvalueCutoff  = 1, #Pvalue 阈值
                qvalueCutoff  = 1, #qvalue 阈值
                pAdjustMethod = "BH", #Pvalue 矫正方法
                readable      = FALSE) #TRUE 则展示SYMBOL，FALSE 则展示原来的ID

#将 ego 对象转换为dataframe，新版本可以用as.data.frame(ego)
ego_results<-as.data.frame(ego)
write.table(ego_results, file = paste(out_prefix, "ego_results.txt", sep = "."),
            quote = F)

pdf(file = paste(out_prefix, "ego_barplot.pdf", sep = "."))
barplot(ego, showCategory=20, x = "GeneRatio")
dev.off()

pdf(file = paste(out_prefix, "ego_dotplot.pdf", sep = "."))
dotplot(ego)
dev.off()

pdf(file = paste(out_prefix, "ego_emapplot.pdf", sep = "."))
emapplot(ego)
dev.off()

pdf(file = paste(out_prefix, "ego_goplot.pdf", sep = "."))
goplot(ego)
dev.off()


# Do Pathway enrich ------------------------------------------------------
# 从 OrgDB 提取 Pathway 和基因的对应关系
pathway2gene <- select(org.My.eg.db, keys = keys(org.My.eg.db), columns = c("Pathway")) %>%
  na.omit() %>%
  dplyr::select(Pathway, GID)

# 导入 Pathway 与名称对应关系
load("tools/emcp/kegg_info.RData")

#KEGG pathway 富集
ekp <- enricher(gene_list, 
                TERM2GENE = pathway2gene, 
                TERM2NAME = pathway2name, 
                pvalueCutoff = 1, 
                qvalueCutoff = 1,
                pAdjustMethod = "BH",
                minGSSize = 1)

ekp_results <- as.data.frame(ekp)
write.table(ekp_results, file = paste(out_prefix, "ekp_results.txt", sep = "."),
            quote = F)

pdf(file = paste(out_prefix, "ekp_barplot.pdf", sep = "."))
barplot(ekp, showCategory=20, x = "GeneRatio")
dev.off()

pdf(file = paste(out_prefix, "ekp_dotplot.pdf", sep = "."))
dotplot(ekp)
dev.off()

pdf(file = paste(out_prefix, "ekp_emapplot.pdf", sep = "."))
emapplot(ekp)
dev.off()
