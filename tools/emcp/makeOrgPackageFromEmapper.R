# parse parameter ---------------------------------------------------------
library(argparser, quietly=TRUE)
# Create a parser
p <- arg_parser("make OrgDb from emapper")

# Add command line arguments
p <- add_argument(p, "annotation", help="emapper annotation result", type="character")
p <- add_argument(p, "all_gene", help="all gene id list", type="character")
p <- add_argument(p, "out_dir", help="output directory", type="character")

# Parse the command line arguments
argv <- parse_args(p)
script_dir <- dirname(strsplit(commandArgs(trailingOnly = FALSE)[4],"=")[[1]][2])

# function:make OrgDB -----------------------------------------------------
library(tidyverse)
library(stringr)
library(AnnotationForge)
library(formattable)
library(clusterProfiler)
#' make OrgDB Package From eggnog-mapper result
#'
#' @param f_emapper_anno eggnog-mapper annotation result
#' @param author Who is the creator of this package? like "xxx <xxx@xxx.xx>"
#' @param tax_id The Taxonomy ID that represents your organism. (NCBI has a nice online browser for finding the one you need)
#' @param genus Single string indicating the genus
#' @param species Single string indicating the species
#'
#' @return OrgDb name
#' @export
#'
#' @examples
makeOrgPackageFromEmapper <- function(f_emapper_anno, 
                                      author, 
                                      tax_id = "0", 
                                      genus = "default", 
                                      species = "default") {
  
  # read emapper result
  emapper <- read_delim(f_emapper_anno,
                        "\t", escape_double = FALSE, trim_ws = TRUE)
  
  # extract gene name from emapper ------------------------------------------
  gene_info <- emapper %>%
    dplyr::select(GID = query_name, GENENAME = `eggNOG annot`) %>%
    na.omit()
  
  # extract go annotation from emapper --------------------------------------
  gos <- emapper %>%
    dplyr::select(query_name, GO_terms) %>%
    na.omit()
  
  gene2go = data.frame(GID = character(),
                       GO = character(),
                       EVIDENCE = character())
  
  df_temp <- list()
  for (row in 1:nrow(gos)) {
    the_gid <- gos[row, "query_name"][[1]]
    the_gos <- str_split(gos[row,"GO_terms"], ",", simplify = FALSE)[[1]]
    
    df_temp[[row]] <- data_frame(GID = rep(the_gid, length(the_gos)),
                                 GO = the_gos,
                                 EVIDENCE = rep("IEA", length(the_gos)))
  }
  
  gene2go <- bind_rows(df_temp)

  # extract kegg pathway annotation from emapper ----------------------------
  gene2ko <- emapper %>%
    dplyr::select(GID = query_name, Ko = KEGG_KOs) %>%
    na.omit()
  
  load(file = paste(script_dir, "kegg_info.RData", sep = "/"))
  gene2pathway <- gene2ko %>% left_join(ko2pathway, by = "Ko") %>%
    left_join(pathway2name, by = "Pathway") %>%
    dplyr::select(GID, Pathway, Pathway_Name, Pathway_Class) %>%
    na.omit()
  
  
  # extract COG annotation from emapper -------------------------------------
  cog_info <- read_delim(paste(script_dir, "cog_funclass.tab", sep = "/"), 
                         "\t", escape_double = FALSE, trim_ws = TRUE)
  
  cogs <- emapper %>%
    dplyr::select(query_name, COG = `COG cat`) %>%
    na.omit()
  
  gene2cog = data.frame(GID = character(),
                        COG = character())
  
  df_temp <- list()
  for (row in 1:nrow(cogs)) {
    the_gid <- cogs[row, "query_name"][[1]]
    the_cogs <- str_trim(str_split(cogs[row,"COG"], ",", simplify = FALSE)[[1]])
    
    df_temp[[row]] <- data_frame(GID = rep(the_gid, length(the_cogs)),
                                 COG = the_cogs)
  }
  gene2cog <- bind_rows(df_temp)
  
  gene2cog <- gene2cog %>% left_join(cog_info, by = "COG")
  
  save(gene_info, gene2go, gene2pathway, gene2cog, file = "gene_annotation.RData")
  # make OrgDb --------------------------------------------------------------
  makeOrgPackage(gene_info=gene_info,
                 go=gene2go,
                 ko=gene2ko,
                 pathway=gene2pathway,
                 cog=gene2cog,
                 maintainer=author,
                 author=author,
                 outputDir=argv$out_dir,
                 tax_id=tax_id,
                 genus=genus,
                 species=species,
                 goTable="go",
                 version="1.0")
  
  my_orgdb <- str_c("org.", str_to_upper(str_sub(genus, 1, 1)) , species, ".eg.db", sep = "")
  return(my_orgdb)
}

# function: annotation statistics and plot --------------------------------
annoStat <- function() {
  
  # read all gene list
  all_gene <- as.character(read.table(argv$all_gene, quote="\"", comment.char="")$V1)
  load("gene_annotation.RData")
  
  # GO statistics and plot --------------------------------------------------
  
  go_bp <- groupGO(gene     = all_gene,
                   OrgDb    = org.My.eg.db,
                   keyType  = "GID",
                   ont      = "BP",
                   level    = 2,
                   readable = FALSE)
  
  go_bp <- as.data.frame(go_bp)
  go_bp$GO_Class <- "Biological Process"
  
  go_cc <- groupGO(gene     = all_gene,
                   OrgDb    = org.My.eg.db,
                   keyType  = "GID",
                   ont      = "CC",
                   level    = 2,
                   readable = FALSE)
  
  go_cc <- as.data.frame(go_cc)
  go_cc$GO_Class <- "Cellular Component"
  
  go_mf <- groupGO(gene     = all_gene,
                   OrgDb    = org.My.eg.db,
                   keyType  = "GID",
                   ont      = "MF",
                   level    = 2,
                   readable = FALSE)
  go_mf <- as.data.frame(go_mf)
  go_mf$GO_Class <- "Molecular Function"
  
  go_all <- rbind(go_bp, go_cc, go_mf)
  
  p <- ggplot(go_all) + 
    geom_histogram(aes(x = Description, 
                       y = Count,
                       fill = GO_Class),
                   stat = "identity") + facet_wrap(~GO_Class, scales = "free_x") + 
    labs(title = "GO function classification", y = "Number of genes") +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5),
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.4),
          axis.title.x = element_blank(),
          legend.position = "none")
  
  ggsave(paste(argv$out_dir, "go.pdf", sep = "/"), p, width = 20, height = 7)
  
  
  # Pathway statistics and plot ---------------------------------------------
  
  
  # COG statistics and plot -------------------------------------------------
  gene2cog$COG_Name = paste("(", gene2cog$COG, ")", gene2cog$COG_Name, sep = " ")
  
  p <- ggplot(data = gene2cog) + 
    geom_bar(aes(x = COG, 
                 fill = COG_Name)) +
    labs(title = "COG/KOG Function Classification ", 
         x = "",
         y = "Number of genes") +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5),
          legend.title = element_blank(),
          legend.key.size=unit(1,"line"),
          legend.text = element_text(size = 7.5)) +
    guides(fill=guide_legend(ncol=1))
  ggsave(paste(argv$out_dir, "cog.pdf", sep = "/"), p, width = 16, height = 7)
  
  
  # number and percentage ---------------------------------------------------
  total_gene = length(all_gene)
  eggnog_anno = length(gene_info$GID)
  go_anno = length(unique(gene2go$GID))
  cog_anno = length(unique(gene2cog$GID))
  pathway_anno = length(unique(gene2pathway$GID))
  
  anno_stat <- data_frame(
    database = c("EggNOG", "GO", "COG/KOG", "KEGG Pathway"),
    number = comma(c(eggnog_anno, go_anno, cog_anno, pathway_anno), digits = 0),
    percentage = percent(c(eggnog_anno, go_anno, cog_anno, pathway_anno)/total_gene)
  )
  
  write.table(anno_stat, paste(argv$out_dir, "anno_stat.txt", sep = "/"), quote = F, row.names = F, sep = "\t")
}

# makeOrgPackage ----------------------------------------------------------
my_orgdb <- makeOrgPackageFromEmapper(argv$annotation, 
                                      "test <test@genek.tv>", 
                                      tax_id = "0000", 
                                      genus = "M", 
                                      species = "y")

is.installed <- function(mypkg) is.element(mypkg, installed.packages()[,1]) 
if (is.installed(my_orgdb))
  remove.packages(my_orgdb)
install.packages(paste(argv$out_dir, my_orgdb, sep = "/"), repos = NULL)

library(my_orgdb, character.only = TRUE)

annoStat()
