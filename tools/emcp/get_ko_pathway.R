library(tibble)
library(stringr)
library(KEGGREST)
kos <- names(keggList("ko"))


ko_pathway_link <- c()
for (ko in kos) {
  v <- keggLink('pathway', "K00844")
  t <- tibble(KO = str_sub(names(info), 4,9), 
              Pathway = str_sub(info,6, 13))
  write.table(t, file = "ko_pathway.txt", 
              row.names = F, 
              quote = F, 
              col.names = F, 
              append = F)
}

#save(ko_pathway_link, file = "ko_pathway_link.RData")
