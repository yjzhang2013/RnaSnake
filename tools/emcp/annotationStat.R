
# load library ------------------------------------------------------------

library(org.Sindicum.eg.db)
columns(org.Sindicum.eg.db)


library(GO.db)
head(keys(GO.db))

columns(GO.db)

go_info <- select(GO.db, keys = head(keys(GO.db)), columns = c("DEFINITION", "GOID","ONTOLOGY","TERM")) %>%
  dplyr::select(GO = GOID, GO_Name = TERM, GO_Class = ONTOLOGY)
