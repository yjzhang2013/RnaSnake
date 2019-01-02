my_orgdb <- makeOrgPackageFromEmapper("input/sesame.emapper.annotations", 
                                      "zhangxudong <zhangxudong@genek.tv>", 
                                      tax_id = "4182", 
                                      genus = "Sesamum", 
                                      species = "indicum")

if (requireNamespace(my_orgdb, quietly = TRUE))
  remove.packages(my_orgdb)
install.packages(my_orgdb, repos = NULL)

library(org.Sindicum.eg.db)
columns(org.Sindicum.eg.db)