# parse parameter ---------------------------------------------------------
library(argparser, quietly=TRUE)
# Create a parser
p <- arg_parser("expression matrix statistics")

# Add command line arguments
p <- add_argument(p, "--dat", help="expression matrix", type="character")
p <- add_argument(p, "--out", help="output prefix", type="character")

# Parse the command line arguments
argv <- parse_args(p)

dat<-argv$dat
outprefix<-argv$out

library(ggplot2)
library(tidyverse)
library(formattable)

## read data and melt

expr <- read_delim(dat, 
                   "\t", escape_double = FALSE, trim_ws = TRUE)

expr <- rename(expr, featureId = X1)

expr_tidy<-gather(expr, key = "sample", value = "expression", -1)

# density
f<-paste(outprefix, "density.pdf", sep = ".")
p<-qplot(log10(expression), data=expr_tidy, geom="density", xlim=c(-2,6), colour=sample, fill=sample, alpha=I(1/10), main = "Density distribution of transcript expression") + 
  theme_bw() + 
  scale_y_continuous(expand = c(0,0))

pdf(file=f)
p
dev.off()

# boxplot
p<-ggplot(expr_tidy, aes(x=sample, y=log10(expression), fill=sample)) + geom_boxplot() +
  theme_bw() +  
  theme(axis.text.x=element_text(hjust=1,angle=60)) + scale_y_continuous(expand = c(0,0))

f<-paste(outprefix, "boxplot.pdf", sep = ".")
pdf(file=f)
p
dev.off()

# expressed gene
total_number = length(expr$featureId)
expressed_tidy <- expr_tidy %>% filter(expression > 1) %>%
  group_by(sample) %>%
  summarise(expressed = n()) %>%
  mutate(percentage = percent(expressed/total_number))

f<-paste(outprefix, "stat.txt", sep = ".")
write.table(expressed_tidy, file = f, quote = F, row.names = F)
