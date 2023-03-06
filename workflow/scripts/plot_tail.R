#!/usr/bin/env Rscript
library(tidyverse)
library(data.table)
library(ggpubr)
library(grid)
library(gridExtra)

## TO DO: Change name of script and check if ok for github

######## ARGUMENTS
args = commandArgs(trailingOnly=TRUE)

#suffix_add_tail=args[1]
suffix_add_tail=".read_info.result.merged.parts.csv"
#dir_add_tail=args[2]
dir_add_tail="/home/jpeter/DATA/FLEPseq/RUN08_Fleurs_v3/4_Tail"
#sample_corresp=args[3]
sample_corresp="/home/jpeter/DATA/FLEPseq/RUN08_Fleurs_v3/barcode_correspondances.tsv"

######## \ ARGUMENTS

######## FUNCTIONS
give.n <- function(x){
  return(c(y = -1, label = length(x)))
  # experiment with the multiplier to find the perfect position
}
######## \ FUNCTIONS

######## DATA IMPORT

samples_infos <- fread(sample_corresp, header = F, col.names = c("code", "sample")) %>%
  mutate(add_tail_path=file.path(dir_add_tail,paste0(code,suffix_add_tail)))

nlist_add_tail <- samples_infos$add_tail_path

names(nlist_add_tail) <- samples_infos$code

REF_genotype <- as.character(samples_infos[1, "sample"])

df_uri <- rbindlist(lapply(nlist_add_tail, fread), idcol = "code") %>%
  left_join(samples_infos, by = "code") %>%
  separate(read_core_id, into=c(NA, "chr", "read_start", "read_end"), sep = ",") %>%
  mutate(U_state= case_when(
    add_tail_pct_T>0.70 ~ "U-tail",
    TRUE ~ "non-U"))

######## \ DATA IMPORT

######## URIDYLATION

df_uri$tail_length <- nchar(df_uri$additional_tail)


mRNA_pctU <- df_uri%>% group_by(sample, U_state, mRNA, .drop=FALSE) %>%
  summarise(nb_reads=n())  %>%
  group_by(sample, mRNA) %>%
  mutate(total=sum(nb_reads), Percent=100*nb_reads/sum(nb_reads)) %>%
  filter(total>=50) %>%
  arrange(sample, mRNA)

global_pctU <- df_uri %>% group_by(sample, U_state, .drop=FALSE) %>%
  summarise(nb_reads=n())  %>%
  group_by(sample) %>%
  mutate(total=sum(nb_reads), Percent=100*nb_reads/sum(nb_reads))

mRNA_pctU$sample <- factor(mRNA_pctU$sample, levels=as.vector(samples_infos$sample))
global_pctU$sample <- factor(global_pctU$sample, levels=as.vector(samples_infos$sample))

my_colors <- c("gray50", "darkblue", "wheat", "darkred")

mRNA_Utails <- mRNA_pctU %>% filter(U_state=="U-tail")
  
global_Utails <- global_pctU%>% filter(U_state=="U-tail")

p1 <- ggplot(mRNA_Utails, aes(x=sample, y=Percent, fill=sample)) + 
  geom_boxplot(outlier.shape = '.', show.legend = FALSE) +
  stat_compare_means(label = "p.format", method = "wilcox.test",ref.group = REF_genotype, size = 2.9) +
  stat_summary(fun.data = give.n, geom = "text", fun = median, size=3) +
  ggtitle("Distribution of genes", subtitle = paste0("Wilcoxon tests of each sample versus ", REF_genotype)) +
  scale_fill_manual(values = my_colors) +
  theme_bw() +
  theme(axis.text.x=element_text(angle=45, hjust=1)) +
  theme(plot.title = element_text(size = 12), plot.subtitle = element_text(size=10))

p2 <- ggplot(global_Utails, aes(x=sample, y=Percent, fill=sample)) + 
  geom_bar(stat="identity", color="black", show.legend = FALSE) +
  ggtitle("Global percentage") +
  scale_fill_manual(values = my_colors) +
  theme_bw() +
  theme(axis.text.x=element_text(angle=45, hjust=1)) +
  theme(plot.title = element_text(size = 12))


p <- grid.arrange(p2, p1, ncol=2, nrow = 1, widths=c(3,5), heights=c(4), top = textGrob(paste("Percent of uridylation of",paste(as.character(as.vector(samples_infos$sample)), collapse=", ")),gp=gpar(fontsize=14,font=2)))

ggsave(filename =file.path(dir_add_tail,"AddTail_Percent_of_uridylation.pdf"), p, width = 8, height = 4, dpi = 300)

df_uri$sample <- factor(df_uri$sample, levels=as.vector(samples_infos$sample))

p <- ggplot(df_uri%>%filter(U_state=="U-tail")) + geom_bar(aes(tail_length, fill=sample), color="black") +
  facet_wrap(~sample, ncol=1, scales="free") +
  #scale_x_continuous(limits=c(0,20), breaks=seq(0,20,by=2)) +
  scale_fill_manual(values = my_colors) +
  theme_bw()+
  ggtitle("Read counts vs Utail length")

p

ggsave(filename =file.path(dir_add_tail,"AddTail_length_barplot.pdf"), p, width = 5, height = 6, dpi = 300)

df_uri_add_tail_long <- df_uri %>% filter(U_state=="U-tail") %>%
  pivot_longer(cols=c("add_tail_pct_A", "add_tail_pct_C", "add_tail_pct_G", "add_tail_pct_T"), names_to = "add_tail_nucl", values_to = "add_tail_pct")

df_uri_add_tail_long_summ <- df_uri_add_tail_long %>%
  group_by(sample, add_tail_nucl) %>%
  summarise(mean_add_tail_pct=mean(add_tail_pct))

df_uri_add_tail_long$sample <- factor(df_uri_add_tail_long$sample, levels=as.vector(samples_infos$sample))
df_uri_add_tail_long_summ$sample <- factor(df_uri_add_tail_long_summ$sample, levels=as.vector(samples_infos$sample))

p <- ggplot(df_uri_add_tail_long_summ, aes(x=sample, y=mean_add_tail_pct, fill=sample)) + 
  geom_bar(stat='identity') + facet_wrap(~add_tail_nucl) +
  scale_fill_manual(values = my_colors) +
  theme_bw() +
  ggtitle("additional tail base composition")

p

ggsave(filename =file.path(dir_add_tail,"Addtail_BaseComposition.pdf"), p, width = 5, height = 6, dpi = 300)


p <- ggplot(df_uri_add_tail_long, aes(x=sample, y=tail_length, fill=sample)) + 
  geom_boxplot(outlier.shape = '.') +
  scale_fill_manual(values = my_colors) +
  theme_bw() +
  coord_cartesian(ylim=c(0, 25))

p
ggsave(filename =file.path(dir_add_tail,"Utail_length_bp.pdf"), p, width = 5, height = 6, dpi = 300)


######## \URIDYLATION


######## PolyA

df_uri$polya_length_nchar <- nchar(df_uri$polytail)

ggplot(df_uri %>% filter(polya_length_nchar<200,
                         polya_length_nchar>10)) +
  geom_density(aes(polya_length_nchar, fill=sample, color=sample), alpha=0.2,  lwd = 1) +
  facet_wrap(~U_state, scales="free") +
  scale_fill_manual(values = my_colors) +
  scale_color_manual(values = my_colors) +
  theme_bw() +
  ggtitle("Distribution of Poly-A tail sizes", subtitle = "polyA size=number of characters PolyA")

ggsave(filename =file.path(dir_add_tail,"PolyA_length_nchar.pdf"),width = 5, height = 6, dpi = 300)


ggplot(df_uri %>% filter(polya_length_nchar<200,
                     polya_length_nchar>10)) +
  geom_density(aes(polya_length, fill=sample, color=sample), alpha=0.2,  lwd = 1) +
  facet_wrap(~U_state, scales="free") +
  scale_fill_manual(values = my_colors) +
  scale_color_manual(values = my_colors) +
  theme_bw() +
  ggtitle("Distribution of Poly-A tail sizes", subtitle = "polyA size=polya_length from FLEPSeq2")

ggsave(filename =file.path(dir_add_tail,"PolyA_length_polyAlength.pdf"), width = 5, height = 6, dpi = 300)

ggplot(df_uri %>% filter(polya_length_nchar<200,
                         polya_length_nchar>10)) +
  geom_density(aes(init_polya_length, fill=sample, color=sample), alpha=0.2,  lwd = 1) +
  facet_wrap(~U_state, scales="free") +
  scale_fill_manual(values = my_colors) +
  scale_color_manual(values = my_colors) +
  theme_bw() +
  ggtitle("Distribution of Poly-A tail sizes", subtitle = "polyA size=init_polya_length from FLEPSeq2")

ggsave(filename =file.path(dir_add_tail,"PolyA_length_initPolyALength.pdf"), width = 5, height = 6, dpi = 300)


######## \ PolyA



