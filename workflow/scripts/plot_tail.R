library(tidyverse)
library(data.table)
library(ggpubr)
library(grid)
library(gridExtra)

## TO DO: REMOVE POLY A SECTION AS ALL THE INFOS WILL BE ON THE URI FILE

######## ARGUMENTS
#suffix_polyA=".read_info.result.merged.txt"
suffix_add_tail=".read_info.result.merged.parts.csv"

#dir_polyA="/home/jpeter/DATA/FLEPseq/RUN08_Fleurs_v2/3_PolyA/"
dir_add_tail="/home/jpeter/DATA/FLEPseq/RUN08_Fleurs_v3/4_Tail/"

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
  mutate(add_tail_path=paste0(dir_add_tail,code,suffix_add_tail ))

nlist_add_tail <- samples_infos$add_tail_path

names(nlist_add_tail) <- samples_infos$code

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
  stat_compare_means(label = "p.format", method = "t.test",ref.group = "WT", size = 2.9) +
  stat_summary(fun.data = give.n, geom = "text", fun = median, size=3) +
  ggtitle("Distribution of genes", subtitle = "T-tests of each sample versus WT") +
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

ggsave(filename =file.path(dir_add_tail,"Percent_of_uridylation.pdf"), p, width = 8, height = 4, dpi = 300)

df_uri$sample <- factor(df_uri$sample, levels=as.vector(samples_infos$sample))


p <- ggplot(df_uri%>%filter(U_state=="U-tail")) + geom_bar(aes(tail_length, fill=sample), color="black") +
  facet_wrap(~sample, ncol=1, scales="free") +
  scale_x_continuous(limits=c(0,20), breaks=seq(0,20,by=2)) +
  scale_fill_manual(values = my_colors) +
  theme_bw()+
  ggtitle("Read counts vs Utail length")

p

ggsave(filename =file.path(dir_add_tail,"Utail_seeds.pdf"), p, width = 5, height = 6, dpi = 300)


######## \URIDYLATION


######## PolyA

ggplot(df_uri) + geom_density(aes(polya_length, fill=Sample, color=Sample), alpha=0.2,  lwd = 1) +
  facet_wrap(~Uri, scales="free") +
  scale_x_continuous(limits=c(10,200), breaks=seq(10,200,by=20)) +
  scale_fill_manual(values = my_colors) +
  scale_color_manual(values = my_colors) +
  theme(
    panel.spacing = unit(0.4, "lines"),
    plot.title = element_text(hjust = 0.5),
    panel.background = element_blank(),
    panel.grid.major = element_line(size = 0.2, linetype = 'dashed', colour = "gray70"),
    panel.border = element_rect(colour="black",fill=NA,size=0.5),
    legend.position = "none" ) 

######## \ PolyA



