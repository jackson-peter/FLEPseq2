library(tidyverse)
library(data.table)
library(ggvenn)
library(rtracklayer)

shift <- function(x, n){
  c(x[-(seq(n))], rep(NA, n))
}


#
repr_genes_ara11 <- fread("~/DATA/ReferenceGenomes/Athaliana/TAIR10/araport11.representative.gene_model.bed") %>%
  select(c(V4))
colnames(repr_genes_ara11) <- c("AGI")

#
repr_genes_tair10 <- fread("~/DATA/ReferenceGenomes/Athaliana/TAIR10/TAIR10_representative_gene_models", col.names = c("AGI"), skip =4 ) %>%
  filter(str_detect(AGI, "^AT"))

#
genes_flepseq <- fread("/home/jpeter/DATA/ReferenceGenomes/Athaliana/TAIR10/exon_intron_pos.repr.bed") %>%
  select(V4)
colnames(genes_flepseq) = c("AGIlong")
genes_flepseq <- genes_flepseq %>%
  separate_wider_delim(AGIlong, "_", names=c("AGI", "type")) %>%
  select(c(AGI))
flepseq_genes <- unique(genes_flepseq$AGI)

# Chart
x = list(ara11_repr=repr_genes_ara11$AGI, tair10_repr=repr_genes_tair10$AGI, FLEPSEQ=flepseq_genes)
vennplot <- ggvenn(
  x,
  fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF"),
  stroke_size = 0.5, set_name_size = 4
)


ggsave(vennplot, file="~/DATA/FLEPseq/Genes_list/overlap_list_reprGeneModel_original.pdf")


only_in_flepseq <- as.data.frame(setdiff(unique(genes_flepseq$AGI),c(repr_genes_ara11$AGI, repr_genes_tair10$AGI)))
write_tsv(as.data.frame(only_in_flepseq), file= "~/DATA/FLEPseq/Genes_list/only_in_flepseq.tsv")

only_in_ara11 <- as.data.frame(setdiff(repr_genes_ara11$AGI,c(unique(genes_flepseq$AGI), repr_genes_tair10$AGI)) )

write_tsv(as.data.frame(only_in_ara11), file= "~/DATA/FLEPseq/Genes_list/only_in_ara11.tsv")

only_in_tair10 <- as.data.frame(setdiff(repr_genes_tair10$AGI,c(unique(genes_flepseq$AGI), repr_genes_ara11$AGI)))
write_tsv(as.data.frame(only_in_ara11), file= "~/DATA/FLEPseq/Genes_list/only_in_tair10.tsv")

gtf_ara11_only <- as_tibble(import("/home/jpeter/DATA/ReferenceGenomes/Athaliana/Araport11/Araport11_GTF_genes_transposons.current.gtf")) %>%
  filter()



#
gtf_ara11_exon <- as_tibble(import("/home/jpeter/DATA/ReferenceGenomes/Athaliana/Araport11/Araport11_GTF_genes_transposons.current.gtf")) %>%
  filter(type %in%  "exon") %>%
  mutate(start=start-1) %>%
  filter(transcript_id %in% repr_genes_ara11$AGI) %>%
  group_by(transcript_id, type) %>%
  mutate(group_id = case_when(strand=="+" ~ row_number(),
                              strand=='-' ~ rev(row_number())),
         nb_exon=max(group_id))



gtf_ara11_exon %>% group_by(gene_id) %>%summarise()

gtf_ara11_intron_plus <- gtf_ara11_exon %>%
  filter(strand=='+') %>%
  mutate(intron_start=end, type="intron")
gtf_ara11_intron_plus$intron_end <- shift(gtf_ara11_intron_plus$start, 1)

gtf_ara11_intron_plus <- gtf_ara11_intron_plus %>%
  filter(group_id!=nb_exon)%>%
  group_by(transcript_id) %>%
  arrange(transcript_id, group_id)


gtf_ara11_intron_minus <- gtf_ara11_exon %>%
  filter(strand=='-') %>%
  mutate(intron_start=end, type="intron")
gtf_ara11_intron_minus$intron_end <- shift(gtf_ara11_intron_minus$start, 1)

gtf_ara11_intron_minus <- gtf_ara11_intron_minus%>%
  mutate(group_id= group_id -1) %>%
  filter(group_id!=0)%>%
  group_by(transcript_id) %>%
  arrange(transcript_id, desc(group_id), type)


gtf_ara11_intron <- rbind(gtf_ara11_intron_minus, gtf_ara11_intron_plus) %>%
  select(-c(start, end)) %>%
  rename("intron_start"="start") %>%
  rename("intron_end"="end")

gtf_ara11_intronexon <- rbind(gtf_ara11_exon, gtf_ara11_intron) %>%
  mutate(type_id=paste0(type, group_id))  %>%
  arrange(transcript_id, start) %>%
  mutate(reg_name=paste(transcript_id, type_id, sep='_'),
         point='.') %>%
  ungroup() %>%
  select(c(seqnames, start, end, point, strand, reg_name, transcript_id))



x = list(ara11_repr=repr_genes_ara11$AGI, tair10_repr=repr_genes_tair10$AGI, intronexon_11_repr=gtf_ara11_intronexon%>%pull(transcript_id))
vennplot <- ggvenn(
  x,
  fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF"),
  stroke_size = 0.5, set_name_size = 4
)

ggsave(vennplot, file="~/DATA/FLEPseq/Genes_list/overlap_list_reprGeneModel.pdf")

