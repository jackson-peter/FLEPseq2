library(tidyverse)
library(data.table)
library(rtracklayer)


shift <- function(x, n){
  c(x[-(seq(n))], rep(NA, n))
}

build_intronexon <- function(gtf, gene_list,out_ie, out_si) {
  gtf_df_exon <- as_tibble(import(gtf)) %>%
    filter(type %in%  "exon") %>%
    print(n=10) %>%
    mutate(start=start-1) %>%
    filter(transcript_id %in% gene_list) %>%
    group_by(transcript_id, type) %>%
    mutate(group_id = case_when(strand=="+" ~ row_number(),
                                strand=='-' ~ rev(row_number())),
           nb_exon=max(group_id))

  gtf_df_intron_plus <- gtf_df_exon %>%
    filter(strand=='+') %>%
    mutate(intron_start=end, type="intron")
  gtf_df_intron_plus$intron_end <- shift(gtf_df_intron_plus$start, 1)

  gtf_df_intron_plus <- gtf_df_intron_plus %>%
    filter(group_id!=nb_exon)%>%
    group_by(transcript_id) %>%
    arrange(transcript_id, group_id)


  gtf_df_intron_minus <- gtf_df_exon %>%
    filter(strand=='-') %>%
    mutate(intron_start=end, type="intron")
  gtf_df_intron_minus$intron_end <- shift(gtf_df_intron_minus$start, 1)

  gtf_df_intron_minus <- gtf_df_intron_minus%>%
    mutate(group_id= group_id -1) %>%
    filter(group_id!=0)%>%
    group_by(transcript_id) %>%
    arrange(transcript_id, desc(group_id), type)


  gtf_df_intron <- rbind(gtf_df_intron_minus, gtf_df_intron_plus) %>%
    print(n=10) %>%
    select(-c(start, end)) %>%
    dplyr::rename("start"="intron_start") %>%
    dplyr::rename("end"="intron_end")



  gtf_df_intronexon <- rbind(gtf_df_exon, gtf_df_intron) %>%
    mutate(type_id=paste0(type, group_id))  %>%
    arrange(transcript_id, start) %>%
    mutate(reg_name=paste(transcript_id, type_id, sep='_'),
           point='.') %>%
    ungroup() %>%
    select(c(seqnames, start, end, reg_name, point, strand))



  write_tsv(gtf_df_intronexon, file=out_ie, col_names = FALSE)
  write_tsv(gtf_df_intron%>%
              ungroup() %>%
              mutate(intron_id=paste0(transcript_id, '_', type, group_id)) %>%
              select(intron_id),
            file=out_si, col_names = FALSE)

}

gtf <- "/shared/home/jpeter/Data/ReferenceGenomes/A_thaliana/Araport11/Araport11_GTF_genes_transposons.current.gtf"
repr_genes_ara11 <- fread("/shared/home/jpeter/Data/ReferenceGenomes/A_thaliana/Araport11/araport11.representative.gene_model.bed") %>%
  select(c(V4))
colnames(repr_genes_ara11) <- c("AGI")

out_ie <- "/shared/home/jpeter/Scripts/FLEPseq2/Genes_list/araport11_repr_gene_model/intron_exon_ara11_repr.bed"
out_si <- "/shared/home/jpeter/Scripts/FLEPseq2/Genes_list/araport11_repr_gene_model/select_intron_ara11_repr.bed"

build_intronexon(gtf, repr_genes_ara11$AGI, out_ie, out_si)

