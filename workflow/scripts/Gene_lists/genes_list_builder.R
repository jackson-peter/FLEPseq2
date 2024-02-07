library(tidyverse)
library(data.table)
library(rtracklayer)

##### ARGUMENTS #####
args = commandArgs(trailingOnly=TRUE)

# annotation file in gtf format
gtf <- args[1]

# subset list of genes (repr gene model or canonical list...)
gene_list_f <- args[2]
gene_list <- fread(gene_list_f) %>%
  select(c(V4))
colnames(gene_list) <- c("AGI")
# output directory
outdir <- args[3]

# basename for gene lists files for flepseq
suffix <- args[4]

# outfiles names
out_ie <- file.path(outdir, paste0("intron_exon_", suffix, ".bed"))
out_si <- file.path(outdir, paste0("select_intron_", suffix, ".bed"))


##### FUNCTIONS #####

# allows to shift up all line from a column n times. this is usefull here because the start on an intron is the end of the previous exon.
# the same for the end of intron: start of next exon
shift <- function(x, n){
  c(x[-(seq(n))], rep(NA, n))
}

# main function here
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




build_intronexon(gtf, gene_list$AGI, out_ie, out_si)

