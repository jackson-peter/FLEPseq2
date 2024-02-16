library(tidyverse)
library(data.table)
library(rtracklayer)

##### ARGUMENTS #####
args = commandArgs(trailingOnly=TRUE)

# annotation file in gtf format
gtf <- args[1]
gtf <- "~/Data/ReferenceGenomes/S_cerevisiae/S288C_reference_genome_R64-1-1_20110203/Saccharomyces_cerevisiae.R64-1-1.111.gtf"
gtf_col <- args[2]
gtf_col <- "gene_id"

# output directory
outdir <- args[3]
outdir="~/Scripts/FLEPseq2/Genes_list/yeast_flepseq/"

# basename for gene lists files for flepseq
suffix <- args[4]
suffix <- "yeast"
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
  
  print("### building exon list from gtf")
  # get all the exons that belong to a gene in gene_list
  gtf_df_exon <- as_tibble(import(gtf)) %>% 
    #filter(seqnames %in% chroms) %>%
    filter(type %in%  "exon") %>% 
    mutate(start=start-1) %>%
    filter(tag=="Ensembl_canonical") %>%
    #    filter(gtf_col %in% gene_list) %>%
    group_by(!!as.name(gtf_col), type) %>%
    mutate(group_id = case_when(strand=="+" ~ row_number(),
                                strand=='-' ~ rev(row_number())),
           nb_exon=max(group_id))
  
  #gtf_df_exon %>% filter(gene_id=="YDR129C")
  

  print("### building intron_plus list from gtf")
  
  # create introns entries using nb exons and strand
  # strand +
  gtf_df_intron_plus <- gtf_df_exon %>%
    filter(strand=='+') %>%
    mutate(intron_start=end, type="intron")
  gtf_df_intron_plus$intron_end <- shift(gtf_df_intron_plus$start, 1)
  print(gtf_df_intron_plus)
  
  #c("seqnames", "start", "end", "width", "strand", "source", "type",  "score", "phase", "gene_id", "gene_name")
  
  #gtf_df_intron_plus %>% filter(gene_id=="YDR376W")
  
  gtf_df_intron_plus <- gtf_df_intron_plus %>%
    filter(group_id!=nb_exon)%>%
    group_by(!!as.name(gtf_col)) %>%
    arrange(!!as.name(gtf_col), group_id)
  print(gtf_df_intron_plus)
  
  print("### building intron_minus list from gtf")
  
  # strand -
  gtf_df_intron_minus <- gtf_df_exon %>%
    filter(strand=='-') %>% 
    arrange(!!as.name(gtf_col), start) %>%
    mutate(intron_start=end, type="intron")
  gtf_df_intron_minus$intron_end <- shift(gtf_df_intron_minus$start, 1)
  
  gtf_df_intron_minus <- gtf_df_intron_minus%>%
    filter(group_id!=nb_exon)%>%
    group_by(!!as.name(gtf_col)) %>%
    arrange(!!as.name(gtf_col), desc(group_id), type)
  #print(gtf_df_intron_minus %>% filter(gene_id=="YDR129C"))
  print("### merging introns + and -")
  
  # bind intron + and - and rename columns to avoid confusion
  gtf_df_intron <- rbind(gtf_df_intron_minus, gtf_df_intron_plus) %>%
    select(-c(start, end)) %>%
    dplyr::rename("start"="intron_start") %>%
    dplyr::rename("end"="intron_end")
  print(gtf_df_intron)
  print("### merging introns and exons")
  
  
  # now bind introns and exons and sort file so it makes sense
  gtf_df_intronexon <- rbind(gtf_df_exon, gtf_df_intron) %>%
    mutate(type_id=paste0(type, group_id))  %>%
    arrange(!!as.name(gtf_col), start) %>%
    mutate(reg_name=paste(!!as.name(gtf_col), type_id, sep='_'),
           point='.') %>%
    ungroup() %>%
    select(c(seqnames, start, end, reg_name, point, strand))
  #test <- gtf_df_intronexon %>% filter(gene_id=="YBL091C")
  
  print(gtf_df_intronexon)
  print("### outputing files")
  
  # write output files
  write_tsv(gtf_df_intronexon, file=out_ie, col_names = FALSE)
  print("...")
  write_tsv(gtf_df_intron%>%
              ungroup() %>%
              mutate(intron_id=paste0(!!as.name(gtf_col), '_', type, group_id)) %>%
              select(intron_id),
            file=out_si, col_names = FALSE)
  print("###DONE.")
  
}

print(out_ie)
print(out_si)


build_intronexon(gtf, gene_list$list, out_ie, out_si)

