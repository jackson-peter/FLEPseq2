library(tidyverse)
library(data.table)
library(rtracklayer)

##### ARGUMENTS #####
args = commandArgs(trailingOnly=TRUE)

# annotation file in gtf format
gtf <- args[1]
gtf_col <- args[2]

# subset list of genes (repr gene model or canonical list...)
gene_list_f <- args[3]
gene_list_col <- args[4]

print(gene_list_col)
cmd=paste0('grep -v "^#" ', gene_list_f)
gene_list <- as_tibble(fread(cmd, na.strings = ""))%>% 
  drop_na(gene_list_col) %>% 
  select({{gene_list_col}})

print(gene_list)

colnames(gene_list) <- c("list")

# output directory
outdir <- args[5]

# basename for gene lists files for flepseq
suffix <- args[6]

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
    filter((!!as.name(gtf_col)) %in% gene_list) %>%
    #    filter(gtf_col %in% gene_list) %>%
    group_by(!!as.name(gtf_col), type) %>%
    mutate(group_id = case_when(strand=="+" ~ row_number(),
                                strand=='-' ~ rev(row_number())),
           nb_exon=max(group_id),
           gene_name = str_replace(!!as.name(gtf_col), "_", "-")) # replacing underscores in gene id to avoid problems during execution
  
  print(gtf_df_exon)
  
  print("### building intron_plus list from gtf")
  
  # create introns entries using nb exons and strand
  # strand +
  gtf_df_intron_plus <- gtf_df_exon %>%
    filter(strand=='+') %>%
    mutate(intron_start=end, type="intron")
  gtf_df_intron_plus$intron_end <- shift(gtf_df_intron_plus$start, 1)
  print(gtf_df_intron_plus)
  
  gtf_df_intron_plus <- gtf_df_intron_plus %>%
    filter(group_id!=nb_exon)%>%
    group_by(gene_name) %>%
    arrange(gene_name, group_id)
  print(gtf_df_intron_plus)
  
  print("### building intron_minus list from gtf")
  
  # strand -
  gtf_df_intron_minus <- gtf_df_exon %>%
    filter(strand=='-') %>%
    arrange(gene_name, start) %>%
    mutate(intron_start=end, type="intron")
  gtf_df_intron_minus$intron_end <- shift(gtf_df_intron_minus$start, 1)
  
  gtf_df_intron_minus <- gtf_df_intron_minus%>%
    filter(group_id!=nb_exon) %>%
    filter(group_id!=0)%>%
    group_by(gene_name) %>%
    arrange(gene_name, desc(group_id), type)
  print(gtf_df_intron_minus)
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
    arrange(gene_name, start) %>%
    mutate(reg_name=paste(gene_name, type_id, sep='_'),
           point='.') %>%
    ungroup() %>%
    select(c(seqnames, start, end, reg_name, point, strand))
  print(gtf_df_intronexon)
  print("### outputing files")
  
  # write output files
  write_tsv(gtf_df_intronexon, file=out_ie, col_names = FALSE)
  print("...")
  write_tsv(gtf_df_intron%>%
              ungroup() %>%
              mutate(intron_id=paste0(gene_name, '_', type, group_id)) %>%
              select(intron_id),
            file=out_si, col_names = FALSE)
  print("###DONE.")
  
}




build_intronexon(gtf, gene_list$list, out_ie, out_si)
print(out_ie)
print(out_si)
