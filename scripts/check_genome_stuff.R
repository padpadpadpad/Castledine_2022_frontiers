# play with PacBio sequencing data

# load in packages
library(tidyverse)
library(Biostrings)
library(rBLAST)
library(ShortRead)
library(vcfR)
library(flextable)

d4 <- readFastq('genome_seq/cleaned_sequences/file4_corrected_reads.fastq.gz') %>%
  sread(.)
length(d4)

# load reference in - rBlast only works on reference files that do not have spaces in the file name
# move the folder if this is the case (as it is for me on Google Drive)

# now to try blast ####
makeblastdb('~/Desktop/reference/GCA_000886135.1_ViralProj42717_genomic.fna', dbtype = 'nucl')
db <- blast('~/Desktop/reference/GCA_000886135.1_ViralProj42717_genomic.fna')

d4_prediction <- predict(db, d4)
unique(d4_prediction$QueryID) %>% length()

#---------------------------------------#
# read in gff file and look for SNPs ####
#---------------------------------------#

tidy_freebayes <- function(freebayes_vcf){
  temp <- vcfR::read.vcfR(freebayes_vcf, verbose = FALSE)
  temp <- vcfR::vcfR2tidy(temp, single_frame = TRUE) %>%
    .$dat %>%
    janitor::clean_names() %>%
    dplyr::filter(dp > mean(dp) - 2*sd(dp) & dp < mean(dp) + 2*sd(dp)) %>%
    dplyr::select(., -c(id, filter, ns))
  temp$file <- tools::file_path_sans_ext(basename(freebayes_vcf))
  return(temp)
}

# function to assign a SNP pos to its position in the genome
# returns the downstream, current, and upstream gene, and the distance to those
get_gene_info <- function(SNP_pos, d_genes, types_to_keep = c('gene'), pos_to_keep = c('downstream', 'current_pos', 'upstream')){
  
  # make numeric just in case it is not - needed for use of tidyr and nest()
  SNP_pos <- as.numeric(SNP_pos)
  
  # keep just gene names in d_genes
  d_genes <- dplyr::filter(d_genes, type %in% types_to_keep)
  
  # get gene name
  temp <- dplyr::filter(d_genes, start <= SNP_pos & stop >= SNP_pos) %>%
    dplyr::pull(., name)
  
  # if it is not there then return NA
  if(length(temp) == 0) temp <- NA
  
  # if is not NA, delete gene name from d_genes
  if(length(temp) > 0){d_genes <- dplyr::filter(d_genes, ! name %in% temp)}
  
  # some repeat regions overlap so just put return the first output of this
  # I dont think any of the gene positions overlap
  if(length(temp > 1)) temp <- temp[1]
  
  # get downstream gene - later position
  d_down <- dplyr::filter(d_genes, start >= SNP_pos) %>%
    arrange(., start) %>%
    dplyr::slice(seq_len(1))
  if(nrow(d_down) == 0){
    d_down <- data.frame(name = NA,
                         start = SNP_pos)
  }
  
  # get upstream_gene - earlier position
  d_up <- dplyr::filter(d_genes, stop <= SNP_pos) %>%
    arrange(., desc(start)) %>%
    dplyr::slice(seq_len(1))
  if(nrow(d_up) == 0){
    d_up <- data.frame(name = NA,
                       stop = SNP_pos)
  }
  
  # get upstream pos
  d_temp <- data.frame(gene_type = c('current_pos', 'downstream', 'upstream'),
                       gene_name = c(temp, d_down$name, d_up$name),
                       distance = c(0, abs(d_down$start - SNP_pos), abs(d_up$stop - SNP_pos)))
  
  # filter based on what to keep
  d_temp <- dplyr::filter(d_temp, gene_type %in% pos_to_keep)
  
  return(d_temp)
  
}

# function to clean string up
clean_gff_note <- function(note){
  # lowercase
  temp <- tolower(note)
  # remove things that arent a letter or number
  temp <- stringr::str_replace_all(temp,"[^[:alnum:]]", " ")
  # change white space to '_'
  temp <- stringr::str_replace_all(temp,"[\\s]+", "_")
  temp <- if(str_split(temp, pattern = '')[[1]][nchar(temp)] == '_'){temp <- substr(temp, 1, nchar(temp) - 1)}
  return(temp)
}

# read in gff file ####
gff <- readr::read_tsv('genome_seq/reference/GCA_000886135.1_ViralProj42717_genomic.gff',
                col_names = c(
                  "seqid",
                  "source",
                  "type",
                  "start",
                  "stop",
                  "score",
                  "strand",
                  "phase",
                  "attr"
                ),
                na        = ".",
                comment   = "#",
                col_types = "ccciidcic"
)

# clean up gff file ####

# keep protein coding
gff1 <- # remove first line (type == region) as it described the whole genome
  filter(gff, type != 'region' & type != 'rRNA' & type != 'CDS') %>%
  # some funky regex
  mutate(., attr = tolower(attr),
         attr_id = gsub(pattern = "(.*id=)(.*?)(;.*)", replacement = "\\2", attr),
         name = gsub(pattern = "(.*name=)(.*?)(;.*)", replacement = "\\2", attr),
         gene_biotype = gsub(pattern = "(.*gene_biotype=)(.*?)(;.*)", replacement = "\\2", attr),
         note = gsub(pattern = "(.*note=)(.*?)(;.*)", replacement = "\\2", attr),
         # make instances where no name is present the type by id
         name = ifelse(grepl('=', name), paste(attr_id, type, sep = '_'), name),
         gene_biotype = ifelse(grepl('=', gene_biotype), NA, gene_biotype),
         note = ifelse(grepl('=', note), NA, note),
         note = clean_gff_note(note)) %>%
  select(., seqid, type, start, stop, strand, attr_id, name, gene_biotype, note)

gff2 <- # remove first line (type == region) as it described the whole genome
  filter(gff, type != 'region' & type != 'rRNA' & type == 'CDS') %>%
  # some funky regex
  mutate(., attr = tolower(attr),
         attr_id = gsub(pattern = "(.*id=)(.*?)(;.*)", replacement = "\\2", attr),
         name = gsub(pattern = "(.*name=)(.*?)(;.*)", replacement = "\\2", attr),
         gene_biotype = gsub(pattern = "(.*gene_biotype=)(.*?)(;.*)", replacement = "\\2", attr),
         note2 = gsub(pattern = "(.*note=)(.*?)(;.*)", replacement = "\\2", attr),
         # make instances where no name is present the type by id
         name = ifelse(grepl('=', name), paste(attr_id, type, sep = '_'), name),
         gene_biotype = ifelse(grepl('=', gene_biotype), NA, gene_biotype),
         note2 = ifelse(grepl('=', note2), NA, note2),
         note2 = clean_gff_note(note2)) %>%
  select(., seqid, type, start, stop, strand, attr_id, name, gene_biotype, note2)

# bind both gff outputs
gff <- merge(gff1, select(gff2, start, note2), by = 'start', all.x = TRUE) %>%
  mutate(note = ifelse(is.na(note), note2, note))%>%
  select(., seqid, type, everything())

# add a column to easily identify tail fiber gene
gff <- mutate(gff, tail_gene = ifelse(grepl('tail', note), 'yes', 'no'))

# read in SNP/indel file
d_snp <- tidy_freebayes('genome_seq/vcf/file4_filter.vcf') %>%
  select(., chrom:pao)

d_snp <- d_snp %>%
  mutate(., n = 1:n()) %>%
  nest(-n) %>%
  mutate(gene_info = purrr::map(data, ~get_gene_info(.x$pos, d_genes = gff))) %>%
  unnest(., data) %>%
  unnest(., gene_info) %>%
  filter(., gene_type == 'current_pos') %>%
  left_join(., select(gff, gene_name = name, tail_gene))

d_snp

# read in ancestral 
# read in SNP/indel file
d_snp2 <- tidy_freebayes('genome_seq/vcf/110-Pf_Phi2_ancestral_171011_L008_minimap2_sorted_freebayes_filter.vcf') %>%
  select(., chrom:pao)

d_snp2 <- d_snp2 %>%
  mutate(., n = 1:n()) %>%
  nest(-n) %>%
  mutate(gene_info = purrr::map(data, ~get_gene_info(.x$pos, d_genes = gff))) %>%
  unnest(., data) %>%
  unnest(., gene_info) %>%
  filter(., gene_type == 'current_pos') %>%
  left_join(., select(gff, gene_name = name, tail_gene))

# lets make a table
table_info_phage1 <- select(d_snp2, pos, ref, alt, tail_gene) %>%
  mutate(phage = 'Phage 1')
table_info_phage2 <- select(d_snp, pos, ref, alt, tail_gene) %>%
  mutate(phage = 'Phage 2') 

d_table <- bind_rows(table_info_phage1, table_info_phage2) %>%
  mutate(in_gene = ifelse(is.na(tail_gene), 'no', 'yes')) %>%
  select(phage, pos, ref, alt, in_gene, tail_gene)

table <- flextable(d_table) %>%
  set_header_labels(phage = 'Phage\nvariant',
                    pos = 'Genome\nposition',
                    ref = 'Reference\nsequence',
                    alt = 'Alternative\nvariant',
                    in_gene = 'Variant in\na putative gene',
                    tail_gene = 'Variant in\na tail fiber gene') %>%
  font(fontname = 'Times', part = 'all') %>% 
  autofit() %>%
  align(align = 'center', part = 'all') %>%
  merge_v(j = 'phage') %>%
  valign(valign = 'top', part = 'body') %>%
  hline(i = c(3), border = fp_border_default()) %>%
  fix_border_issues() %>%
  bold(~pos == 37265, j = c(2:6))
table

save_as_image(table, 'plots/phi2_differences.png', webshot = 'webshot2')
