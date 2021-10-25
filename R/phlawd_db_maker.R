#!/usr/bin/env Rscript

library(tidyverse)
library(RSQLite)
library(rentrez)
library(XML)
library(dplyr)

# get the labridae database
r_search <- entrez_search(db = "nuccore", term = "labridae[ORGN] AND cytb[GENE]")
r_search <- entrez_search(db = "nucleotide", term = "labridae[Organism] AND cytb[Gene]", retmax = r_search$count)

# prepare for data fetching
RECORDS_PER_ITERATION <- 200
iterations <- (r_search$count / RECORDS_PER_ITERATION) + 1
labridae_df <- data.frame()

# fetch the data and wrangle the XML database
for (i in 1:iterations) {
  print(i)
  upload <- entrez_post(db = "nucleotide", r_search$ids[RECORDS_PER_ITERATION*(i-1)+1:RECORDS_PER_ITERATION*i])
  entrez_XML <- entrez_fetch(db = "nucleotide", rettype = "xml", retmode = "xml", 
                             web_history = upload, retmax = RECORDS_PER_ITERATION, parsed = TRUE)
  
  # extract the id
  print('extracting id...')
  id_df <- xmlToDataFrame(entrez_XML, nodes = getNodeSet(entrez_XML, "//GBSeqid"))
  id_df <- id_df[grep("gi", id_df$text),]
  id_df <- str_replace(id_df, "gi\\|", "")
  id_df <- as.numeric(id_df)
  id_df <- cbind(id_df, id_df)
  
  # extract the title
  print('extracting title...')
  title_df <- xmlToDataFrame(entrez_XML, nodes = getNodeSet(entrez_XML, "//GBReference"))
  title_df <- filter(title_df, GBReference_reference == 1)
  title_df <- title_df[, c('GBReference_title')]
  
  # get the rest of the contents
  print('extracting contents...')
  content_df <- xmlToDataFrame(entrez_XML, nodes = getNodeSet(entrez_XML, "//GBSeq"))
  
  # combine three dataframes at once
  print('combining...')
  tmp_df <- cbind(id_df, cbind(content_df, title_df))
  labridae_df <- rbind(labridae_df, tmp_df)
}

# select and rename columns according to the SQL prompt
names(labridae_df)[1] <- "id"
names(labridae_df)[2] <- "ncbi_id"

labridae_filtered <- labridae_df[, c('id', 'ncbi_id',
                                     'GBSeq_locus', 'GBSeq_primary-accession', 
                                     'GBSeq_accession-version', 'GBSeq_definition',
                                     'title_df', 'GBSeq_sequence')]

names(labridae_filtered)[3] <- 'locus'
names(labridae_filtered)[4] <- 'accession_id'
names(labridae_filtered)[5] <- 'version_id'
names(labridae_filtered)[6] <- 'description'
names(labridae_filtered)[7] <- 'title'
names(labridae_filtered)[8] <- 'seq'

# SQL Part
mydb <- dbConnect(RSQLite::SQLite(), "")
dbWriteTable(mydb, "sequences", labridae_filtered, overwrite = TRUE)
