# install.packages("rentrez")
# install.packages("glue")
# install.packages("RSQLite")

library(tidyverse)
library(rentrez)
library(glue)
library(XML)
library(dplyr)
library(RSQLite)

# change these variables so it fits with your computer
pump_dir <- "./UCLA/Research/pump"   # filepath to the pump repository


# functions
# join two columns while keeping NA entries
cbind_with_na <- function(df1, df2, colname, new_name) {
  i <- 1
  j <- 1
  
  df1[new_name] <- NA
  col <- df1[colname][[1]]
  col2 <- df1[new_name][[1]]
  
  while (i <= length(col)) {
    if (!is.na(col[i])) {
      col2[i] <- df2[j]
      j <- j + 1
    }
    i <- i + 1
  }
  
  df1[new_name] <- col2
  return(df1)
}

# get taxonomy ids to get values from the taxonomy database
get_tax_ids <- function(sequence_XML) {
  qual_df <- xmlToDataFrame(sequence_XML, nodes = getNodeSet(sequence_XML, "//GBSeq_feature-table/GBFeature/*/GBQualifier"))
  qual_df <- filter(qual_df, GBQualifier_name == "db_xref")
  qual_df <- qual_df[grepl('taxon:', qual_df$GBQualifier_value),]
  qual_df <- qual_df[, c('GBQualifier_value')]
  qual_df <- str_replace(qual_df, "taxon:", "")
  qual_df <- na.omit(as.numeric(qual_df))
  
  return(as.list(qual_df))
}

# build the two SQL tables
sqlAdd <- function(sequence_XML, taxonomy_XML) {
  # extract the id
  print('extracting contents...')
  id_df <- xmlToDataFrame(sequence_XML, nodes = getNodeSet(sequence_XML, "//GBSeqid"))
  id_df <- id_df[grep("gi", id_df$text),]
  id_df <- str_replace(id_df, "gi\\|", "")
  id_df <- as.numeric(id_df)
  id_df <- cbind(id_df, id_df)
  
  # extract the title
  title_df <- xmlToDataFrame(sequence_XML, nodes = getNodeSet(sequence_XML, "//GBReference"))
  title_df <- filter(title_df, GBReference_reference == 1)
  title_df <- title_df[, c('GBReference_title')]
  
  # get the rest of the contents
  content_df <- xmlToDataFrame(sequence_XML, nodes = getNodeSet(sequence_XML, "//GBSeq"))
  
  # combine three dataframes at once
  print('combining...')
  mid_df <- cbind_with_na(content_df, title_df, "GBSeq_references", "title")
  labridae_df <- cbind(id_df, mid_df)
  
  # repeat for taxonomy data
  print('repeating for taxonomy...')
  taxon_df <- xmlToDataFrame(tax_record, nodes = getNodeSet(tax_record, "//Taxon"))
  taxon_df <- taxon_df[!(is.na(taxon_df$ParentTaxId)), ] # get only rows that have full contents
  rownames(taxon_df) <- NULL
  taxon_df <- cbind(id_df, taxon_df)
  
  # filter database after including everything
  print('finalizing sql addition...')
  names(labridae_df)[1] <- "id"
  names(labridae_df)[2] <- "ncbi_id"
  
  seq_df <- labridae_df[, c('id', 'ncbi_id',
                            'GBSeq_locus', 'GBSeq_primary-accession', 
                            'GBSeq_accession-version', 'GBSeq_definition',
                            'title', 'GBSeq_sequence')]
  
  names(seq_df)[3] <- 'locus'
  names(seq_df)[4] <- 'accession_id'
  names(seq_df)[5] <- 'version_id'
  names(seq_df)[6] <- 'description'
  names(seq_df)[7] <- 'title'
  names(seq_df)[8] <- 'seq'
  
  names(taxon_df)[1] <- "id"
  names(taxon_df)[2] <- "ncbi_id"
  
  taxon_df["name_class"] <- "scientific name"
  taxon_df["left_value"] <- 0
  taxon_df["right_value"] <- 0
  
  tax_df <- taxon_df[, c('id', 'ncbi_id', 'ScientificName', 
                         'name_class', 'Rank', 'ParentTaxId',
                         'ScientificName', 'left_value', 'right_value')]
  
  names(tax_df)[3] <- "name"
  names(tax_df)[4] <- "name_class"
  names(tax_df)[5] <- "node_rank"
  names(tax_df)[6] <- "parent_ncbi_id"
  names(tax_df)[7] <- "edited_name"
  
  # add data to SQL table
  # print('adding data to SQL...')
  dbWriteTable(mydb, "sequences", seq_df, append = TRUE)
  dbWriteTable(taxdb, "taxonomy", tax_df, append = TRUE)
}

# search for fish entries
# anything with vrt?
gene_names <- c("12s", "16s", "4c4", "coi", "cytb", "enc1", "ficd", "glyt", "hoxc6a", "kiaa1239", "myh6", "panx2", "plagl2", "ptr", "rag1", "rag2", "rhodopsin", "ripk4", "sh3px3", "sidkey", "sreb2", "svep1", "tbr1", "vcpip", "zic1")
r_search1 <- list()
r_search2 <- list()

# search for actinopterygii, not labridae
search_term_frame = "labridae[Organism] AND "

# fill in search logs
for(i in 1:length(gene_names))
{
  # append search terms
  r_glue <-  glue("{search_term_frame}{gene_names[i]}[Gene]")
  r_glue2 <-  glue("{search_term_frame}{gene_names[i]}")
  print(r_glue)
  
  # rentrez search
  r_search1[[i]] <- entrez_search(db="nucleotide", term = r_glue)
  r_search2[[i]] <- entrez_search(db="nucleotide", term = r_glue, retmax = r_search1[[i]]$count)
  
  if(r_search1[[i]]$count == 0)
  {
    print("nucleotide search had 0 results, researching without gene field")
    r_search1[[i]] <- entrez_search(db="nucleotide", term = r_glue2)
    r_search2[[i]] <- entrez_search(db="nucleotide", term = r_glue2, retmax = r_search1[[i]]$count)
    print(r_search2[[i]]$count)
  }
  else
  {
    print(r_search2[[i]]$count)
  }
}

# Set number of records per iteration
RECORDS_PER_ITERATION <- 200

# Initialize SQL database
mydb <- dbConnect(RSQLite::SQLite(), "")
taxdb <- dbConnect(RSQLite::SQLite(), "")

# Fetch and add fish entries to SQL
# real test is 1:length(r_search2)
for(j in 1:length(r_search2)){
  print(paste("fetching ", gene_names[j]))
  iterations = as.integer((r_search2[[j]]$count / RECORDS_PER_ITERATION) + 1)
  
  # real data is 1:iterations
  for (i in 1:iterations) {
    print(paste(i, "of", iterations, "iterations"))
    
    # fetch records for sequence data
    print("fetching sequence data...")
    upload <- entrez_post(db="nucleotide", r_search2[[j]]$ids[RECORDS_PER_ITERATION*(i-1)+1:RECORDS_PER_ITERATION*i])
    curr_record <- entrez_fetch(db="nucleotide", rettype="xml", retmode = "xml", web_history=upload,
                                retmax=RECORDS_PER_ITERATION, parsed = TRUE)
    
    # fetch records for taxonomy data
    print("obtaining taxonomy ids...")
    tax_ids <- get_tax_ids(curr_record)
    
    print("fetching taxonomy data...")
    tax_upload <- entrez_post(db="taxonomy", tax_ids)
    tax_record <- entrez_fetch(db="taxonomy", rettype="xml", retmode = "xml", web_history=tax_upload,
                               retmax=RECORDS_PER_ITERATION, parsed = TRUE)
    
    # build sql tables
    sqlAdd(curr_record, tax_record)
    print("\n")
    
    # TODO: catch any invalid records
    # grab those records and store it into a logfile
    # apparently, we haven't found any invalid records yet.
  }
}

# Extract sql information
data_sql <- dbGetQuery(mydb, 'SELECT * FROM sequences')
tax_sql <- dbGetQuery(taxdb, 'SELECT * FROM taxonomy')

# Export sql
curr_dir <- getwd()
setwd(pump_dir)
write.csv(data_sql,"sequences.csv", row.names = FALSE)
write.csv(tax_sql,"taxonomy.csv", row.names = FALSE)
setwd(curr_dir)

# phlawd_db_maker vrt vrt.db

# DEBUG
# upload <- entrez_post(db="nucleotide", r_search2[[1]]$ids[RECORDS_PER_ITERATION*(1-1)+1:RECORDS_PER_ITERATION*1])
# curr_record <- entrez_fetch(db="nucleotide", rettype="xml", retmode = "xml", web_history=upload,
#                             retmax=RECORDS_PER_ITERATION, parsed = TRUE)
# 
# content_df_2 <- xmlToDataFrame(curr_record, nodes = getNodeSet(curr_record, "//GBSeq"))
# curr_record
# 
# trip_df <- data.frame()
# trip_df <- xmlToDataFrame(curr_record, nodes = getNodeSet(curr_record, "//GBSeq_feature-table/GBFeature/*/GBQualifier"))
# trip_df <- filter(trip_df, GBQualifier_name == "db_xref")
# trip_df <- trip_df[!grepl('GeneID', trip_df$GBQualifier_value),]
# trip_df <- trip_df[, c('GBQualifier_value')]
# trip_df <- str_replace(trip_df, "taxon:", "")
# 
# trip_df <- na.omit(as.numeric(trip_df))
# 
# tax_ids <- as.list(trip_df)
# 
# tax_upload <- entrez_post(db="taxonomy", tax_ids[RECORDS_PER_ITERATION*(1-1)+1:RECORDS_PER_ITERATION*1])
# tax_record <- entrez_fetch(db="taxonomy", rettype="xml", retmode = "xml", web_history=tax_upload,
#                             retmax=RECORDS_PER_ITERATION, parsed = TRUE)
# 
# tax_df <- xmlToDataFrame(tax_record, nodes = getNodeSet(tax_record, "//Taxon"))
# tax_df <- tax_df[!(is.na(tax_df$ParentTaxId)), ]
# rownames(tax_df) <- NULL
# 
# taxes_df <- tax_df[, c('ScientificName', 'ScientificName')]
# taxes_df['class'] <- "scientific name"
# taxes_df['left'] <- 0
# 
# tax_record

j = 2
k = 3
paste(j, "iterations", k)
