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

sqlAdd <- function(entrez_XML) {
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
  tax_df <- content_df[, c('GBSeq_organism', 'GBSeq_taxonomy')]
  
  # combine three dataframes at once
  print('combining...')
  mid_df <- cbind_with_na(content_df, title_df, "GBSeq_references", "title")
  labridae_df <- cbind(id_df, mid_df)
  tax_df <- cbind(id_df, tax_df)
  
  # filter database after including everything
  print('filtering...')
  names(labridae_df)[1] <- "id"
  names(labridae_df)[2] <- "ncbi_id"
  
  batch_filtered <- labridae_df[, c('id', 'ncbi_id',
                                    'GBSeq_locus', 'GBSeq_primary-accession', 
                                    'GBSeq_accession-version', 'GBSeq_definition',
                                    'title', 'GBSeq_sequence')]
  
  names(batch_filtered)[3] <- 'locus'
  names(batch_filtered)[4] <- 'accession_id'
  names(batch_filtered)[5] <- 'version_id'
  names(batch_filtered)[6] <- 'description'
  names(batch_filtered)[7] <- 'title'
  names(batch_filtered)[8] <- 'seq'
  
  names(tax_df)[1] <- "id"
  names(tax_df)[2] <- "ncbi_id"
  names(tax_df)[3] <- "name"
  names(tax_df)[4] <- "taxonomy"
  # names(tax_df)[4] <- "name_class"
  # names(tax_df)[5] <- "node_rank"
  # names(tax_df)[6] <- "parent_ncbi_id"
  # names(tax_df)[7] <- "edited_name"
  # names(tax_df)[8] <- "left_value"
  # names(tax_df)[9] <- "right_value"
  
  # add data to SQL table
  print('adding data to SQL...')
  dbWriteTable(mydb, "sequences", batch_filtered, append = TRUE)
  dbWriteTable(taxdb, "taxonomy", tax_df, append = TRUE)
}

# search for fish entries
# anything with vrt?
gene_names <- c("12s", "16s", "4c4", "coi", "cytb", "enc1", "ficd", "glyt", "hoxc6a", "kiaa1239", "myh6", "panx2", "plagl2", "ptr", "rag1", "rag2", "rhodopsin", "ripk4", "sh3px3", "sidkey", "sreb2", "svep1", "tbr1", "vcpip", "zic1")
r_search1 <- list()
r_search2 <- list()

# search for actinopterygii, not labridae
search_term_frame = "labridae[Organism] AND "


for(i in 1:length(gene_names))
{
  # append search terms
  r_glue <-  glue("{search_term_frame}{gene_names[i]}[Gene]")
  r_glue2 <-  glue("{search_term_frame}{gene_names[i]}")
  print(r_glue)
  
  # rentrez search nucleotides
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
  
  # rentrez search taxonomy
  r_search_tax1[[i]] <- entrez_search(db="nucleotide", term = r_glue)
  r_search_tax2[[i]] <- entrez_search(db="nucleotide", term = r_glue, retmax = r_search1[[i]]$count)
  
  if(r_search_tax1[[i]]$count == 0)
  {
    print("taxonomy search had 0 results, researching without gene field")
    r_search_tax1[[i]] <- entrez_search(db="nucleotide", term = r_glue2)
    r_search_tax2[[i]] <- entrez_search(db="nucleotide", term = r_glue2, retmax = r_search1[[i]]$count)
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
  iterations = (r_search2[[j]]$count / RECORDS_PER_ITERATION) + 1
  
  # real data is 1:iterations
  for (i in 1:iterations) {
    print(i)
    upload <- entrez_post(db="nucleotide", r_search2[[j]]$ids[RECORDS_PER_ITERATION*(i-1)+1:RECORDS_PER_ITERATION*i])
    curr_record <- entrez_fetch(db="nucleotide", rettype="xml", retmode = "xml", web_history=upload,
                               retmax=RECORDS_PER_ITERATION, parsed = TRUE)
    sqlAdd(curr_record)
    
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