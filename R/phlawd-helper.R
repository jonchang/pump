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
