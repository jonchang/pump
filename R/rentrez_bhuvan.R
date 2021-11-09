install.packages("rentrez")
install.packages("glue")
library(rentrez)
library(glue)

# maybe: change entrezXML to record
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
  
  # combine three dataframes at once
  print('combining...')
  tmp_df <- cbind(id_df, cbind(content_df, title_df))
  labridae_df <- rbind(labridae_df, tmp_df)
}

gene_names <- c("12s", "16s", "4c4", "coi", "cytb", "enc1", "ficd", "glyt", "hoxc6a", "kiaa1239", "myh6", "panx2", "plagl2", "ptr", "rag1", "rag2", "rhodopsin", "ripk4", "sh3px3", "sidkey", "sreb2", "svep1", "tbr1", "vcpip", "zic1")
r_search1 <- list()
r_search2 <- list()
search_term_frame = "labridae[Organism] AND "
for(i in 1: length(gene_names))
{
  r_glue <-  glue("{search_term_frame}{gene_names[i]}[Gene]")
  r_glue2 <-  glue("{search_term_frame}{gene_names[i]}")
  print(r_glue)
  r_search1[[i]] <- entrez_search(db="nucleotide", term = r_glue)
  r_search2[[i]] <- entrez_search(db="nucleotide", term = r_glue, retmax = r_search1[[i]]$count)
  r_search1
  if(r_search1[[i]]$count == 0)
  {
    print("search had 0 results, researching without gene field")
    r_search1[[i]] <- entrez_search(db="nucleotide", term = r_glue2)
    r_search2[[i]] <- entrez_search(db="nucleotide", term = r_glue2, retmax = r_search1[[i]]$count)
    print(r_search2[[i]]$count)
  }
  else
  {
    print(r_search2[[i]]$count)
  }
}
x = 1
upload <- c()
first <- list()
for(j in 1:1){
  iterations = (r_search2[[j]]$count / x) + 1
  for (i in 1:5) {
    print(i)
    upload <- entrez_post(db="nucleotide", r_search2[[j]]$ids[x*(i-1)+1:x*i])
    curr_record <- entrez_fetch(db="nucleotide", rettype="xml", retmode = "xml", web_history=upload,
                               retmax=x)
    sqlAdd(curr_record)
    # catch any invalid records
    # grab those records and store it into a logfile
  }
}

# rename labridae_df into something more general, like fish_df
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
