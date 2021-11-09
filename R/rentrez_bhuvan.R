install.packages("rentrez")
install.packages("glue")
library(rentrez)
library(glue)
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
x = 200
upload <- c()
first <- list()
for(j in 1:length(r_search2)){
  iterations = (r_search2[[j]]$count / x) + 1
  for (i in 1:iterations) {
    print(i)
    upload <- entrez_post(db="nucleotide", r_search2[[j]]$ids[x*(i-1)+1:x*i])
    first[[i]] <- entrez_fetch(db="nucleotide", rettype="xml", retmode = "xml", web_history=upload,
                               retmax=x)
  }
}
r_search2[[1]]