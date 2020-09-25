###this script should give me the matrix that has gene_id and TPM 

library(data.table)
library(dplyr)
library(tidyverse)


# In the datasets: Attributes column is a collection of bits of information we are interested in 
#in order to pull out the information I want based on gene_id and TPM, 
#this is the function I use: 

extract_attributes <- function(gtf_attributes, att_of_interest){
  att <- unlist(strsplit(gtf_attributes, " "))
  if(att_of_interest %in% att){
    return(gsub("\"|;","", att[which(att %in% att_of_interest)+1]))
  }else{
    return(NA)}
}


#### I need a loop to go to the folder and read the gtf files and do these for each file: 

filenames <- list.files(path ="F:/Internship/workingdirectoryR", 
                        pattern="*.gtf")
#sample <- basename(filenames)

combined_matrix <- matrix(ncol = 3, nrow = 0) # empty matrix to store the loop results in

  for (i in 1:length(filenames)){
    file <- fread(filenames[i]) %>%
      dplyr::filter(V3=='transcript') 
      file<- setnames(file, names(file), c("chr","source","type","start","end","score","strand","phase","attributes"))
      file$gene_id <- unlist(lapply(file$attributes, extract_attributes, "gene_id"))
      file$TPM <- unlist(lapply(file$attributes, extract_attributes, "TPM"))
      file<- dplyr::mutate(file, Sample= gsub(".gtf","",filenames[i]))
      combined_matrix <- rbind(combined_matrix, file, fill=TRUE)
     }

  
          filtered_data <- combined_matrix %>% 
          select(Sample, gene_id, TPM) %>%
            mutate(j = row_number()) %>%
            spread(Sample, TPM)%>%
            select(-j)
          
          write.csv(filtered_data, file = '/mnt/nemahomes/shaky005/RNAseq_data/stringtie_rosto/filtereddata_rosto.csv')
            
          
          
  




