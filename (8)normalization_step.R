aggregated_geneid_pal <- read.csv("aggregated_geneid_pal - Copy.csv")
aggregated_geneid_pal$`Unnamed: 0`<- NULL 
head(aggregated_geneid_pal)


ag_geneid <- as.data.frame(aggregated_geneid_pal)

##to remove .gtf from sample names 
for (col in 1:ncol(ag_geneid)){ 
  colnames(ag_geneid)[col] <- sub(".gtf", "", colnames(ag_geneid)[col])
}
ag_geneid$`Unnamed: 0`<- NULL 
head(ag_geneid)

########################normalize time series data####################################################
######divide each row with its mean####################################################################
ag_geneid <- cbind(ag_geneid[1], ag_geneid[-1]/rowMeans(ag_geneid[-1]))

##########################################log2 transormation##########################

isnum <- sapply(ag_geneid, is.numeric)
ag_geneid[,isnum] <- lapply(ag_geneid[,isnum], log2)

############################replacing infinity and NA with 0#########################################
ag_geneid <- do.call(data.frame, lapply(ag_geneid, function(x){
  replace(x, is.infinite(x)| is.na(x), 0)
})
)

write.csv(ag_geneid, "actual_normalized_pal.csv")
