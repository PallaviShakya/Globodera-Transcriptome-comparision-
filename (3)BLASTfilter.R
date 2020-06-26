          ####best reciprocal blast hits
		  
		  ##1. loading the reciprocal blast results 
                rosto_on_pal <- read.delim("F:/Internship/Reciprocal BLAST results/qrosto_onpallida.blastn", header = F)
                    colnames(rosto_on_pal) <- c("rostoid","qacc","pal_hit","evalue","bitscore","length","pident","qstart","qend","sstart","send"); head(rosto_on_pal)
                
				pal_on_rosto <- read.delim("F:/Internship/Reciprocal BLAST results/qpal_on_rosto.blastn", header = F)
                    colnames(pal_on_rosto) <- c("palid","qacc","rosto_hit","evalue","bitscore","length","pident","qstart","qend","sstart","send"); head(pal_on_rosto)

		##2. Filtering the result by length_hit>100 and grouping under pal_id, rosto_hit, rosto_id and pal_hit 
                    library(tidyverse)
                    rosto_on_pal2 <- group_by(rosto_on_pal, rostoid, pal_hit) %>%
                                     summarise(length_hit=sum(length), maxhit=max(evalue)) %>%
                                      
                                     filter(length_hit>100) %>%
                                     data.frame() %>%
                                     mutate(rostoid=as.character(unlist(rostoid)), pal_hit=as.character(unlist(pal_hit)))
                
                
                    
                    pal_on_rosto2 <- group_by(pal_on_rosto, palid, rosto_hit) %>%
                                     summarise(length_hit=sum(length), maxhit=max(evalue)) %>%
                                     filter(length_hit>100) %>%
                                     data.frame() %>%
                                     mutate(palid=as.character(unlist(palid)), rosto_hit=as.character(unlist(rosto_hit)))
                
        ##3. used python to filter again by evalue<0.001 and merge (reason: low memory efficiency of R) See section 3.1
               
                      
                             #test <- merge(pal_on_rosto2,rosto_on_pal22,by.x=2,by.y=1)
                    
                    
         ##4. filtered best hits based on evalue             
                    oneonone2 <- filter(pal_palid_merged, palid==pal_hit) %>%
								 separate(palid, into = c("palid", "bla"), sep = "\\.") %>%
								 dplyr::select(-bla) %>%
                                 group_by(palid) %>%
                                 mutate(palbest= (maxhit_pal_on_rosto == min(maxhit_pal_on_rosto))) %>%
                                 group_by(rosto_hit) %>%
                                 mutate(rostobest= (maxhit_rosto_n_pal  == min(maxhit_rosto_n_pal))) %>%
                                 data.frame() %>%
                                 filter(palbest, rostobest, !duplicated(palid))


									
                      write.csv(oneonone2, 'F:/Internship/outputR/oneonone.csv', row.names = TRUE)
                  
                                
                    