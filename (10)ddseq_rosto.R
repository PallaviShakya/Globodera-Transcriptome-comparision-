##setting working directory

setwd("F:/Internship/outputR/Rosto_Results")


####acquiring needed packages 
install.packages("factoextra")
install.packages("pheatmap")
if (!requireNamespace('BiocManager', quietly = TRUE))
  install.packages('BiocManager')
BiocManager::install("DEGreport")
BiocManager::install('EnhancedVolcano')

library(DESeq2)
library(factoextra) 
library(pheatmap)
library(RColorBrewer)
library(EnhancedVolcano)
library(DEGreport)

####getting transcript count matrix made by using prepDE.py after stringtie step

    exp_data <- read.csv("transcript_count_matrix_rosto.csv", row.names="transcript_id")
    exp_data <- as.matrix(exp_data)
    dim(exp_data) ##21037 genes
    apply(exp_data, 2, sum)

#####cleaning the data 
#####removing genes that are not informative 
#####genes that have in none of the samples more than 0 counts 

    mx = apply(exp_data, 1, max) ## will create a list containing maximum read count(ver 65 samples)

    ##for each gene creating a list containing rows for which the maximum count is greater than 0: 

      exp_data <- exp_data[ mx>0, ]

      dim(exp_data) ##20343 

#########making col_data that contains stages of the samples as factors 
          a <- read.csv("stagesrosto.csv")
          a <- as.factor(a$Stages)
          levels(a)
          stages<- factor(c("females", "females", "na", "na", "dry cysts", "hydrated cysts",
                            "hydrated cysts(15m PRD)", "hydrated cysts(1h PRD)", "hydrated cysts(8h PRD)",
                            "hydrated cysts(24h PRD)","hydrated cysts(48h PRD)", "hydrated cysts(7d PRD)", 
                            "J2", "dry cysts", "J2", "J2", "J2", "J2", "J2", 
                            "J2", "J2", "J2", "J2", "J2"), c("dry cysts", "females", "J2", "na", "hydrated cysts", 
                                                             "hydrated cysts(15m PRD)", "hydrated cysts(1h PRD)", "hydrated cysts(8h PRD)",
                                                             "hydrated cysts(24h PRD)","hydrated cysts(48h PRD)", "hydrated cysts(7d PRD)" ))

          col_data <- data.frame(stages)


##########DESeq2 
          dds = DESeqDataSetFromMatrix(exp_data, col_data, ~stages)

######## Normalization Step 

        #calculate the linear correction factors for each sample:
        dds <- estimateSizeFactors(dds)
        sizeFactors(dds)## gives the correction factors for each sample 
        normalized_dds <- counts(dds, normalized = TRUE)
        
        write.csv(normalized_dds, file="F:/Internship/outputR/Rosto_Results/normalized_counts_pal.csv", quote=F, col.names=NA)


######Hierarchical Clustering 
        rld	= vst(dds) #transform the counts to get a more normal distribution, 
        
#######avoids highly expressed genes dominating the results. 
#######THIS SHOWS THE EFFECT OF TRANSFORMATION 
        par( mfrow = c( 1, 2 ) )
        
        plot( log2( 1 + counts(dds, normalized=TRUE)[ , 1:2] ),
              col=rgb(0,0,0,.2), pch=16, cex=0.3 )
        plot( assay(rld)[ , 1:2],
              col=rgb(0,0,0,.2), pch=16, cex=0.3 )
        
        dists = dist(t(assay(rld))) #calculating euclidean distances between the samples using their gene expression values
        save(dists, file = "hclust_rosto.RData")
        plot(hclust(dists))     

#######heatmap 
        rld_mat <- assay(rld)
        rld_cor <- cor(rld_mat)
        save(rld_cor, file = "correlation_rosto.RData")
        head(rld_cor)
        pheatmap(rld_cor)
        heat.colors <- brewer.pal(6, "Blues")
        pheatmap(rld_cor, color = heat.colors, border_color=NA, fontsize = 10, 
                 fontsize_row = 10, height=20)
        
        
#########PCA and PCO 

        plotPCA(rld, intgroup = "stages")
        
        #########PCO plot using the correlation matrix
        pco <- prcomp(t(rld_cor), scale = TRUE)
        save(pco, file = "pco_ros.RData")
        eigenvalue <- get_eigenvalue(pco)
        fviz_eig(pco, addlabels = TRUE, ylim = c(0, 100))
        var <- get_pca_var(pco)
        var
        #coordinates 
        head(var$coord)
        #quality on the factor map 
        head(var$cos2)
        #contributions to the principal components
        head(var$contrib)
        
        #correlation between a variable and a principal component is used as the 
        #coordinates of the variable on the PC. Variables are represented by correlations
        fviz_pca_var(pco, col.var = stages, palette = "default", repel = TRUE)
        fviz_cos2(pco, choice = "var", axes = 1:2, repel = TRUE)
        
        fviz_pca_ind(pco,
                     col.ind = "cos2", # Color by the quality of representation
                     gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                     repel = TRUE     # Avoid text overlapping
        )
        
        fviz_pca_ind(pco, 
                     geom = "point", 
                     col.ind = stages, 
                     palette = "default", 
                     addEllipses = TRUE, 
                     legend.title = "Stages"
        )
        
        
        #########PCA plot using the normalized expression matrix 
        pca <- prcomp(t(rld_mat), scale = TRUE)
        save(pca, file = "pca_rosto.RData")
        
        fviz_eig(pca) ##visualizing scree plot of eigen values (percentage of variances explained by each principal component)
        
        var <- get_pca_var(pca)
        var
        #coordinates 
        head(var$coord)
        #quality on the factor map 
        head(var$cos2)
        #contributions to the principal components
        head(var$contrib)
        
        #fviz_pca_var(pca, col.var = stages, palette = "default", repel = TRUE)
        #fviz_cos2(pca, choice = "var", axes = 1:2, repel = TRUE)
        
        
        ##individual samples with similar profile grouped together 
        fviz_pca_ind(pca,
                     col.ind = "cos2", # Color by the quality of representation
                     gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                     repel = TRUE     # Avoid text overlapping
        )
        fviz_pca_ind(pca, 
                     geom = "point", 
                     col.ind = stages, 
                     palette = "default", 
                     addEllipses = TRUE, 
                     legend.title = "Stages"
        )


#######default analysis for DESeq2 and results table
      dd1 <- DESeq(dds) ##performs nbinomWaldTest 
      res <- results(dd1)
      mcols(res, use.names = TRUE)
      summary(res)
      head(res)
      #padj column contains p-values that are adjusted for multiple testing. 
      #basemean lists the mean count values for the samples. 
      #log2FoldChange is the log2 of fold change 
      
      res$padj	= ifelse(is.na(res$padj),	1,	res$padj) #setting the NAs in padj to 1 to avoid problem later 
      
      write.table(res,	col.names=NA,	row.names=T,	file	="expressions_rosto_hydratedcyst7d_drycysts.tsv",	
                  sep ="\t")
      plotMA(res,	main="MA	plot",ylim=c(-8,8),alpha=0.01)
      
####volcano plot 
####for differential expression
    EnhancedVolcano(res,
                    lab = rownames(res),
                    x = 'log2FoldChange',
                    y = 'pvalue',
                    xlim = c(-8, 8),
                    title = 'Hydrated cysts in 7d PRD vs Dry Cysts(G. rostochiensis)',
                    subtitle = 'Differential Expression', 
                    caption = 'FC cutoff, 1.5; p-value cutoff, 10e-16',
                    pCutoff = 10e-16,
                    FCcutoff = 1.5,
                    pointSize = 3.0,
                    labSize = 3.0)
    
    topGene <- rownames(res)[which.min(res$padj)]
    topGene
    data <- plotCounts(dds , gene = topGene, intgroup = c("stages"), returnData = TRUE)
        
    ggplot(data, aes(x = stages, y = count, color = stages)) + 
          scale_y_log10() + ggtitle("Plot count for Topgene g10251.t1")+
          geom_point(position = position_jitter(width = .1, height = 0)) + 
          theme(axis.text.x= element_text(angle = 90, vjust = 0.5, size = 10))

####for J2 vs female 
          res1 <- results(dd1, 
                          contrast = c("stages", "females", "J2"))
          mcols(res1, use.names = TRUE)
          summary(res1)
          
          head(res1)
          write.table(res1,	col.names=NA,	row.names=T,	file	="expressions_rosto_femalesvsJ2.tsv",	
                      sep ="\t")
          EnhancedVolcano(res1,
                          lab = rownames(res1),
                          x = 'log2FoldChange',
                          y = 'pvalue',
                          xlim = c(-8, 8),
                          title = 'Females vs J2 (G. rostochiensis)',
                          subtitle = 'Differential Expression ', 
                          caption = 'FC cutoff, 1.5; p-value cutoff, 10e-16',
                          pCutoff = 10e-16,
                          FCcutoff = 1.5,
                          pointSize = 3.0,
                          labSize = 3.0)
          
          topGene1 <- rownames(res1)[which.min(res1$padj)]
          topGene1
          data1 <- plotCounts(dds , gene = topGene1, intgroup = c("stages"), returnData = TRUE)
          ggplot(data1, aes(x = stages, y = count, color = stages)) + 
            scale_y_log10() + ggtitle("Plot count for Topgene g1307.t1")+
            geom_point(position = position_jitter(width = .1, height = 0))+ 
            theme(axis.text.x= element_text(angle = 90, vjust = 0.5, size = 10))

####for drycysts vs females
          res2 <- results(dd1, 
                          contrast = c("stages", "dry cysts", "females"))
          head(res2)
          mcols(res2, use.names = TRUE)
          summary(res2)
          write.table(res2,	col.names=NA,	row.names=T,	file	="expressions_rosto_drycystsvsfemales.tsv",	
                      sep ="\t")
          EnhancedVolcano(res2,
                          lab = rownames(res2),
                          x = 'log2FoldChange',
                          y = 'pvalue',
                          xlim = c(-8, 9),
                          title = 'Dry Cysts vs Females (G. rostochiensis)',
                          subtitle = 'Differential Expression ', 
                          caption = 'FC cutoff, 1.5; p-value cutoff, 10e-16',
                          pCutoff = 10e-16,
                          FCcutoff = 1.5,
                          pointSize = 3.0,
                          labSize = 3.0)
          
          topGene2 <- rownames(res2)[which.min(res2$padj)]
          topGene2
          data2 <- plotCounts(dds , gene = topGene2, intgroup = c("stages"), returnData = TRUE)
          ggplot(data2, aes(x = stages, y = count, color = stages)) + 
              scale_y_log10() + ggtitle("Plot count for Topgene g1307.t1")+
              geom_point(position = position_jitter(width = .1, height = 0)) + 
              theme(axis.text.x= element_text(angle = 90, vjust = 0.5, size = 10))

####for J2 vs hydrated cysts
          res3 <- results(dd1, 
                          contrast = c("stages", "J2", "hydrated cysts"))
          head(res3)
          mcols(res3, use.names = TRUE)
          summary(res3)
          write.table(res3,	col.names=NA,	row.names=T,	file	="expressions_rosto_J2vsmales.tsv",	
                      sep ="\t")
          EnhancedVolcano(res3,
                          lab = rownames(res3),
                          x = 'log2FoldChange',
                          y = 'pvalue',
                          xlim = c(-9, 9),
                          title = 'J2 vs hydrated cysts (G. rostochiensis)',
                          subtitle = 'Differential Expression ', 
                          caption = 'FC cutoff, 1.5; p-value cutoff, 10e-16',
                          pCutoff = 10e-16,
                          FCcutoff = 1.5,
                          pointSize = 3.0,
                          labSize = 3.0)
          topGene3 <- rownames(res3)[which.min(res3$padj)]
          topGene3
          data3 <- plotCounts(dds , gene = topGene3, intgroup = c("stages"), returnData = TRUE)
          ggplot(data3, aes(x = stages, y = count, color = stages)) + 
            scale_y_log10() + ggtitle("Plot count for Topgene g14414.t2.t1")+
            geom_point(position = position_jitter(width = .1, height = 0))+ 
            theme(axis.text.x= element_text(angle = 90, vjust = 0.5, size = 10))

####for dry cysts vs J2 
          res5 <- results(dd1, 
                          contrast = c("stages", "J2", "dry cysts"))
          head(res5)
          mcols(res5, use.names = TRUE)
          summary(res5)
          write.table(res5,	col.names=NA,	row.names=T,	file	="expressions_rost0_J2_vsdrycysts.tsv",	
                      sep ="\t")
          EnhancedVolcano(res5,
                          lab = rownames(res5),
                          x = 'log2FoldChange',
                          y = 'pvalue',
                          xlim = c(-9, 9),
                          title = 'J2 vs dry cysts',
                          subtitle = 'Differential Expression', 
                          caption = 'FC cutoff, 1.5; p-value cutoff, 10e-5',
                          pCutoff = 10e-5,
                          FCcutoff = 1.5,
                          pointSize = 3.0,
                          labSize = 3.0)
          topGene5 <- rownames(res5)[which.min(res5$padj)]
          topGene5
          data5 <- plotCounts(dds , gene = topGene5, intgroup = c("stages"), returnData = TRUE)
          ggplot(data5, aes(x = stages, y = count, color = stages)) + 
            scale_y_log10() + ggtitle("Plot count for Topgene g10251.t1")+
            geom_point(position = position_jitter(width = .1, height = 0))+ 
            theme(axis.text.x= element_text(angle = 90, vjust = 0.5, size = 10))

          
          
          ##Investigating hatching transcriptome 
          res4 <- results(dd1, 
                          contrast = c("stages", "hydrated cysts(48h PRD)", "J2"))
          head(res4)
          mcols(res4, use.names = TRUE)
          summary(res4)
          write.table(res4,	col.names=NA,	row.names=T,	file	="expressions_rosto_hydratedcysts(48h)_J2.tsv",	
                      sep ="\t")
          res4$padj	= ifelse(is.na(res4$padj),	1,	res4$padj)
          EnhancedVolcano(res4,
                          lab = rownames(res4),
                          x = 'log2FoldChange',
                          y = 'pvalue',
                          xlim = c(-9, 9),
                          title = 'Hydrated cysts (48h) vs J2',
                          subtitle = 'Differential Expression', 
                          caption = 'FC cutoff, 1.5; p-value cutoff, 10e-5',
                          pCutoff = 10e-5,
                          FCcutoff = 1.5,
                          pointSize = 3.0,
                          labSize = 3.0)
          topGene4 <- rownames(res4)[which.min(res4$padj)]
          topGene4
          data4 <- plotCounts(dds , gene = topGene4, intgroup = c("stages"), returnData = TRUE)
          ggplot(data4, aes(x = stages, y = count, color = stages)) + 
            scale_y_log10() + ggtitle("Plot count for Topgene g15280.t2")+
            geom_point(position = position_jitter(width = .1, height = 0))+ 
            theme(axis.text.x= element_text(angle = 90, vjust = 0.5, size = 10))


          #########Gene Clustering 
          topVarGenes <- head(order(-rowVars(assay(rld))), 2)
          topVarGenes
          colors <- colorRampPalette( rev(brewer.pal(9, "PuOr")) )(255)
          sidecols <- c(brewer.pal(name = "Set2", n = 6), (brewer.pal(name = "Set2", n = 3)))[ rld$stages ]
          mat <- assay(rld)[ topVarGenes, ]
          mat <- mat - rowMeans(mat)
          colnames(mat) <- paste0(rld$stages,"-",rld$sizeFactor)
          heatmap(mat, trace="none", col=colors, ColSideColors=sidecols,
                  labRow=FALSE, mar=c(10,2), scale="row")



          dd2 <- DESeq(dds, test = "LRT", reduced = ~ 1 ) ##likelihood test when the full design has only one factor(stages)
          head(dd2)
          counts <- counts(dd2, normalized = TRUE)
          design <- as.data.frame(colData(dd2))
          ##LRT compares full model to the reduced model to identify significant genes. 
          ##p-values are determined by difference in the deviance between full an reduced model formula 
          ##LRT test tests whether the terms removed in the reduced model explains a significant amount of variation in the data 
          res_LRT<- results(dd2)
          head(res_LRT)
          summary(res_LRT)
          
          EnhancedVolcano(res_LRT,
                          lab = rownames(res_LRT),
                          x = 'log2FoldChange',
                          y = 'pvalue',
                          xlim = c(-9, 9),
                          title = 'Differential expression across all life stages',
                          subtitle = 'DESeq2 Analysis', 
                          caption = 'FC cutoff, 1.5; p-value cutoff, 10e-6',
                          pCutoff = 1.5,
                          FCcutoff = 10e-6,
                          pointSize = 3.0,
                          labSize = 3.0)

#########subsetting LRT results to return genes with padj < 0.0001
          sig_res_LRT <- res_LRT %>% 
            data.frame() %>% 
            rownames_to_column(var = "gene") %>%
            as_tibble() %>%
            filter(padj < 0.001)
            

              ##getting significant gene list 
              sigLRT_genes <- sig_res_LRT %>% 
                pull(gene)
              
              write.csv(sigLRT_genes, file = "significant_genes_acc_to_LRT_pallida.csv")
              length(sigLRT_genes) ### 6153 differentially expressed genes across all life stages 
              
              nrow(exp_data)
              nrow(cluster_rlog)
              
              
              ##clustering trial 
              clustering_sig_genes <- sig_res_LRT %>% 
                arrange(padj)
              length(clustering_sig_genes)
              
              ##obtaining log values for significant genes
              cluster_rlog <- rld_mat[clustering_sig_genes$gene, ]
              length(cluster_rlog)
              
              clusters <- degPatterns(cluster_rlog, metadata = design, time = "stages", col = NULL)
              
              ####THIS IS WHERE I ENCOUNTER ERRORS 
              
              
