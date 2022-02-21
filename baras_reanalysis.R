

#Attempting to replicate results from Mi et al.,  2021 (PMID 34622225)

library(readxl)
library(survival)
library(survminer)
library(tableone)
library(tidyr)
library(dplyr)
library(ggplot2)
library(Hmisc)
library(corrplot)


hopkins_clinical <- read_excel('/Users/hridaybhambhvani/Downloads/Clinicopathologic Features-1.xlsx')
hopkins_clinical$patientID <- hopkins_clinical$Case
hopkins_clinical$age <- hopkins_clinical$`AGE AT OPERATION`
hopkins_clinical$clin_t_stage <- hopkins_clinical$`CLINICAL T-STAGE`
hopkins_clinical$path_t_stage <- hopkins_clinical$`PATHOLOGY T-STAGE`
hopkins_clinical$path_n_stage <- hopkins_clinical$`PATHOLOGY N-STAGE`
hopkins_clinical$response <- hopkins_clinical$`(Nonresponder ypT2,3,4) = 1`
hopkins_clinical <- mutate(hopkins_clinical, response = (response * -1) + 1) 

hopkins_clinical <- hopkins_clinical[, c('Case', 'patientID', 'TMA', 'age', 'RACE', 'SEX', 'clin_t_stage', 'path_t_stage', 'path_n_stage', 'CHEMO_CYCLES', 'response')]

hopkins_cell_pop <- read_excel('/Users/hridaybhambhvani/Downloads/CellPopulation_features963.xlsx')
hopkins_cell_pop1042 <- hopkins_cell_pop[hopkins_cell_pop$TMA.set == 1042,]
hopkins_cell_pop963 <- hopkins_cell_pop[hopkins_cell_pop$TMA.set == 963,]

hopkins_clustering1042 <- read_excel('/Users/hridaybhambhvani/Downloads/Clustering_features1042.xlsx')

hopkins_clustering963 <- read_excel('/Users/hridaybhambhvani/Downloads/Clustering_features963.xlsx')

hopkins_nuc_morph1042 <- read_excel('/Users/hridaybhambhvani/Downloads/NucleusMorphology_features1042.xlsx')

hopkins_nuc_morph963 <- read_excel('/Users/hridaybhambhvani/Downloads/NucleusMorphology_features963.xlsx')

hopkins_image_texture_1042 <- read_excel('/Users/hridaybhambhvani/Downloads/ImageTexture_features1042.xlsx')

hopkins_image_texture_963 <- read_excel('/Users/hridaybhambhvani/Downloads/ImageTexture_features963.xlsx')





#note: can be multiple TMA's for a singular patient
#merging with nuclear morphology features

test1042 <- hopkins_clinical[hopkins_clinical$TMA == 1042,]

test1042 <- merge(test1042, hopkins_image_texture_1042[,c('patientID', 'TMA.core')], by = 'patientID', all.x = TRUE, sort = FALSE)

test1042$uniqueID <- paste(test1042$patientID, test1042$TMA, test1042$TMA.core, sep =  '-')

clin_plus_nuc_morph1042 <- merge(test1042, hopkins_nuc_morph1042, by = 'TMA.core', all.x = TRUE, sort = FALSE)


test963 <- hopkins_clinical[hopkins_clinical$TMA == 963,]

test963 <- merge(test963, hopkins_image_texture_963[,c('patientID', 'TMA.core')], by = 'patientID', all.x = TRUE, sort = FALSE)

test963$uniqueID <- paste(test963$patientID, test963$TMA, test963$TMA.core, sep =  '-')

clin_plus_nuc_morph963 <- merge(test963, hopkins_nuc_morph963, by = 'TMA.core', all.x = TRUE, sort = FALSE)



hopkins_nuc_morph_merged <- rbind(clin_plus_nuc_morph1042, clin_plus_nuc_morph963)

#univarite logistic regression analysis for nuclear morphology features

univ_formulas <- sapply(names(hopkins_nuc_morph_merged)[14:73],function(x)as.formula(paste('response~',x)))

univ_models <- lapply(univ_formulas, function(x){glm(x,family = 'binomial', data=hopkins_nuc_morph_merged)})

univ_results <- lapply(univ_models,function(x){return(cbind(exp(coef(x)), exp(confint(x)),summary(x)$coefficients[,4]))})

#merging with clustering features

clin_plus_clustering_1042 <- merge(test1042, hopkins_clustering1042, by = 'TMA.core', all.x = TRUE, sort = FALSE)

clin_plus_clustering_963 <- merge(test963, hopkins_clustering963, by = 'TMA.core', all.x = TRUE, sort = FALSE)

hopkins_clustering_merged <- rbind(clin_plus_clustering_1042, clin_plus_clustering_963)
hopkins_clustering_merged$Density_of_nucleus_CoV <- as.numeric(hopkins_clustering_merged$Density_of_nucleus_CoV)

#univariate logistic regrssion analysis for clustering features 

univ_formulas <- sapply(names(hopkins_clustering_merged)[14:66],function(x)as.formula(paste('response~',x)))

univ_models <- lapply(univ_formulas, function(x){glm(x,family = 'binomial', data=hopkins_clustering_merged)})

univ_results <- lapply(univ_models,function(x){return(cbind(exp(coef(x)), exp(confint(x)),summary(x)$coefficients[,4]))})


#merging with cell population features

clin_plus_cellpop_1042 <- merge(test1042, hopkins_cell_pop1042, by = 'TMA.core', all.x = TRUE, sort = FALSE)

clin_plus_cellpop_963 <- merge(test963, hopkins_cell_pop963, by = 'TMA.core', all.x = TRUE, sort = FALSE)

hopkins_cellpop_merged <- rbind(clin_plus_cellpop_1042, clin_plus_cellpop_963)

#univariate logistic regrssion analysis for cell population features 

univ_formulas <- sapply(names(hopkins_cellpop_merged)[15:63],function(x)as.formula(paste('response~',x)))

univ_models <- lapply(univ_formulas, function(x){glm(x,family = 'binomial', data=hopkins_cellpop_merged)})

univ_results <- lapply(univ_models,function(x){return(cbind(exp(coef(x)), exp(confint(x)),summary(x)$coefficients[,4]))})

#merging with image texture features

clin_plus_imgfeat_1042 <- merge(test1042, hopkins_image_texture_1042, by = 'TMA.core', all.x = TRUE, sort = FALSE)

clin_plus_imgfeat_963 <- merge(test963, hopkins_image_texture_963, by = 'TMA.core', all.x = TRUE, sort = FALSE)

hopkins_imgfeat_merged <- rbind(clin_plus_imgfeat_1042, clin_plus_imgfeat_963)

#univariate logistic regrssion analysis for image texture features 

univ_formulas <- sapply(names(hopkins_imgfeat_merged)[16:500],function(x)as.formula(paste('response~',x)))

univ_models <- lapply(univ_formulas, function(x){glm(x,family = 'binomial', data=hopkins_imgfeat_merged)})

univ_results <- lapply(univ_models,function(x){return(cbind(exp(coef(x)), exp(confint(x)),summary(x)$coefficients[,4]))})









hopkins_clustering_merged


myVars <- names(hopkins_clustering_merged)
catVars <- c("response")
tab3 <- CreateTableOne(vars = myVars, data = hopkins_clustering_merged, factorVars = catVars)
print(tab3, showAllLevels = TRUE, quote = TRUE, noSpaces = TRUE)




#making correlation matrices for nuclear morphology features
#making correlation matrices for nuclear morphology features
#making correlation matrices for nuclear morphology features


nuc_cor_dta <- hopkins_nuc_morph_merged

nuc_cor_dta <- nuc_cor_dta[,14:73]

nuc_cor_dta <- subset(nuc_cor_dta, select = -c(Solidity_Max, n_indentation_Min, n_protrusion_Min))

nuc_cor <-  cor(nuc_cor_dta, use = 'complete.obs')

p1 <- { # Prepare the Corrplot 
  corrplot(nuc_cor, method = 'color', type = "upper", order = "AOE", 
           tl.col = "black", tl.cex = 0.3);
  # Call the recordPlot() function to record the plot
  recordPlot()
}

ggsave(filename = "p1.pdf", plot = replayPlot(p1))

# ++++++++++++++++++++++++++++
# flattenCorrMatrix
# ++++++++++++++++++++++++++++
# cormat : matrix of the correlation coefficients
# pmat : matrix of the correlation p-values
flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut],
    p = pmat[ut]
  )
}

res2 <- rcorr(as.matrix(nuc_cor_dta))

nuc_cor_p <- flattenCorrMatrix(res2$r, res2$P)

write.csv(nuc_cor_p, "nuc_morph_cor.csv", row.names = FALSE)

#remaking nuclear feature correlation matrix for Spearman's coefficient

nuc_cor_spearman <-  cor(nuc_cor_dta, use = 'complete.obs', method = 'spearman')

p1 <- { # Prepare the Corrplot 
  corrplot(nuc_cor_spearman, method = 'color', type = "upper", order = "AOE", 
           tl.col = "black", tl.cex = 0.3);
  # Call the recordPlot() function to record the plot
  recordPlot()
}

ggsave(filename = "p1.pdf", plot = replayPlot(p1))

res2 <- rcorr(as.matrix(nuc_cor_dta), type = 'spearman')

nuc_cor_p_sp <- flattenCorrMatrix(res2$r, res2$P)

write.csv(nuc_cor_p_sp, "nuc_morph_cor_spearman.csv", row.names = FALSE)


#making correlation matrices for clustering  features
#making correlation matrices for clustering  features
#making correlation matrices for clustering  features


hopkins_clustering_dta <- hopkins_clustering_merged

clustering_cor_dta <- hopkins_clustering_dta[,14:66]

cluster_cor <-  cor(clustering_cor_dta, use = 'complete.obs')

p1 <- { # Prepare the Corrplot 
  corrplot(cluster_cor, method = 'color', type = "upper", order = "AOE", 
           tl.col = "black", tl.cex = 0.3);
  # Call the recordPlot() function to record the plot
  recordPlot()
}

ggsave(filename = "cluster_cor.pdf", plot = replayPlot(p1))

cluster_cor_with_p <- rcorr(as.matrix(clustering_cor_dta))

cluster_cor_p <- flattenCorrMatrix(cluster_cor_with_p$r, cluster_cor_with_p$P)

write.csv(cluster_cor_p, "cluster_cor.csv", row.names = FALSE)


#remaking cluster feature correlation matrix for Spearman's coefficient

cluster_cor_spearman <-  cor(clustering_cor_dta, use = 'complete.obs', method = 'spearman')

p1 <- { # Prepare the Corrplot 
  corrplot(cluster_cor_spearman, method = 'color', type = "upper", order = "AOE", 
           tl.col = "black", tl.cex = 0.3);
  # Call the recordPlot() function to record the plot
  recordPlot()
}

ggsave(filename = "p1.pdf", plot = replayPlot(p1))

cluster_cor_with_p_sp <- rcorr(as.matrix(clustering_cor_dta), type = 'spearman')

cluster_cor_p_sp <- flattenCorrMatrix(res2$r, res2$P)

write.csv(cluster_cor_p_sp, "cluster_cor_spearman.csv", row.names = FALSE)


#making correlation matrices for cell population features
#making correlation matrices for cell population features
#making correlation matrices for cell population features

hopkins_cellpop_merged

cellpop_cor_dta <- hopkins_cellpop_merged

cellpop_cor_dta <- cellpop_cor_dta[,15:63]

cellpop_cor <-  cor(cellpop_cor_dta, use = 'complete.obs')

p1 <- { # Prepare the Corrplot 
  corrplot(cellpop_cor, method = 'color', type = "upper", order = "AOE", 
           tl.col = "black", tl.cex = 0.3);
  # Call the recordPlot() function to record the plot
  recordPlot()
}

ggsave(filename = "cellpop_cor.pdf", plot = replayPlot(p1))

cellpop_cor_p <- rcorr(as.matrix(cellpop_cor_dta))

cellpop_cor_with_p <- flattenCorrMatrix(cellpop_cor_p$r, cellpop_cor_p$P)

write.csv(cellpop_cor_with_p, "cellpop_cor.csv", row.names = FALSE)

#remaking cell population feature correlation matrix for Spearman's coefficient

cellpop_cor_spearman <-  cor(cellpop_cor_dta, use = 'complete.obs', method = 'spearman')

p1 <- { # Prepare the Corrplot 
  corrplot(cellpop_cor_spearman, method = 'color', type = "upper", order = "AOE", 
           tl.col = "black", tl.cex = 0.3);
  # Call the recordPlot() function to record the plot
  recordPlot()
}

ggsave(filename = "p1.pdf", plot = replayPlot(p1))

cellpop_cor_p_sp <- rcorr(as.matrix(cellpop_cor_dta), type = 'spearman')

cellpop_cor_with_p_sp <- flattenCorrMatrix(cellpop_cor_p_sp$r, cellpop_cor_p_sp$P)

write.csv(cellpop_cor_with_p_sp, "cellpop_cor_spearman.csv", row.names = FALSE)


#repeating logistic regression analysis for nuclear morphology features with z-scored values
#repeating logistic regression analysis for nuclear morphology features with z-scored values
#repeating logistic regression analysis for nuclear morphology features with z-scored values

hopkins_nuc_morph_merged_zscore <- 
  hopkins_nuc_morph_merged %>% 
  mutate_at(vars(Area.um.2_Max:MOI_CoV), scale)

hopkins_nuc_morph_merged_zscore <- subset(hopkins_nuc_morph_merged_zscore,
                                          select = -c(Solidity_Max, n_indentation_Min, n_protrusion_Min))



univ_formulas <- sapply(names(hopkins_nuc_morph_merged_zscore)[14:70],function(x)as.formula(paste('response~',x)))

univ_models <- lapply(univ_formulas, function(x){glm(x,family = 'binomial', data=hopkins_nuc_morph_merged_zscore)})

univ_results <- lapply(univ_models,function(x){return(cbind(exp(coef(x)), exp(confint(x)),summary(x)$coefficients[,4]))})


#repeating logistic regression analysis for clustering morphology features with z-scored values
#repeating logistic regression analysis for clustering morphology features with z-scored values
#repeating logistic regression analysis for clustering morphology features with z-scored values

hopkins_clustering_merged_zscore <- 
  hopkins_clustering_merged %>% 
  mutate_at(vars(Density_of_nucleus_mean:COrE_SumVariance_min), scale)

univ_formulas <- sapply(names(hopkins_clustering_merged_zscore)[14:66],function(x)as.formula(paste('response~',x)))

univ_models <- lapply(univ_formulas, function(x){glm(x,family = 'binomial', data=hopkins_clustering_merged_zscore)})

univ_results <- lapply(univ_models,function(x){return(cbind(exp(coef(x)), exp(confint(x)),summary(x)$coefficients[,4]))})




#repeating logistic regression analysis for cell population features with z-scored values
#repeating logistic regression analysis for cell population features with z-scored values
#repeating logistic regression analysis for cell population features with z-scored values

hopkins_cellpop_merged_zscore <- 
  hopkins_cellpop_merged %>% 
  mutate_at(vars(Kcross.Cancer2lympho:localized_LS), scale)

univ_formulas <- sapply(names(hopkins_cellpop_merged_zscore)[15:63],function(x)as.formula(paste('response~',x)))

univ_models <- lapply(univ_formulas, function(x){glm(x,family = 'binomial', data=hopkins_cellpop_merged_zscore)})

univ_results <- lapply(univ_models,function(x){return(cbind(exp(coef(x)), exp(confint(x)),summary(x)$coefficients[,4]))})



x <- iris[,c(1, 3)]

stats::kmeans(x, centers = 3, nstart = 10)










