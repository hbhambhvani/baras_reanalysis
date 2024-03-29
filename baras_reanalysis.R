

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


#unsupervised clustering analysis, starting with nuclear morph respondners
#unsupervised clustering analysis, starting with nuclear morph respondners
#unsupervised clustering analysis, starting with nuclear morph respondners

nuc_morph_resp_clustering <- hopkins_nuc_morph_merged[hopkins_nuc_morph_merged$response == 1,]

nuc_morph_resp_clustering <- t(nuc_morph_resp_clustering[,14:66])

nuc_morph_resp_clustering <- nuc_morph_resp_clustering[,colSums(is.na(nuc_morph_resp_clustering))<nrow(nuc_morph_resp_clustering)]

nuc_morph_resp_clustering <- scale(nuc_morph_resp_clustering)

# Dissimilarity matrix
d <- dist(nuc_morph_resp_clustering, method = "euclidean")

# Hierarchical clustering using Complete Linkage
hc1 <- hclust(d, method = "complete")

pdf("nuc_morph_resp_clustering.pdf", width=40, height=35)

# Do some plotting
plot(hc1)

# Close the PDF file's associated graphics device (necessary to finalize the output)
dev.off()

#moving on to nuclear morphology non responders

nuc_morph_nonresp_clustering <- hopkins_nuc_morph_merged[hopkins_nuc_morph_merged$response == 0,]

nuc_morph_nonresp_clustering <- t(nuc_morph_nonresp_clustering[,14:66])

nuc_morph_nonresp_clustering <- nuc_morph_nonresp_clustering[,colSums(is.na(nuc_morph_nonresp_clustering))<nrow(nuc_morph_nonresp_clustering)]

nuc_morph_nonresp_clustering <- scale(nuc_morph_nonresp_clustering)

d <- dist(nuc_morph_nonresp_clustering, method = "euclidean")

hc1 <- hclust(d, method = "complete")

pdf("nuc_morph_nonresp_clustering.pdf", width=40, height=35)

plot(hc1)

dev.off()

# Moving on to clustering features among responders

clusterFeat_resp_clustering <- hopkins_clustering_merged[hopkins_clustering_merged$response == 1,]

clusterFeat_resp_clustering <- t(clusterFeat_resp_clustering[,14:66])

clusterFeat_resp_clustering <- clusterFeat_resp_clustering[,colSums(is.na(clusterFeat_resp_clustering))<nrow(clusterFeat_resp_clustering)]

clusterFeat_resp_clustering <- scale(clusterFeat_resp_clustering)

d <- dist(clusterFeat_resp_clustering, method = "euclidean")

hc1 <- hclust(d, method = "complete")

plot(hc1)

pdf("clusterFeat_resp_clustering.pdf", width=40, height=35)

plot(hc1)

dev.off()

#moving on to clustering features among nonresponders

clusterFeat_nonresp_clustering <- hopkins_clustering_merged[hopkins_clustering_merged$response == 0,]

clusterFeat_nonresp_clustering <- t(clusterFeat_nonresp_clustering[,14:66])

clusterFeat_nonresp_clustering <- clusterFeat_nonresp_clustering[,colSums(is.na(clusterFeat_nonresp_clustering))<nrow(clusterFeat_nonresp_clustering)]

clusterFeat_nonresp_clustering <- scale(clusterFeat_nonresp_clustering)

d <- dist(clusterFeat_nonresp_clustering, method = "euclidean")

hc1 <- hclust(d, method = "complete")

plot(hc1)

pdf("clusterFeat_nonresp_clustering.pdf", width=40, height=35)

plot(hc1)

dev.off()


#moving on to cell population features among responders

cellpop_resp_clustering <- hopkins_cellpop_merged[hopkins_cellpop_merged$response == 1,]

cellpop_resp_clustering <- t(cellpop_resp_clustering[,15:63])

cellpop_resp_clustering <- cellpop_resp_clustering[,colSums(is.na(cellpop_resp_clustering))<nrow(cellpop_resp_clustering)]

cellpop_resp_clustering <- scale(cellpop_resp_clustering)

d <- dist(cellpop_resp_clustering, method = "euclidean")

hc1 <- hclust(d, method = "complete")

plot(hc1)

pdf("cellpop_resp_clustering.pdf", width=40, height=35)

plot(hc1)

dev.off()

#moving on to cell population features among nonresponders

cellpop_nonresp_clustering <- hopkins_cellpop_merged[hopkins_cellpop_merged$response == 0,]

cellpop_nonresp_clustering <- t(cellpop_nonresp_clustering[,15:63])

cellpop_nonresp_clustering <- cellpop_nonresp_clustering[,colSums(is.na(cellpop_nonresp_clustering))<nrow(cellpop_nonresp_clustering)]

cellpop_nonresp_clustering <- scale(cellpop_nonresp_clustering)

d <- dist(cellpop_nonresp_clustering, method = "euclidean")

hc1 <- hclust(d, method = "complete")

plot(hc1)

pdf("cellpop_nonresp_clustering.pdf", width=40, height=35)

plot(hc1)

dev.off()

### combining all of the datasets in order to develop a Cox model

test <- merge(hopkins_nuc_morph_merged, hopkins_clustering_merged, by = 'uniqueID', all.x = TRUE, sort = FALSE)

test2 <- merge(test, hopkins_cellpop_merged, by = 'uniqueID', all.x = TRUE, sort = FALSE)

test3 <- test2[is.na(test2$Area.um.2_Max) == FALSE,]



### trying models to predict response now with all features 
### trying models to predict response now with all features 
### trying models to predict response now with all features 


baras_dta <- select(test3, c('response',
                             'elongation_CoV', 
                             'Min.diameter.um_CoV',
                             'elongation_Min',
                             'Circularity_CoV',
                             'elongation_Mean',
                             'Circularity_Mean',
                             'chullArea_Min',
                             'chullLength_Min',
                             'Length.um_CoV',
                             'chullLength_Mean',
                             'chullLength_CoV',
                             'chullArea_Mean',
                             'Max.diameter.um_CoV',
                             'convexity_Max',
                             'convexity_Mean',
                             'curvMeanStat',
                             'MOI_Max',
                             'Density_of_nucleus_mean',
                             'Density_of_nucleus_min',
                             'Number_of_triangles',
                             'Perimeter_of_triangle_mean',
                             'Area_of_triangle_mean',
                             'Area_of_triangle_min',
                             'stromal_ratio',
                             'lymphocyte.dens',
                             'cancer.dens',
                             'cancer.to.lympho.avgCount',
                             'cancer.to.stroma.avgCount',
                             'ShannonH',
                             'DoC_TL_min',
                             'DoC_LT_min',
                             'DoC_TS_max'))

baras_dta <- baras_dta %>% mutate(response = as.factor(response))

#split into train and test data sets

baras_dta_split <- initial_split(baras_dta, prop = 4/5, strata = response)

baras_dta_training <- training(baras_dta_split)
baras_dta_testing  <- testing(baras_dta_split)

#create cross validation folds

baras_dta_folds <- vfold_cv(baras_dta_training, v = 5)

### random forest modeling
### random forest modeling
### random forest modeling

rf_rec <- 
  recipe(response ~ ., data = baras_dta_training) %>% 
  step_zv(all_predictors())

prep(rf_rec)

#random forest model specs

rf_mod <- 
  rand_forest(trees = 1000,
              mtry = tune(),
              min_n = tune()) %>% 
  set_engine("ranger") %>% 
  set_mode("classification")

#create a workflow to pair recipe and the model

rf_workflow <- workflow() %>%
  add_model(rf_mod) %>%
  add_recipe(rf_rec)

#tune hyperparameters across folds 

rf_tune <- 
  rf_workflow %>% 
  tune_grid(resamples = baras_dta_folds,
            grid = 11)

#get tuning results 

rf_tune_results <- collect_metrics(rf_tune)

#finalize the workflow

final_rf <- rf_workflow %>%
  finalize_workflow(select_best(rf_tune, metric = 'roc_auc'))

#develop final model and assess on hold-out test set

rf_final_fit <- last_fit(final_rf, baras_dta_split)

collect_metrics(rf_final_fit)

#trying out SVM


baras_dta_rec_svm <- 
  recipe(response ~ ., data = baras_dta_training) %>% 
  step_zv(all_predictors())

prep(baras_dta_rec_svm)

svm_mod <-
  svm_rbf(cost = tune(), rbf_sigma = tune()) %>%
  set_mode("classification") %>%
  set_engine("kernlab")

svm_workflow <- 
  workflow() %>% 
  add_model(svm_mod) %>% 
  add_recipe(baras_dta_rec_svm)

svm_tune <- 
  svm_workflow %>% 
  tune_grid(resamples = baras_dta_folds,
            grid = 10)

test <- collect_metrics(svm_tune)

final_svm <- svm_workflow %>%
  finalize_workflow(select_best(svm_tune, metric = 'roc_auc'))

#develop final model

svm_final_fit <- last_fit(final_svm, baras_dta_split)

collect_metrics(svm_final_fit)

test <- collect_predictions(nn_final_fit)
















