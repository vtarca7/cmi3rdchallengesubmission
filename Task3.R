suppressPackageStartupMessages({
  library(dplyr)
  library(randomForest)
  library(readr)
  library(tidyr)
  library(purrr)
})

#submission_template = readr::read_tsv("3rdChallengeSubmissionTemplate_revised.tsv")
master_harmonized_data = read_rds("master_harmonized_data.RDS")
excludedcolumns <- c("year_of_birth","actual_day_relative_to_boost","specimen_type", "date_of_boost", "dataset","specimen_id","timepoint","race","ethnicity","visit")

master_harmonized_data[["training"]][["subject_specimen"]]$birth_year <- as.numeric(substr(master_harmonized_data[["training"]][["subject_specimen"]]$year_of_birth, 1, 4))
master_harmonized_data[["training"]][["subject_specimen"]]$dataset_year <- as.numeric(sub('_dataset', '', master_harmonized_data[["training"]][["subject_specimen"]]$dataset))
master_harmonized_data[["training"]][["subject_specimen"]]$age <- master_harmonized_data[["training"]][["subject_specimen"]]$dataset_year - master_harmonized_data[["training"]][["subject_specimen"]]$birth_year

master_harmonized_data[["challenge"]][["subject_specimen"]]$birth_year <- as.numeric(substr(master_harmonized_data[["challenge"]][["subject_specimen"]]$year_of_birth, 1, 4))
master_harmonized_data[["challenge"]][["subject_specimen"]]$dataset_year <- as.numeric(sub('_dataset', '', master_harmonized_data[["challenge"]][["subject_specimen"]]$dataset))
master_harmonized_data[["challenge"]][["subject_specimen"]]$age <- master_harmonized_data[["challenge"]][["subject_specimen"]]$dataset_year - master_harmonized_data[["challenge"]][["subject_specimen"]]$birth_year


load("selected_genes.RData")


dataframes <- list(master_harmonized_data[["training"]][["subject_specimen"]], master_harmonized_data[["training"]][["pbmc_gene_expression"]][["wide_tpm"]][selected_genes],master_harmonized_data[["training"]][["plasma_antibody_levels"]][["wide"]])
gene_expression_raw <- reduce(dataframes, full_join, by = "specimen_id") %>% drop_na()
dataframes <- list(master_harmonized_data[["challenge"]][["subject_specimen"]], master_harmonized_data[["challenge"]][["pbmc_gene_expression"]][["wide_tpm"]][selected_genes],master_harmonized_data[["challenge"]][["plasma_antibody_levels"]][["wide"]])
gene_expression_challenge <- reduce(dataframes, full_join, by = "specimen_id") %>% drop_na()





target_pred <- rlang::sym("ENSG00000277632.1")


baseline_df <- gene_expression_raw %>%
  filter(planned_day_relative_to_boost == 0)


actual_df <- gene_expression_raw %>%
  filter(planned_day_relative_to_boost == 3)


colnames(actual_df)[colnames(actual_df) == rlang::as_string(target_pred)] <- "pred"

actual_df <- actual_df %>% filter(subject_id %in% intersect(actual_df$subject_id, baseline_df$subject_id))
baseline_df <- baseline_df %>% filter(subject_id %in% intersect(actual_df$subject_id, baseline_df$subject_id))

gene_expression_training <- baseline_df %>%
  left_join(actual_df %>% select(subject_id, pred), by = "subject_id")



excludedcolumns <- c("year_of_birth","actual_day_relative_to_boost","specimen_type", "date_of_boost", "dataset","specimen_id","timepoint","race","ethnicity","visit")



gene_expression_training <- gene_expression_training[, !names(gene_expression_training) %in% excludedcolumns]
gene_expression_challenge <- gene_expression_challenge[, !names(gene_expression_challenge) %in% excludedcolumns]


gene_expression_challenge <- gene_expression_challenge %>%
  filter(planned_day_relative_to_boost == 0)


modb = randomForest(gene_expression_training[ , !(colnames(gene_expression_training) %in% c("pred", "subject_id"))], 
                    y = gene_expression_training$pred, 
                    ntree = 2000, 
                    importance = TRUE)
pred=predict(modb,gene_expression_challenge)

gene_expression_challenge["pred"] <- pred
gene_expression_challenge["rank"] <- rank(gene_expression_challenge$pred)



submission_template[,"3.1) CCL3-D3-Rank"]<-gene_expression_challenge[match(submission_template$SubjectID,gene_expression_challenge$subject_id),"rank"]



###########################################################################################


dataframes <- list(master_harmonized_data[["training"]][["subject_specimen"]], master_harmonized_data[["training"]][["pbmc_gene_expression"]][["wide_tpm"]])
gene_expression_raw <- reduce(dataframes, full_join, by = "specimen_id") %>% drop_na()
dataframes <- list(master_harmonized_data[["challenge"]][["subject_specimen"]], master_harmonized_data[["challenge"]][["pbmc_gene_expression"]][["wide_tpm"]])
gene_expression_challenge <- reduce(dataframes, full_join, by = "specimen_id") %>% drop_na()





target_pred <- rlang::sym("ENSG00000277632.1")


baseline_df <- gene_expression_raw %>%
  filter(planned_day_relative_to_boost == 0)


actual_df <- gene_expression_raw %>%
  filter(planned_day_relative_to_boost == 3)


colnames(actual_df)[colnames(actual_df) == rlang::as_string(target_pred)] <- "pred"

actual_df <- actual_df %>% filter(subject_id %in% intersect(actual_df$subject_id, baseline_df$subject_id))
baseline_df <- baseline_df %>% filter(subject_id %in% intersect(actual_df$subject_id, baseline_df$subject_id))

gene_expression_training <- baseline_df %>%
  left_join(actual_df %>% select(subject_id, pred), by = "subject_id")

gene_expression_training["pred"] <- gene_expression_training["pred"] / gene_expression_training["ENSG00000277632.1"]







gene_expression_training <- gene_expression_training[, !names(gene_expression_training) %in% excludedcolumns]
gene_expression_challenge <- gene_expression_challenge[, !names(gene_expression_challenge) %in% excludedcolumns]


gene_expression_challenge <- gene_expression_challenge %>%
  filter(planned_day_relative_to_boost == 0)


modb = randomForest(gene_expression_training[ , !(colnames(gene_expression_training) %in% c("pred", "subject_id"))], 
                    y = gene_expression_training$pred, 
                    ntree = 2000, 
                    importance = TRUE)
pred=predict(modb,gene_expression_challenge)

gene_expression_challenge["pred"] <- pred
gene_expression_challenge["rank"] <- rank(gene_expression_challenge$pred)



submission_template[,"3.2) CCL3-D3-FC-Rank"]<-gene_expression_challenge[match(submission_template$SubjectID,gene_expression_challenge$subject_id),"rank"]




write.csv(submission_template,"last_predictions.csv")

