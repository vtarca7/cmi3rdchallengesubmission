suppressPackageStartupMessages({
  library(dplyr)
  library(randomForest)
  library(readr)
  library(tidyr)
  library(purrr)
})

#submission_template = readr::read_tsv("3rdChallengeSubmissionTemplate_revised.tsv")
master_harmonized_data = read_rds("master_harmonized_data.RDS")
master_harmonized_datai <- master_harmonized_data
excludedcolumns <- c("year_of_birth","actual_day_relative_to_boost","specimen_type", "date_of_boost", "dataset","specimen_id","timepoint","race","ethnicity","visit")

master_harmonized_data[["training"]][["subject_specimen"]]$birth_year <- as.numeric(substr(master_harmonized_data[["training"]][["subject_specimen"]]$year_of_birth, 1, 4))
master_harmonized_data[["training"]][["subject_specimen"]]$dataset_year <- as.numeric(sub('_dataset', '', master_harmonized_data[["training"]][["subject_specimen"]]$dataset))
master_harmonized_data[["training"]][["subject_specimen"]]$age <- master_harmonized_data[["training"]][["subject_specimen"]]$dataset_year - master_harmonized_data[["training"]][["subject_specimen"]]$birth_year

master_harmonized_data[["challenge"]][["subject_specimen"]]$birth_year <- as.numeric(substr(master_harmonized_data[["challenge"]][["subject_specimen"]]$year_of_birth, 1, 4))
master_harmonized_data[["challenge"]][["subject_specimen"]]$dataset_year <- as.numeric(sub('_dataset', '', master_harmonized_data[["challenge"]][["subject_specimen"]]$dataset))
master_harmonized_data[["challenge"]][["subject_specimen"]]$age <- master_harmonized_data[["challenge"]][["subject_specimen"]]$dataset_year - master_harmonized_data[["challenge"]][["subject_specimen"]]$birth_year


dataframes <- list(master_harmonized_data[["training"]][["subject_specimen"]], master_harmonized_data[["training"]][["pbmc_cell_frequency"]][["wide"]])
monocytes_freq_raw <- reduce(dataframes, full_join, by = "specimen_id") %>% drop_na()
dataframes <- list(master_harmonized_data[["challenge"]][["subject_specimen"]], master_harmonized_data[["challenge"]][["pbmc_cell_frequency"]][["wide"]])
monocytes_freq_challenge <- reduce(dataframes, full_join, by = "specimen_id") %>% drop_na()



target_pred <- rlang::sym("Monocytes")


baseline_df <- monocytes_freq_raw %>%
  filter(planned_day_relative_to_boost == 0)


actual_df <- monocytes_freq_raw %>%
  filter(planned_day_relative_to_boost == 1)


colnames(actual_df)[colnames(actual_df) == rlang::as_string(target_pred)] <- "pred"

actual_df <- actual_df %>% filter(subject_id %in% intersect(actual_df$subject_id, baseline_df$subject_id))
baseline_df <- baseline_df %>% filter(subject_id %in% intersect(actual_df$subject_id, baseline_df$subject_id))

monocytes_freq_training <- baseline_df %>%
  left_join(actual_df %>% select(subject_id, pred), by = "subject_id")







monocytes_freq_training <- monocytes_freq_training[, !names(monocytes_freq_training) %in% excludedcolumns]
monocytes_freq_challenge <- monocytes_freq_challenge[, !names(monocytes_freq_challenge) %in% excludedcolumns]


monocytes_freq_challenge <- monocytes_freq_challenge %>%
  filter(planned_day_relative_to_boost == 0)


modb = randomForest(monocytes_freq_training[ , !(colnames(monocytes_freq_training) %in% c("pred", "subject_id"))], 
                    y = monocytes_freq_training$pred, 
                    ntree = 2000)
pred=predict(modb,monocytes_freq_challenge)

monocytes_freq_challenge["pred"] <- pred
monocytes_freq_challenge["rank"] <- rank(monocytes_freq_challenge$pred)



submission_template[,"2.1) Monocytes-D1-Rank"]<-monocytes_freq_challenge[match(submission_template$SubjectID,monocytes_freq_challenge$subject_id),"rank"]



#######################################################################2.2#####################################################################



dataframes <- list(master_harmonized_data[["training"]][["subject_specimen"]], master_harmonized_data[["training"]][["pbmc_cell_frequency"]][["wide"]])
monocytes_freq_raw <- reduce(dataframes, full_join, by = "specimen_id") %>% drop_na()
dataframes <- list(master_harmonized_data[["challenge"]][["subject_specimen"]], master_harmonized_data[["challenge"]][["pbmc_cell_frequency"]][["wide"]])
monocytes_freq_challenge <- reduce(dataframes, full_join, by = "specimen_id") %>% drop_na()



target_pred <- rlang::sym("Monocytes")


baseline_df <- monocytes_freq_raw %>%
  filter(planned_day_relative_to_boost == 0)


actual_df <- monocytes_freq_raw %>%
  filter(planned_day_relative_to_boost == 1)


colnames(actual_df)[colnames(actual_df) == rlang::as_string(target_pred)] <- "pred"

actual_df <- actual_df %>% filter(subject_id %in% intersect(actual_df$subject_id, baseline_df$subject_id))
baseline_df <- baseline_df %>% filter(subject_id %in% intersect(actual_df$subject_id, baseline_df$subject_id))

monocytes_freq_training <- baseline_df %>%
  left_join(actual_df %>% select(subject_id, pred), by = "subject_id")

monocytes_freq_training["pred"] <- monocytes_freq_training["pred"] / monocytes_freq_training["Monocytes"]


monocytes_freq_training <- monocytes_freq_training[, !names(monocytes_freq_training) %in% excludedcolumns]
monocytes_freq_challenge <- monocytes_freq_challenge[, !names(monocytes_freq_challenge) %in% excludedcolumns]


monocytes_freq_challenge <- monocytes_freq_challenge %>%
  filter(planned_day_relative_to_boost == 0)


modb = randomForest(monocytes_freq_training[ , !(colnames(monocytes_freq_training) %in% c("pred", "subject_id"))], 
                    y = monocytes_freq_training$pred, 
                    ntree = 2000, 
                    importance = TRUE)
pred=predict(modb,monocytes_freq_challenge)

monocytes_freq_challenge["pred"] <- pred
monocytes_freq_challenge["rank"] <- rank(monocytes_freq_challenge$pred)



submission_template[,"2.2) Monocytes-D1-FC-Rank"]<-monocytes_freq_challenge[match(submission_template$SubjectID,monocytes_freq_challenge$subject_id),"rank"]
