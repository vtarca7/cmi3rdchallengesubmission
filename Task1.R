suppressPackageStartupMessages({
  library(dplyr)
  library(randomForest)
  library(readr)
  library(tidyr)
  library(purrr)
})

submission_template = readr::read_tsv("3rdChallengeSubmissionTemplate_revised.tsv")
master_harmonized_data = read_rds("master_harmonized_data.RDS")
excludedcolumns <- c("year_of_birth", "date_of_boost", "dataset","specimen_id","subject_id","subject_id.1","timepoint","race","ethnicity","visit")



master_harmonized_data[["training"]][["subject_specimen"]]$birth_year <- as.numeric(substr(master_harmonized_data[["training"]][["subject_specimen"]]$year_of_birth, 1, 4))
master_harmonized_data[["training"]][["subject_specimen"]]$dataset_year <- as.numeric(sub('_dataset', '', master_harmonized_data[["training"]][["subject_specimen"]]$dataset))
master_harmonized_data[["training"]][["subject_specimen"]]$age <- master_harmonized_data[["training"]][["subject_specimen"]]$dataset_year - master_harmonized_data[["training"]][["subject_specimen"]]$birth_year

master_harmonized_data[["challenge"]][["subject_specimen"]]$birth_year <- as.numeric(substr(master_harmonized_data[["challenge"]][["subject_specimen"]]$year_of_birth, 1, 4))
master_harmonized_data[["challenge"]][["subject_specimen"]]$dataset_year <- as.numeric(sub('_dataset', '', master_harmonized_data[["challenge"]][["subject_specimen"]]$dataset))
master_harmonized_data[["challenge"]][["subject_specimen"]]$age <- master_harmonized_data[["challenge"]][["subject_specimen"]]$dataset_year - master_harmonized_data[["challenge"]][["subject_specimen"]]$birth_year




dataframes <- list(master_harmonized_data[["training"]][["subject_specimen"]], master_harmonized_data[["training"]][["plasma_antibody_levels"]][["wide"]])
plasma_antibody_levels_raw <- reduce(dataframes, full_join, by = "specimen_id") %>% drop_na()
dataframes <- list(master_harmonized_data[["challenge"]][["subject_specimen"]], master_harmonized_data[["challenge"]][["plasma_antibody_levels"]][["wide"]])
plasma_antibody_challenge <- reduce(dataframes, full_join, by = "specimen_id") %>% drop_na()



target_pred <- rlang::sym("IgG_PT")


baseline_df <- plasma_antibody_levels_raw %>%
  filter(planned_day_relative_to_boost == 0)


actual_df <- plasma_antibody_levels_raw %>%
  filter(planned_day_relative_to_boost == 14)


colnames(actual_df)[colnames(actual_df) == rlang::as_string(target_pred)] <- "pred"

actual_df <- actual_df %>% filter(subject_id %in% intersect(actual_df$subject_id, baseline_df$subject_id))
baseline_df <- baseline_df %>% filter(subject_id %in% intersect(actual_df$subject_id, baseline_df$subject_id))

plasma_antibody_training <- baseline_df %>%
  left_join(actual_df %>% select(subject_id, pred), by = "subject_id")


plasma_antibody_training <- plasma_antibody_training[, !names(plasma_antibody_training) %in% excludedcolumns]
plasma_antibody_challenge <- plasma_antibody_challenge[, !names(plasma_antibody_challenge) %in% excludedcolumns]


plasma_antibody_challenge <- plasma_antibody_challenge %>%
  filter(planned_day_relative_to_boost == 0)


modb = randomForest(plasma_antibody_training[ , !(colnames(plasma_antibody_training) %in% c("pred", "subject_id"))], 
                    y = plasma_antibody_training$pred, 
                    ntree = 2000, 
                    importance = TRUE)
pred=predict(modb,plasma_antibody_challenge)

plasma_antibody_challenge["pred"] <- pred
plasma_antibody_challenge["rank"] <- rank(plasma_antibody_challenge$pred)



submission_template[,"1.1) IgG-PT-D14-titer-Rank"]<-plasma_antibody_challenge[match(submission_template$SubjectID,plasma_antibody_challenge$subject_id),"rank"]



#################################################################1.2############################################################################

dataframes <- list(master_harmonized_data[["training"]][["subject_specimen"]], master_harmonized_data[["training"]][["plasma_antibody_levels"]][["wide"]])
plasma_antibody_levels_raw <- reduce(dataframes, full_join, by = "specimen_id") %>% drop_na()
dataframes <- list(master_harmonized_data[["challenge"]][["subject_specimen"]], master_harmonized_data[["challenge"]][["plasma_antibody_levels"]][["wide"]])
plasma_antibody_challenge <- reduce(dataframes, full_join, by = "specimen_id") %>% drop_na()

target_pred <- rlang::sym("IgG_PT")


baseline_df <- plasma_antibody_levels_raw %>%
  filter(planned_day_relative_to_boost == 0)


actual_df <- plasma_antibody_levels_raw %>%
  filter(planned_day_relative_to_boost == 14)


colnames(actual_df)[colnames(actual_df) == rlang::as_string(target_pred)] <- "pred"

actual_df <- actual_df %>% filter(subject_id %in% intersect(actual_df$subject_id, baseline_df$subject_id))
baseline_df <- baseline_df %>% filter(subject_id %in% intersect(actual_df$subject_id, baseline_df$subject_id))

plasma_antibody_training <- baseline_df %>%
  left_join(actual_df %>% select(subject_id, pred), by = "subject_id")

plasma_antibody_training["pred"] <- plasma_antibody_training["pred"] / plasma_antibody_training["IgG_PT"]



plasma_antibody_training <- plasma_antibody_training[, !names(plasma_antibody_training) %in% excludedcolumns]
plasma_antibody_challenge <- plasma_antibody_challenge[, !names(plasma_antibody_challenge) %in% excludedcolumns]


plasma_antibody_challenge <- plasma_antibody_challenge %>%
  filter(planned_day_relative_to_boost == 0)


modb = randomForest(plasma_antibody_training[ , !(colnames(plasma_antibody_training) %in% c("pred", "subject_id"))], 
                    y = plasma_antibody_training$pred, 
                    ntree = 2000, 
                    importance = TRUE)
pred=predict(modb,plasma_antibody_challenge)

plasma_antibody_challenge["pred"] <- pred
plasma_antibody_challenge["rank"] <- rank(plasma_antibody_challenge$pred)



submission_template[,"1.2) IgG-PT-D14-FC-Rank"]<-plasma_antibody_challenge[match(submission_template$SubjectID,plasma_antibody_challenge$subject_id),"rank"]

