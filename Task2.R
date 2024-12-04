#Figure out how to impute
suppressPackageStartupMessages({
  library(verification)
  library(pROC)
  library(discretefit)
  library(dplyr)
  library(ggplot2)
  library(tidyr)
  library(tidyverse)
  library(corrplot)
  library(data.table)
  library(readr)
  library(kableExtra)
  library(glmnet)
  library(here)
  library(kableExtra)
  library(formattable)
  library(colorRamp2)
  #library(ComplexHeatmap)
  library(GetoptLong)
  library(ellipse)
  library("ggplot2")                     
  library("GGally")
  library(randomForest)
})
#Reading in raw data
#submission_template = readr::read_tsv("3rdChallengeSubmissionTemplate_revised.tsv")
master_harmonized_data = read_rds("master_harmonized_data.RDS")
master_harmonized_datai <- master_harmonized_data

master_harmonized_data[["training"]][["subject_specimen"]]$birth_year <- as.numeric(substr(master_harmonized_data[["training"]][["subject_specimen"]]$year_of_birth, 1, 4))
master_harmonized_data[["training"]][["subject_specimen"]]$dataset_year <- as.numeric(sub('_dataset', '', master_harmonized_data[["training"]][["subject_specimen"]]$dataset))
master_harmonized_data[["training"]][["subject_specimen"]]$age <- master_harmonized_data[["training"]][["subject_specimen"]]$dataset_year - master_harmonized_data[["training"]][["subject_specimen"]]$birth_year

master_harmonized_data[["challenge"]][["subject_specimen"]]$birth_year <- as.numeric(substr(master_harmonized_data[["challenge"]][["subject_specimen"]]$year_of_birth, 1, 4))
master_harmonized_data[["challenge"]][["subject_specimen"]]$dataset_year <- as.numeric(sub('_dataset', '', master_harmonized_data[["challenge"]][["subject_specimen"]]$dataset))
master_harmonized_data[["challenge"]][["subject_specimen"]]$age <- master_harmonized_data[["challenge"]][["subject_specimen"]]$dataset_year - master_harmonized_data[["challenge"]][["subject_specimen"]]$birth_year


#Merge SubjectSpecimin table with plasma levels df
dataframes <- list(master_harmonized_data[["training"]][["subject_specimen"]], master_harmonized_data[["training"]][["pbmc_cell_frequency"]][["wide"]])
monocytes_freq_raw <- reduce(dataframes, full_join, by = "specimen_id") %>% drop_na()
dataframes <- list(master_harmonized_data[["challenge"]][["subject_specimen"]], master_harmonized_data[["challenge"]][["pbmc_cell_frequency"]][["wide"]])
monocytes_freq_challenge <- reduce(dataframes, full_join, by = "specimen_id") %>% drop_na()



target_pred <- rlang::sym("Monocytes")

# Filter for baseline day
baseline_df <- monocytes_freq_raw %>%
  filter(planned_day_relative_to_boost == 0)

# Filter for actual day and select the target_pred column and subject_id
actual_df <- monocytes_freq_raw %>%
  filter(planned_day_relative_to_boost == 1)

# Rename target_pred column to "pred"
colnames(actual_df)[colnames(actual_df) == rlang::as_string(target_pred)] <- "pred"

actual_df <- actual_df %>% filter(subject_id %in% intersect(actual_df$subject_id, baseline_df$subject_id))
baseline_df <- baseline_df %>% filter(subject_id %in% intersect(actual_df$subject_id, baseline_df$subject_id))

monocytes_freq_training <- baseline_df %>%
  left_join(actual_df %>% select(subject_id, pred), by = "subject_id")



#Do we filter out different types of Monocyctes?

excludedcolumns <- c("year_of_birth","actual_day_relative_to_boost","specimen_type", "date_of_boost", "dataset","specimen_id","timepoint","race","ethnicity","visit")


#Remove Columns that won't be 7used
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




#######test
set.seed(1);
monocytes_freq_training=data.frame(monocytes_freq_training)
tr=sample(1:dim(monocytes_freq_training)[1],dim(monocytes_freq_training)[1]/2)
te=setdiff(1:dim(monocytes_freq_training)[1],tr)

modb = randomForest(monocytes_freq_training[tr , !(colnames(monocytes_freq_training) %in% c("pred", "subject_id"))], 
                    y = monocytes_freq_training$pred[tr], 
                    ntree = 2000)
pred=predict(modb,monocytes_freq_training[te , !(colnames(monocytes_freq_training) %in% c("pred", "subject_id"))])
plot(monocytes_freq_training[te,"pred"],pred)
cor(monocytes_freq_training[te,"pred"],pred)



#######################################################################2.2#####################################################################


#Merge SubjectSpecimin table with plasma levels df
dataframes <- list(master_harmonized_data[["training"]][["subject_specimen"]], master_harmonized_data[["training"]][["pbmc_cell_frequency"]][["wide"]])
monocytes_freq_raw <- reduce(dataframes, full_join, by = "specimen_id") %>% drop_na()
dataframes <- list(master_harmonized_data[["challenge"]][["subject_specimen"]], master_harmonized_data[["challenge"]][["pbmc_cell_frequency"]][["wide"]])
monocytes_freq_challenge <- reduce(dataframes, full_join, by = "specimen_id") %>% drop_na()



target_pred <- rlang::sym("Monocytes")

# Filter for baseline day
baseline_df <- monocytes_freq_raw %>%
  filter(planned_day_relative_to_boost == 0)

# Filter for actual day and select the target_pred column and subject_id
actual_df <- monocytes_freq_raw %>%
  filter(planned_day_relative_to_boost == 1)

# Rename target_pred column to "pred"
colnames(actual_df)[colnames(actual_df) == rlang::as_string(target_pred)] <- "pred"

actual_df <- actual_df %>% filter(subject_id %in% intersect(actual_df$subject_id, baseline_df$subject_id))
baseline_df <- baseline_df %>% filter(subject_id %in% intersect(actual_df$subject_id, baseline_df$subject_id))

monocytes_freq_training <- baseline_df %>%
  left_join(actual_df %>% select(subject_id, pred), by = "subject_id")

monocytes_freq_training["pred"] <- monocytes_freq_training["pred"] / monocytes_freq_training["Monocytes"]

#Remove Columns that won't be 7used
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
