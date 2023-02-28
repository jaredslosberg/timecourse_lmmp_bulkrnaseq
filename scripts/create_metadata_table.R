samples <- read.csv(here("samples.txt"), header = F) 
colnames(samples) <- "sample_name"
data <- samples$sample_name %>% str_split("-|_") %>% do.call(rbind,.) #split at - or _
colnames(data) <- c("age","sex","rep","cell","lane")                  
metadata <- cbind(samples, data)
write.csv(metadata, here("metadata.csv"), quote = F, row.names = F)

metadata_merged <- metadata %>% dplyr::select(-c(sample_name, lane)) %>%
  mutate(sample_name = paste(age, sex,rep, sep = "-")) %>%
  unique() %>%
  relocate(sample_name)

write.csv(metadata_merged, here("metadata_merged.csv"), quote = F, row.names = F)

write_csv(metadata_merged[,"sample_name",drop = F], here("samples_merged.txt"), col_names = F)
