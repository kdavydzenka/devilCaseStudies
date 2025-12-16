
rm(list = ls())
library(tidyverse)

conditions <- c("age_only", "age_type1", "age_type2", "interaction")
subs_full = c("subsampled", "full")

for (c in conditions) {
  for (sf in subs_full) {
    
    rna_devil_name <- paste0("results/MuscleRNA/",sf,"/devil_", c, ".RDS")
    rna_devil <- readRDS(rna_devil_name)
    
    rna_nebula_name <- paste0("results/MuscleRNA/",sf,"/nebula_", c, ".RDS")
    rna_nebula <- readRDS(rna_nebula_name)
    
    if (any(dim(rna_devil) != dim(rna_nebula))) {
      stop()
    }
  
    rna_nebula$name = rna_devil$name
    saveRDS(rna_nebula, file = rna_nebula_name)
  }
}
