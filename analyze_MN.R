## analyse images

library(tidyverse)
library(bioimagetools)

## define input

main_dir <- "C:/Users/macrh/Dropbox/Studium/Master/Praktikum_de_Stefano/confocal_images/final_analysis/"
output_dir <- file.path(main_dir, "intermediate")
channel_dir <- file.path(main_dir, "raw")
selection_dir <- file.path(main_dir, "mn_masks")

pixel_shift <- 4.16
cutoff_ym <- 50
nrow_img <- 1024
radius_ym <- 2


## processing
source("C:/Users/macrh/repos/MN_analysis/fncts.R")
dir.create(output_dir)
cutoff_px <- cutoff_ym*pixel_shift
mns <- list.files(selection_dir, pattern=".tif$", full.names = T)
radius <- pixel_shift*radius_ym


analyzed_mns <- lapply(mns, measure_mn,
                       radius, pixel_shift, channel_dir, output_dir) %>%
  bind_rows()

write_tsv(analyzed_mns, file=opt$output_file)
