## code to prepare `data_6VXX` dataset goes here
library(bio3d)
library(tidyverse)
pdb <- read.pdb("6vxx")

data_6VXX <- torsion.pdb(pdb)
data_6VXX <- data_6VXX$tbl / 180 * pi

usethis::use_data(data_6VXX, overwrite = TRUE)
