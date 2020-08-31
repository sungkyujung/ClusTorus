## code to prepare `SARS_CoV_2` dataset goes here
library(bio3d)

pdb <- read.pdb("6vxx")
SARS_CoV_2 <- torsion.pdb(pdb)
usethis::use_data(SARS_CoV_2, overwrite = TRUE)
