## code to prepare `protein_6a32` dataset goes here
library(bio3d)

pdb <- read.pdb("6a32")
protein_6a32 <- torsion.pdb(pdb)
usethis::use_data(protein_6a32, overwrite = TRUE, version = 3)
