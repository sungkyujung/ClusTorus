## code to prepare `protein_4q21` dataset goes here
library(bio3d)

pdb <- read.pdb("4q21")
protein_4q21 <- torsion.pdb(pdb)
usethis::use_data(protein_4q21, overwrite = TRUE, version = 3)
