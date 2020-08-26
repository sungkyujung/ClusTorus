## code to prepare `protein_6K7G` dataset goes here
library(bio3d)

pdb <- read.pdb("6K7G")
protein_6K7G <- torsion.pdb(pdb)
usethis::use_data(protein_6K7G, overwrite = TRUE, version = 3)
