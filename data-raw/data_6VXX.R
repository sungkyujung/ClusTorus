## code to prepare `data_6VXX` dataset goes here
library(bio3d)
pdb <- bio3d::read.pdb("6vxx")
data_6VXX <- bio3d::torsion.pdb(pdb)

usethis::use_data(data_6VXX, overwrite = TRUE)
