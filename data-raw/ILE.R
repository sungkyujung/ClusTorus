## prepare `ILE` dataset goes here
library(bio3d)
library(tidyverse)
devtools::load_all()

# high quality protein data obtained from PISCES server
# resolution : 1.6A(angstrom) or better
# R-factor : 0.22 or better
# Sequence percentage identity: <= 25%
path <- paste0(getwd(), "/data-raw/PISCES_data.txt")
protein <- read.table(path, header = TRUE)

# function for obtaining diheral angles of specified amino acid
# In Mardia (2012), the article uses ILE only data.
amino.data <- function(data, amino, THRES = 1000){
  IDs <- data$IDs
  ang <- data.frame()
  n <- 0
  for (ID in IDs){
    id <- substr(ID, 1, 4)
    type <-substr(ID, 5, 5)

    tbl <- bio3d::torsion.pdb(bio3d::read.pdb(id))$tbl
    tbl <- data.frame(tbl, id = rownames(tbl)) %>%
      separate(id,into = c(',','position','type','rest'))

    tbl[,1:7] <- tbl[,1:7]/180*pi
    tbl[,1:7] <- on.torus(tbl[,1:7])

    tbl <- tbl[(tbl$type == type) & (tbl$rest == amino), ]
    ang <- rbind(ang, tbl)
    n <- n + 1
    if (n == THRES) {break}
  }

  ang
}
ILE <- amino.data(protein, "ILE")
ILE <- na.omit(ILE[, 1:4])
usethis::use_data(ILE, overwrite = TRUE)
