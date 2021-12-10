## code to prepare `SARS_CoV_2` dataset goes here
library(bio3d)
library(tidyverse)
pdb <- read.pdb("6vxx")

SARS_CoV_2 <- torsion.pdb(pdb)
SARS_CoV_2 <- SARS_CoV_2$tbl[, 1:2] / 180 * pi
SARS_CoV_2 <- SARS_CoV_2 %>%
  data.frame(id = rownames(.)) %>%
  mutate(id = trimws(id)) %>%
  separate(id, into = c('position', 'type', 'rest')) %>%
  filter(type %in% c("A", "B", "C"))
SARS_CoV_2[, 1:2] <- on.torus(as.matrix(SARS_CoV_2[, 1:2]))
SARS_CoV_2 <- SARS_CoV_2[, c(1, 2, 4)]

usethis::use_data(SARS_CoV_2, overwrite = TRUE)
