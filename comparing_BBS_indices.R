library(tidyverse)


b17 <- read.csv('Rosenberg et al annual indices of abundance.csv',stringsAsFactors = F) 



b18a <- read.csv('inde_best_1966-2018_core.csv',stringsAsFactors = F) 
b18 <- b18a %>% filter(Region == "SU1")





