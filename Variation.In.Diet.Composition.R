
### load packages

library(readxl)
library(tidyverse)
library(vegan)

### Import data

data<-read_excel("SM.xlsx")

### For family level

# Long table to wide table 

sm <- dplyr::select(data, sample, family, clean.reads)
sm<-aggregate(clean.reads~sample+family, data=sm, FUN=sum) 

sm<-pivot_wider(sm, 
                names_from = family, 
                values_from = clean.reads, 
                values_fill = 0)
samples <- sm$sample
sm<-dplyr::select(sm, -sample)
sm[sm > 0] <- 1 # change number of reads to binary (presence/absence)
sm<-cbind(samples, sm)
colnames(sm)[1] <- "sample"

sm_sam <- dplyr::group_by(data, sample) %>% filter(row_number()==1) %>% 
  dplyr::select(sample, site, month.cat, "C/R", ring, age.cat, sex)
sm<- merge(sm_sam, sm, by.x = "sample", all.x = FALSE, all.y = FALSE)

# PERMANOVA

sm <- filter(sm, sex != "unk")
PRM.FAM <- adonis2(sm[,8:ncol(sm)] ~ sm$month.cat + sm$sex + sm$age.cat, permutations = 9999)

TukeyHSD(betadisper(dist(sm[,8:ncol(sm)]), sm$month.cat, type = c("median")))
anova(betadisper(dist(sm[,8:ncol(sm)]), sm$month.cat, type = c("median")))

# SIMPER: Percentage of similarity

simper(sm[,8:ncol(sm)],sm$month.cat, permutations = 9999)

### For OTU level

# Long table to wide table 

sm <- dplyr::select(data, sample, otu.id, clean.reads)
sm<-aggregate(clean.reads~sample+otu.id, data=sm, FUN=sum) 

sm<-pivot_wider(sm, 
                names_from = otu.id, 
                values_from = clean.reads, 
                values_fill = 0)
samples <- sm$sample
sm<-dplyr::select(sm, -sample)
sm[sm > 0] <- 1 # change number of reads to binary (presence/absence)
sm<-cbind(samples, sm)
colnames(sm)[1] <- "sample"

sm_sam <- group_by(data, sample) %>% filter(row_number()==1) %>% 
  dplyr::select(sample, site, month.cat, "C/R", ring, age.cat, sex)
sm<- merge(sm_sam, sm, by.x = "sample", all.x = FALSE, all.y = FALSE)

# PERMANOVA

sm <- filter(sm, sex != "unk")
PER.OTU <- adonis2(sm[,8:ncol(sm)] ~ sm$month.cat + sm$sex + sm$age.cat, permutations = 9999)

TukeyHSD(betadisper(dist(sm[,11:ncol(sm)]), sm$month.cat, type = c("median")))
anova(betadisper(dist(sm[,11:ncol(sm)]), sm$month.cat, type = c("median")))

# SIMPER: Percentage of similarity

simper(sm[,8:ncol(sm)],sm$month.cat, permutations = 9999)
