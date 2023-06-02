
### load packages

library(readxl)
library(tidyverse)
library(vegan)

### Import data

data<-read_excel("SM.xlsx")

### For family level

# Long table to wide table 

sm <- dplyr::select(data, sample.new, family, final.reads)
sm<-aggregate(final.reads~sample.new+family, data=sm, FUN=sum) 

sm<-pivot_wider(sm, 
                names_from = family, 
                values_from = final.reads, 
                values_fill = 0)
samples <- sm$sample.new
sm<-dplyr::select(sm, -sample.new)
sm[sm > 0] <- 1 # change number of reads to binary (presence/absence)
sm<-cbind(samples, sm)
colnames(sm)[1] <- "sample.new"

sm_sam <- group_by(data, sample.new) %>% filter(row_number()==1) %>%
 dplyr::select(sample.new, Local, year, data, month, month2, "N/R", ring, age2, sex)
sm<- merge(sm_sam, sm, by.x = "sample.new", all.x = FALSE, all.y = FALSE)

# PERMANOVA

sm <- filter(sm, sex != "unk")
PRM.FAM <- adonis2(sm[,11:ncol(sm)] ~ sm$month2 + sm$sex + sm$age2, permutations = 9999)

TukeyHSD(betadisper(dist(sm[,11:ncol(sm)]), sm$month2, type = c("median")))
anova(betadisper(dist(sm[,11:ncol(sm)]), sm$month2, type = c("median")))

# SIMPER: Percentage of similarity

simper(sm[,11:ncol(sm)],sm$month2, permutations = 9999)

### For OTU level

# Long table to wide table 

sm <- dplyr::select(data, sample.new, final_id, final.reads)
sm<-aggregate(final.reads~sample.new+final_id, data=sm, FUN=sum) 

sm<-pivot_wider(sm, 
                names_from = final_id, 
                values_from = final.reads, 
                values_fill = 0)
samples <- sm$sample.new
sm<-dplyr::select(sm, -sample.new)
sm[sm > 0] <- 1 # change number of reads to binary (presence/absence)
sm<-cbind(samples, sm)
colnames(sm)[1] <- "sample.new"

sm_sam <- group_by(data, sample.new) %>% filter(row_number()==1) %>% 
dplyr::select(sample.new, Local, year, data, month, month2, "N/R", ring, age2, sex)
sm<- merge(sm_sam, sm, by.x = "sample.new", all.x = FALSE, all.y = FALSE)

# PERMANOVA

sm <- filter(sm, sex != "unk")
PER.OTU <- adonis2(sm[,11:ncol(sm)] ~ sm$month2 + sm$sex + sm$age2, permutations = 9999)

TukeyHSD(betadisper(dist(sm[,11:ncol(sm)]), sm$month2, type = c("median")))
anova(betadisper(dist(sm[,11:ncol(sm)]), sm$month2, type = c("median")))

# SIMPER: Percentage of similarity

simper(sm[,11:ncol(sm)],sm$month2, permutations = 9999)



