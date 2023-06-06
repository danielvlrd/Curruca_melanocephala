
### load packages

library(iNEXT)
library(tidyverse)
library(readxl)

### import data

data<-read_excel("SM.xlsx")

### Estimate population's total diet richness 

# For OTU

OTU<-as.data.frame(with(data, tapply(sample, otu.id, 
                                        FUN = function(x) length(unique(x)))))
names(OTU)[1] <- "N_of_samples"

OTU <- list(Diet_Otu = c(length(unique(data$sample)), 
                         sort(OTU$N_of_samples, decreasing = TRUE)))
iNEXT(OTU, datatype = "incidence_freq", nboot = 9999)

# For Family

FAM<-as.data.frame(with(data, tapply(sample, family, 
                                        FUN = function(x) length(unique(x)))))
names(FAM)[1] <- "N_of_samples"

FAM <- list(Diet_FAM = c(length(unique(data$sample)), 
                         sort(FAM$N_of_samples, decreasing = TRUE)))
iNEXT(FAM, datatype = "incidence_freq", nboot = 9999)

# For Order

ORD<-as.data.frame(with(data, tapply(sample, order, 
                                        FUN = function(x) length(unique(x)))))
names(ORD)[1] <- "N_of_samples"

ORD <- list(Diet_ORD = c(length(unique(data$sample)), 
                         sort(ORD$N_of_samples, decreasing = TRUE)))
iNEXT(ORD, datatype = "incidence_freq", nboot = 9999)

### Comparisons between pair of months, sexes and ages 

### For OTU 

## Between pair of months

OTU.AM<-as.data.frame.table(with(filter(data, month.cat== "Apr-May"), tapply(sample, otu.id, FUN = function(x) length(unique(x)))))
names(OTU.AM)[2] <- "N_of_samples"
OTU.JJ<-as.data.frame.table(with(filter(data, month.cat== "Jun-Jul"), tapply(sample, otu.id, FUN = function(x) length(unique(x)))))
names(OTU.JJ)[2] <- "N_of_samples"
OTU.AS<-as.data.frame.table(with(filter(data, month.cat== "Aug-Sep"), tapply(sample, otu.id, FUN = function(x) length(unique(x)))))
names(OTU.AS)[2] <- "N_of_samples"
OTU.ON<-as.data.frame.table(with(filter(data, month.cat== "Oct-Nov"), tapply(sample, otu.id, FUN = function(x) length(unique(x)))))
names(OTU.ON)[2] <- "N_of_samples"
OTU.DJ<-as.data.frame.table(with(filter(data, month.cat== "Dec-Jan"), tapply(sample, otu.id, FUN = function(x) length(unique(x)))))
names(OTU.DJ)[2] <- "N_of_samples"
OTU.FM<-as.data.frame.table(with(filter(data, month.cat== "Feb-Mar"), tapply(sample, otu.id, FUN = function(x) length(unique(x)))))
names(OTU.FM)[2] <- "N_of_samples"

OTU.MONTHS <-list(Apr_May = c(nrow(unique(data[which(data$month.cat == "Apr-May"), "sample"])), 
                             sort(OTU.AM$N_of_samples, decreasing = TRUE)),
                  Jun_Jul = c(nrow(unique(data[which(data$month.cat == "Jun-Jul"), "sample"])), 
                              sort(OTU.JJ$N_of_samples, decreasing = TRUE)),
                  Aug_Sep = c(nrow(unique(data[which(data$month.cat == "Aug-Sep"), "sample"])), 
                              sort(OTU.AS$N_of_samples, decreasing = TRUE)),
                  Oct_Nov = c(nrow(unique(data[which(data$month.cat == "Oct-Nov"), "sample"])), 
                              sort(OTU.ON$N_of_samples, decreasing = TRUE)),
                  Dec_Jan = c(nrow(unique(data[which(data$month.cat == "Dec-Jan"), "sample"])), 
                              sort(OTU.DJ$N_of_samples, decreasing = TRUE)),
                  Feb_Mar = c(nrow(unique(data[which(data$month.cat == "Feb-Mar"), "sample"])), 
                              sort(OTU.FM$N_of_samples, decreasing = TRUE)))
rm(OTU.AM, OTU.JJ, OTU.AS, OTU.ON, OTU.DJ, OTU.FM)

iNEXT(OTU.MONTHS, datatype = "incidence_freq", nboot = 9999) 
# Lowest sample Coverage at two time the observed t (Oct-Nov) is 0.646
OTU_Months <- estimateD(OTU.MONTHS, datatype = "incidence_freq", base = "coverage", level = 0.64, conf = 0.95, nboot = 9999, q = 0) 

## Between sexes

OTU.M<-as.data.frame.table(with(filter(data, sex== "M"), tapply(sample, otu.id, FUN = function(x) length(unique(x)))))
names(OTU.M)[2] <- "N_of_samples"
OTU.F<-as.data.frame.table(with(filter(data, sex== "F"), tapply(sample, otu.id, FUN = function(x) length(unique(x)))))
names(OTU.F)[2] <- "N_of_samples"

OTU.SEX <-list(Male = c(nrow(unique(data[which(data$sex == "M"), "sample"])), 
                              sort(OTU.M$N_of_samples, decreasing = TRUE)),
                  Female = c(nrow(unique(data[which(data$sex == "F"), "sample"])), 
                              sort(OTU.F$N_of_samples, decreasing = TRUE)))
rm(OTU.M, OTU.F)

iNEXT(OTU.SEX, datatype = "incidence_freq", nboot = 9999) 
# Lowest sample coverage at two time the observed t (Female) is 0.75
OTU_Sex <- estimateD(OTU.SEX, datatype = "incidence_freq", base = "coverage", level = 0.75, conf = 0.95) 


## Between age categories

OTU.1<-as.data.frame.table(with(filter(data, age.cat== "1cy"), tapply(sample, otu.id, FUN = function(x) length(unique(x)))))
names(OTU.1)[2] <- "N_of_samples"
OTU.2<-as.data.frame.table(with(filter(data, age.cat== "2cy+"), tapply(sample, otu.id, FUN = function(x) length(unique(x)))))
names(OTU.2)[2] <- "N_of_samples"

OTU.AGE <-list(year1 = c(nrow(unique(data[which(data$age.cat == "1cy"), "sample"])), 
                        sort(OTU.1$N_of_samples, decreasing = TRUE)),
               year2 = c(nrow(unique(data[which(data$age.cat == "2cy+"), "sample"])), 
                          sort(OTU.2$N_of_samples, decreasing = TRUE)))
rm(OTU.1, OTU.2)

iNEXT(OTU.AGE, datatype = "incidence_freq", nboot = 9999) 
# Lowest sample coverage at two time the observed t (year2) is 0.747
OTU_Age <- estimateD(OTU.AGE, datatype = "incidence_freq", base = "coverage", level = 0.747, conf = 0.95, nboot = 999)

### For Families 

## Between pair of months

FAM.AM<-as.data.frame.table(with(filter(data, month.cat== "Apr-May"), tapply(sample, family, FUN = function(x) length(unique(x)))))
names(FAM.AM)[2] <- "N_of_samples"
FAM.JJ<-as.data.frame.table(with(filter(data, month.cat== "Jun-Jul"), tapply(sample, family, FUN = function(x) length(unique(x)))))
names(FAM.JJ)[2] <- "N_of_samples"
FAM.AS<-as.data.frame.table(with(filter(data, month.cat== "Aug-Sep"), tapply(sample, family, FUN = function(x) length(unique(x)))))
names(FAM.AS)[2] <- "N_of_samples"
FAM.ON<-as.data.frame.table(with(filter(data, month.cat== "Oct-Nov"), tapply(sample, family, FUN = function(x) length(unique(x)))))
names(FAM.ON)[2] <- "N_of_samples"
FAM.DJ<-as.data.frame.table(with(filter(data, month.cat== "Dec-Jan"), tapply(sample, family, FUN = function(x) length(unique(x)))))
names(FAM.DJ)[2] <- "N_of_samples"
FAM.FM<-as.data.frame.table(with(filter(data, month.cat== "Feb-Mar"), tapply(sample, family, FUN = function(x) length(unique(x)))))
names(FAM.FM)[2] <- "N_of_samples"

FAM.MONTHS <-list(Apr_May = c(nrow(unique(data[which(data$month.cat == "Apr-May"), "sample"])), 
                              sort(FAM.AM$N_of_samples, decreasing = TRUE)),
                  Jun_Jul = c(nrow(unique(data[which(data$month.cat == "Jun-Jul"), "sample"])), 
                              sort(FAM.JJ$N_of_samples, decreasing = TRUE)),
                  Aug_Sep = c(nrow(unique(data[which(data$month.cat == "Aug-Sep"), "sample"])), 
                              sort(FAM.AS$N_of_samples, decreasing = TRUE)),
                  Oct_Nov = c(nrow(unique(data[which(data$month.cat == "Oct-Nov"), "sample"])), 
                              sort(FAM.ON$N_of_samples, decreasing = TRUE)),
                  Dec_Jan = c(nrow(unique(data[which(data$month.cat == "Dec-Jan"), "sample"])), 
                              sort(FAM.DJ$N_of_samples, decreasing = TRUE)),
                  Feb_Mar = c(nrow(unique(data[which(data$month.cat == "Feb-Mar"), "sample"])), 
                              sort(FAM.FM$N_of_samples, decreasing = TRUE)))

rm(FAM.AM, FAM.JJ, FAM.AS, FAM.ON, FAM.DJ, FAM.FM)

iNEXT(FAM.MONTHS, datatype = "incidence_freq", nboot = 9999) 
# Lowest sample Coverage at two time the observed t (Oct-Nov) is 0.95
FAM_Months <- estimateD(FAM.MONTHS, datatype = "incidence_freq", 
                        base = "coverage", level = 0.90, conf = 0.95, nboot = 9999, q = 0) 

## Between sexes

FAM.M<-as.data.frame.table(with(filter(data, sex== "M"), tapply(sample, family, FUN = function(x) length(unique(x)))))
names(FAM.M)[2] <- "N_of_samples"
FAM.F<-as.data.frame.table(with(filter(data, sex== "F"), tapply(sample, family, FUN = function(x) length(unique(x)))))
names(FAM.F)[2] <- "N_of_samples"

FAM.SEX <-list(Male = c(nrow(unique(data[which(data$sex == "M"), "sample"])), 
                        sort(FAM.M$N_of_samples, decreasing = TRUE)),
               Female = c(nrow(unique(data[which(data$sex == "F"), "sample"])), 
                          sort(FAM.F$N_of_samples, decreasing = TRUE)))
rm(FAM.M, FAM.F)

iNEXT(FAM.SEX, datatype = "incidence_freq", nboot = 9999) 
# Lowest sample coverage at two time the observed t (Female) is 0.895
FAM_Sex <-estimateD(FAM.SEX, datatype = "incidence_freq", base = "coverage", level = 0.895, conf = 0.95) 

## Between age categories

FAM.1<-as.data.frame.table(with(filter(data, age.cat== "1cy"), tapply(sample, family, FUN = function(x) length(unique(x)))))
names(FAM.1)[2] <- "N_of_samples"
FAM.2<-as.data.frame.table(with(filter(data, age.cat== "2cy+"), tapply(sample, family, FUN = function(x) length(unique(x)))))
names(FAM.2)[2] <- "N_of_samples"

FAM.AGE <-list(year1 = c(nrow(unique(data[which(data$age.cat == "1cy"), "sample"])), 
                         sort(FAM.1$N_of_samples, decreasing = TRUE)),
               year2 = c(nrow(unique(data[which(data$age.cat == "2cy+"), "sample"])), 
                         sort(FAM.2$N_of_samples, decreasing = TRUE)))
rm(FAM.1, FAM.2)

iNEXT(FAM.AGE, datatype = "incidence_freq", nboot = 9999) 
# Lowest sample coverage at two time the observed t (year2) is 0.932
FAM_Age <- estimateD(FAM.AGE, datatype = "incidence_freq", base = "coverage", level = 0.932, conf = 0.95)

