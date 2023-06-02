
### load packages

library(iNEXT)
library(tidyverse)
library(readxl)

### import data

data<-read_excel("SM.xlsx")

### Estimate population's total diet richness 

# For OTU

OTU<-as.data.frame(with(data, tapply(sample.new, final_id, 
                                        FUN = function(x) length(unique(x)))))
names(OTU)[1] <- "N_of_samples"

OTU <- list(Diet_Otu = c(length(unique(data$sample.new)), 
                         sort(OTU$N_of_samples, decreasing = TRUE)))
iNEXT(OTU, datatype = "incidence_freq")

# For Family

FAM<-as.data.frame(with(data, tapply(sample.new, family, 
                                        FUN = function(x) length(unique(x)))))
names(FAM)[1] <- "N_of_samples"

FAM <- list(Diet_FAM = c(length(unique(data$sample.new)), 
                         sort(FAM$N_of_samples, decreasing = TRUE)))
iNEXT(FAM, datatype = "incidence_freq")

# For Order

ORD<-as.data.frame(with(data, tapply(sample.new, order, 
                                        FUN = function(x) length(unique(x)))))
names(ORD)[1] <- "N_of_samples"

ORD <- list(Diet_ORD = c(length(unique(data$sample.new)), 
                         sort(ORD$N_of_samples, decreasing = TRUE)))
iNEXT(ORD, datatype = "incidence_freq")

### Comparisons between pair of months, sexes and ages 

### For OTU 

## Between pair of months

OTU.AM<-as.data.frame.table(with(filter(data, month2== "Apr-May"),
            tapply(sample.new, final_id, FUN = function(x) length(unique(x)))))
names(OTU.AM)[2] <- "N_of_samples"
OTU.JJ<-as.data.frame.table(with(filter(data, month2== "Jun-Jul"),
            tapply(sample.new, final_id, FUN = function(x) length(unique(x)))))
names(OTU.JJ)[2] <- "N_of_samples"
OTU.AS<-as.data.frame.table(with(filter(data, month2== "Aug-Sep"),
            tapply(sample.new, final_id, FUN = function(x) length(unique(x)))))
names(OTU.AS)[2] <- "N_of_samples"
OTU.ON<-as.data.frame.table(with(filter(data, month2== "Oct-Nov"), 
            tapply(sample.new, final_id, FUN = function(x) length(unique(x)))))
names(OTU.ON)[2] <- "N_of_samples"
OTU.DJ<-as.data.frame.table(with(filter(data, month2== "Dec-Jan"),
            tapply(sample.new, final_id, FUN = function(x) length(unique(x)))))
names(OTU.DJ)[2] <- "N_of_samples"
OTU.FM<-as.data.frame.table(with(filter(data, month2== "Feb-Mar"), 
            tapply(sample.new, final_id, FUN = function(x) length(unique(x)))))
names(OTU.FM)[2] <- "N_of_samples"

OTU.MONTHS <-list(Apr_May = c(nrow(unique(data[which(data$month2 == "Apr-May"), "sample.new"])), 
                             sort(OTU.AM$N_of_samples, decreasing = TRUE)),
                  Jun_Jul = c(nrow(unique(data[which(data$month2 == "Jun-Jul"), "sample.new"])), 
                              sort(OTU.JJ$N_of_samples, decreasing = TRUE)),
                  Aug_Sep = c(nrow(unique(data[which(data$month2 == "Aug-Sep"), "sample.new"])), 
                              sort(OTU.AS$N_of_samples, decreasing = TRUE)),
                  Oct_Nov = c(nrow(unique(data[which(data$month2 == "Oct-Nov"), "sample.new"])), 
                              sort(OTU.ON$N_of_samples, decreasing = TRUE)),
                  Dec_Jan = c(nrow(unique(data[which(data$month2 == "Dec-Jan"), "sample.new"])), 
                              sort(OTU.DJ$N_of_samples, decreasing = TRUE)),
                  Feb_Mar = c(nrow(unique(data[which(data$month2 == "Feb-Mar"), "sample.new"])), 
                              sort(OTU.FM$N_of_samples, decreasing = TRUE)))
rm(OTU.AM, OTU.JJ, OTU.AS, OTU.ON, OTU.DJ, OTU.FM)

iNEXT(OTU.MONTHS, datatype = "incidence_freq") 
# Lowest sample Coverage at two time the observed t (Oct-Nov) is 0.646
estimateD(OTU.MONTHS, datatype = "incidence_freq", base = "coverage", level = 0.646, conf = 0.95) 

## Between sexes

OTU.M<-as.data.frame.table(with(filter(data, sex== "M"),
            tapply(sample.new, final_id, FUN = function(x) length(unique(x)))))
names(OTU.M)[2] <- "N_of_samples"
OTU.F<-as.data.frame.table(with(filter(data, sex== "F"), 
            tapply(sample.new, final_id, FUN = function(x) length(unique(x)))))
names(OTU.F)[2] <- "N_of_samples"

OTU.SEX <-list(Male = c(nrow(unique(data[which(data$sex == "M"), "sample.new"])), 
                              sort(OTU.M$N_of_samples, decreasing = TRUE)),
                  Female = c(nrow(unique(data[which(data$sex == "F"), "sample.new"])), 
                              sort(OTU.F$N_of_samples, decreasing = TRUE)))
rm(OTU.M, OTU.F)

iNEXT(OTU.SEX, datatype = "incidence_freq") 
# Lowest sample coverage at two time the observed t (Female) is 0.75
estimateD(OTU.SEX, datatype = "incidence_freq", base = "coverage", level = 0.75, conf = 0.95) 

## Between age categories

OTU.1<-as.data.frame.table(with(filter(data, age2== "1cy"), 
            tapply(sample.new, final_id, FUN = function(x) length(unique(x)))))
names(OTU.1)[2] <- "N_of_samples"
OTU.2<-as.data.frame.table(with(filter(data, age2== "2cy+"),
            tapply(sample.new, final_id, FUN = function(x) length(unique(x)))))
names(OTU.2)[2] <- "N_of_samples"

OTU.AGE <-list(year1 = c(nrow(unique(data[which(data$age2 == "1cy"), "sample.new"])), 
                        sort(OTU.1$N_of_samples, decreasing = TRUE)),
               year2 = c(nrow(unique(data[which(data$age2 == "2cy+"), "sample.new"])), 
                          sort(OTU.2$N_of_samples, decreasing = TRUE)))
rm(OTU.1, OTU.2)

iNEXT(OTU.AGE, datatype = "incidence_freq") 
# Lowest sample coverage at two time the observed t (year2) is 0.747
estimateD(OTU.AGE, datatype = "incidence_freq", base = "coverage", level = 0.747, conf = 0.95)

### For Families 

## Between pair of months

FAM.AM<-as.data.frame.table(with(filter(data, month2== "Apr-May"), 
              tapply(sample.new, family, FUN = function(x) length(unique(x)))))
names(FAM.AM)[2] <- "N_of_samples"
FAM.JJ<-as.data.frame.table(with(filter(data, month2== "Jun-Jul"),
              tapply(sample.new, family, FUN = function(x) length(unique(x)))))
names(FAM.JJ)[2] <- "N_of_samples"
FAM.AS<-as.data.frame.table(with(filter(data, month2== "Aug-Sep"),
              tapply(sample.new, family, FUN = function(x) length(unique(x)))))
names(FAM.AS)[2] <- "N_of_samples"
FAM.ON<-as.data.frame.table(with(filter(data, month2== "Oct-Nov"), 
              tapply(sample.new, family, FUN = function(x) length(unique(x)))))
names(FAM.ON)[2] <- "N_of_samples"
FAM.DJ<-as.data.frame.table(with(filter(data, month2== "Dec-Jan"), 
              tapply(sample.new, family, FUN = function(x) length(unique(x)))))
names(FAM.DJ)[2] <- "N_of_samples"
FAM.FM<-as.data.frame.table(with(filter(data, month2== "Feb-Mar"), 
              tapply(sample.new, family, FUN = function(x) length(unique(x)))))
names(FAM.FM)[2] <- "N_of_samples"

FAM.MONTHS <-list(Apr_May = c(nrow(unique(data[which(data$month2 == "Apr-May"), "sample.new"])), 
                              sort(FAM.AM$N_of_samples, decreasing = TRUE)),
                  Jun_Jul = c(nrow(unique(data[which(data$month2 == "Jun-Jul"), "sample.new"])), 
                              sort(FAM.JJ$N_of_samples, decreasing = TRUE)),
                  Aug_Sep = c(nrow(unique(data[which(data$month2 == "Aug-Sep"), "sample.new"])), 
                              sort(FAM.AS$N_of_samples, decreasing = TRUE)),
                  Oct_Nov = c(nrow(unique(data[which(data$month2 == "Oct-Nov"), "sample.new"])), 
                              sort(FAM.ON$N_of_samples, decreasing = TRUE)),
                  Dec_Jan = c(nrow(unique(data[which(data$month2 == "Dec-Jan"), "sample.new"])), 
                              sort(FAM.DJ$N_of_samples, decreasing = TRUE)),
                  Feb_Mar = c(nrow(unique(data[which(data$month2 == "Feb-Mar"), "sample.new"])), 
                              sort(FAM.FM$N_of_samples, decreasing = TRUE)))

rm(FAM.AM, FAM.JJ, FAM.AS, FAM.ON, FAM.DJ, FAM.FM)

iNEXT(FAM.MONTHS, datatype = "incidence_freq") 
# Lowest sample Coverage at two time the observed t (Apr-May) is 0.83
estimateD(FAM.MONTHS, datatype = "incidence_freq", base = "coverage", level = 0.83, conf = 0.95) 

## Between sexes

FAM.M<-as.data.frame.table(with(filter(data, sex== "M"), 
              tapply(sample.new, family, FUN = function(x) length(unique(x)))))
names(FAM.M)[2] <- "N_of_samples"
FAM.F<-as.data.frame.table(with(filter(data, sex== "F"), 
              tapply(sample.new, family, FUN = function(x) length(unique(x)))))
names(FAM.F)[2] <- "N_of_samples"

FAM.SEX <-list(Male = c(nrow(unique(data[which(data$sex == "M"), "sample.new"])), 
                        sort(FAM.M$N_of_samples, decreasing = TRUE)),
               Female = c(nrow(unique(data[which(data$sex == "F"), "sample.new"])), 
                          sort(FAM.F$N_of_samples, decreasing = TRUE)))
rm(FAM.M, FAM.F)

iNEXT(FAM.SEX, datatype = "incidence_freq") 
# Lowest sample coverage at two time the observed t (Female) is 0.895
estimateD(FAM.SEX, datatype = "incidence_freq", base = "coverage", level = 0.895, conf = 0.95) 

## Between age categories

FAM.1<-as.data.frame.table(with(filter(data, age2== "1cy"),
              tapply(sample.new, family, FUN = function(x) length(unique(x)))))
names(FAM.1)[2] <- "N_of_samples"
FAM.2<-as.data.frame.table(with(filter(data, age2== "2cy+"), 
              tapply(sample.new, family, FUN = function(x) length(unique(x)))))
names(FAM.2)[2] <- "N_of_samples"

FAM.AGE <-list(year1 = c(nrow(unique(data[which(data$age2 == "1cy"), "sample.new"])), 
                         sort(FAM.1$N_of_samples, decreasing = TRUE)),
               year2 = c(nrow(unique(data[which(data$age2 == "2cy+"), "sample.new"])), 
                         sort(FAM.2$N_of_samples, decreasing = TRUE)))
rm(FAM.1, FAM.2)

iNEXT(FAM.AGE, datatype = "incidence_freq") 
# Lowest sample coverage at two time the observed t (year2) is 0.932
estimateD(FAM.AGE, datatype = "incidence_freq", base = "coverage", level = 0.932, conf = 0.95)

