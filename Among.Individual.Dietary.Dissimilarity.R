
### load packages

library(readxl)
library(vegan)
library(tidyverse)

### import data

data<-read_excel("SM.xlsx")

sm <- dplyr::select(data, sample, family, clean.reads)
sm<-aggregate(clean.reads~sample+family, data=sm, FUN=sum) 

### Move from long table to wide table (list to matrix)

sm<-pivot_wider(sm, 
                names_from = family, 
                values_from = clean.reads, 
                values_fill = 0)
which(is.na(sm)) # no NA
samples <- sm$sample
sm<-dplyr::select(sm, -sample)
sm[sm > 0] <- 1
sm<-cbind(samples, sm)
colnames(sm)[1] <- "sample"

### add information

sm_sam <- group_by(data, sample) %>% filter(row_number()==1) %>% 
  dplyr::select(sample, site, month.cat, "C/R", ring, age.cat, sex)
sm<- merge(sm_sam, sm, by.x = "sample", all.x = FALSE, all.y = FALSE)
rm(sm_sam, samples)

### simulate null models and calculate z-scores

meandist <- function(x) mean(vegdist(x, "jaccard"))

## Apr-May

AM<-dplyr::filter(sm, month.cat == "Apr-May") 
AM <- cbind(AM[,1:7], AM[,(7 + which(colSums(AM[,8:ncol(AM)]) > 0))])
Null.AM.r1 <- vegan::nullmodel(AM[,8:ncol(AM)], method="r1")
Null.AM.c0 <- vegan::nullmodel(AM[,8:ncol(AM)], method="c0")
Null.AM.r1 <- simulate(Null.AM.r1, nsim = 9999)
Null.AM.c0 <- simulate(Null.AM.c0, nsim = 9999)
oecosimu(Null.AM.r1, meandist)
oecosimu(Null.AM.c0, meandist)

## Jun-Jul

JJ<-dplyr::filter(sm, month.cat == "Jun-Jul") 
JJ <- cbind(JJ[,1:7], JJ[,(7 + which(colSums(JJ[,8:ncol(JJ)]) > 0))])
Null.JJ.r1 <- vegan::nullmodel(JJ[,8:ncol(JJ)], method="r1")
Null.JJ.c0 <- vegan::nullmodel(JJ[,8:ncol(JJ)], method="c0")
Null.JJ.r1 <- simulate(Null.JJ.r1, nsim = 9999)
Null.JJ.c0 <- simulate(Null.JJ.c0, nsim = 9999)
oecosimu(Null.JJ.r1, meandist)
oecosimu(Null.JJ.c0, meandist)

## Aug-Sep

AS<-dplyr::filter(sm, month.cat == "Aug-Sep") 
AS <- cbind(AS[,1:7], AS[,(7 + which(colSums(AS[,8:ncol(AS)]) > 0))])
Null.AS.r1 <- vegan::nullmodel(AS[,8:ncol(AS)], method="r1")
Null.AS.c0 <- vegan::nullmodel(AS[,8:ncol(AS)], method="c0")
Null.AS.r1 <- simulate(Null.AS.r1, nsim = 9999)
Null.AS.c0 <- simulate(Null.AS.c0, nsim = 9999)
oecosimu(Null.AS.r1, meandist)
oecosimu(Null.AS.c0, meandist)

## Oct-Nov

ON<-dplyr::filter(sm, month.cat == "Oct-Nov") 
ON <- cbind(ON[,1:7], ON[,(7 + which(colSums(ON[,8:ncol(ON)]) > 0))])
Null.ON.r1 <- vegan::nullmodel(ON[,8:ncol(ON)], method="r1")
Null.ON.c0 <- vegan::nullmodel(ON[,8:ncol(ON)], method="c0")
Null.ON.r1 <- simulate(Null.ON.r1, nsim = 9999)
Null.ON.c0 <- simulate(Null.ON.c0, nsim = 9999)
oecosimu(Null.ON.r1, meandist)
oecosimu(Null.ON.c0, meandist)

## Dec-Jan

DJ<-dplyr::filter(sm, month.cat == "Dec-Jan") 
DJ <- cbind(DJ[,1:7], DJ[,(7 + which(colSums(DJ[,8:ncol(DJ)]) > 0))])
Null.DJ.r1 <- vegan::nullmodel(DJ[,8:ncol(DJ)], method="r1")
Null.DJ.c0 <- vegan::nullmodel(DJ[,8:ncol(DJ)], method="c0")
Null.DJ.r1 <- simulate(Null.DJ.r1, nsim = 9999)
Null.DJ.c0 <- simulate(Null.DJ.c0, nsim = 9999)
oecosimu(Null.DJ.r1, meandist)
oecosimu(Null.DJ.c0, meandist)


## Feb-Mar

FM<-dplyr::filter(sm, month.cat == "Feb-Mar") 
FM <- cbind(FM[,1:7], FM[,(7 + which(colSums(FM[,8:ncol(FM)]) > 0))])
Null.FM.r1 <- vegan::nullmodel(FM[,8:ncol(FM)], method="r1")
Null.FM.c0 <- vegan::nullmodel(FM[,8:ncol(FM)], method="c0")
Null.FM.r1 <- simulate(Null.FM.r1, nsim = 9999)
Null.FM.c0 <- simulate(Null.FM.c0, nsim = 9999)
oecosimu(Null.FM.r1, meandist)
oecosimu(Null.FM.c0, meandist)

### Checking variation of replacement, richness difference and overlap

### function of replacement, richness difference and overlap

rep <- function(x) {
  mat.b <- ifelse(x > 0, 1, 0)
  a <- mat.b %*% t(mat.b)
  b <- mat.b %*% (1 - t(mat.b))
  c <- (1 - mat.b) %*% t(mat.b)
  min.bc <- pmin(b, c)
  
  repl <- 2 * min.bc
  return(sum(repl)/length(repl))
}


rich <- function(x) {
  mat.b <- ifelse(x > 0, 1, 0)
  a <- mat.b %*% t(mat.b)
  b <- mat.b %*% (1 - t(mat.b))
  c <- (1 - mat.b) %*% t(mat.b)
  min.bc <- pmin(b, c)
  
  rich <- abs(b - c)
  return(sum(rich)/length(rich))
}

over <- function(x) {
  mat.b <- ifelse(x > 0, 1, 0)
  a <- mat.b %*% t(mat.b)
  
  return(sum(a)/length(a))
}

### Comparison of observed with null models



beta.div.part.r1 <- data.frame(month = c("Apr-May", "Jun-Jul", "Aug-Sep",
                                      "Oct-Nov", "Dec-Jan", "Feb-Mar"),
                               z.rep = c(oecosimu(Null.AM.r1, rep)[[2]]$z,
                                         oecosimu(Null.JJ.r1, rep)[[2]]$z,
                                         oecosimu(Null.AS.r1, rep)[[2]]$z,
                                         oecosimu(Null.ON.r1, rep)[[2]]$z,
                                         oecosimu(Null.DJ.r1, rep)[[2]]$z,
                                         oecosimu(Null.FM.r1, rep)[[2]]$z),
                               z.rich = c(oecosimu(Null.AM.r1, rich)[[2]]$z,
                                          oecosimu(Null.JJ.r1, rich)[[2]]$z,
                                          oecosimu(Null.AS.r1, rich)[[2]]$z,
                                          oecosimu(Null.ON.r1, rich)[[2]]$z,
                                          oecosimu(Null.DJ.r1, rich)[[2]]$z,
                                          oecosimu(Null.FM.r1, rich)[[2]]$z),
                               z.over = c(oecosimu(Null.AM.r1, over)[[2]]$z,
                                          oecosimu(Null.JJ.r1, over)[[2]]$z,
                                          oecosimu(Null.AS.r1, over)[[2]]$z,
                                          oecosimu(Null.ON.r1, over)[[2]]$z,
                                          oecosimu(Null.DJ.r1, over)[[2]]$z,
                                          oecosimu(Null.FM.r1, over)[[2]]$z)
                            )

beta.div.part.c0 <- data.frame(month = c("Apr-May", "Jun-Jul", "Aug-Sep",
                                         "Oct-Nov", "Dec-Jan", "Feb-Mar"),
                               z.rep = c(oecosimu(Null.AM.c0, rep)[[2]]$z,
                                         oecosimu(Null.JJ.c0, rep)[[2]]$z,
                                         oecosimu(Null.AS.c0, rep)[[2]]$z,
                                         oecosimu(Null.ON.c0, rep)[[2]]$z,
                                         oecosimu(Null.DJ.c0, rep)[[2]]$z,
                                         oecosimu(Null.FM.c0, rep)[[2]]$z),
                               z.rich = c(oecosimu(Null.AM.c0, rich)[[2]]$z,
                                          oecosimu(Null.JJ.c0, rich)[[2]]$z,
                                          oecosimu(Null.AS.c0, rich)[[2]]$z,
                                          oecosimu(Null.ON.c0, rich)[[2]]$z,
                                          oecosimu(Null.DJ.c0, rich)[[2]]$z,
                                          oecosimu(Null.FM.c0, rich)[[2]]$z),
                               z.over = c(oecosimu(Null.AM.c0, over)[[2]]$z,
                                          oecosimu(Null.JJ.c0, over)[[2]]$z,
                                          oecosimu(Null.AS.c0, over)[[2]]$z,
                                          oecosimu(Null.ON.c0, over)[[2]]$z,
                                          oecosimu(Null.DJ.c0, over)[[2]]$z,
                                          oecosimu(Null.FM.c0, over)[[2]]$z)
)
