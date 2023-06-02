
### load packages

library(readxl)
library(vegan)
library(tidyverse)

### import data

data<-read_excel("SM.xlsx")

sm <- dplyr::select(data, sample.new, family, final.reads)
sm<-aggregate(final.reads~sample.new+family, data=sm, FUN=sum) 

### Move from long table to wide table (list to matrix)

sm<-pivot_wider(sm, 
                names_from = family, 
                values_from = final.reads, 
                values_fill = 0)
which(is.na(sm)) # no NA
samples <- sm$sample.new
sm<-dplyr::select(sm, -sample.new)
sm[sm > 0] <- 1
sm<-cbind(samples, sm)
colnames(sm)[1] <- "sample.new"

### add information

sm_sam <- group_by(data, sample.new) %>% filter(row_number()==1) %>% 
dplyr::select(sample.new, Local, year, data, month, month2, "N/R", ring, age2, sex)
sm<- merge(sm_sam, sm, by.x = "sample.new", all.x = FALSE, all.y = FALSE)

### simulate null models and calculate z-scores

meandist <- function(x) mean(vegdist(x, "jaccard"))

## Apr-May

AM<-dplyr::filter(sm, month2 == "Apr-May") 
AM <- cbind(AM[,1:10], AM[,(10 + which(colSums(AM[,11:ncol(AM)]) > 0))])
Null.AM.r1 <- vegan::nullmodel(AM[,11:ncol(AM)], method="r1")
Null.AM.c0 <- vegan::nullmodel(AM[,11:ncol(AM)], method="c0")
Null.AM.r1 <- simulate(Null.AM.r1, nsim = 9999)
Null.AM.c0 <- simulate(Null.AM.c0, nsim = 9999)
oecosimu(Null.AM.r1, meandist)
oecosimu(Null.AM.c0, meandist)

## Jun-Jul

JJ<-dplyr::filter(sm, month2 == "Jun-Jul") 
JJ <- cbind(JJ[,1:10], JJ[,(10 + which(colSums(JJ[,11:ncol(JJ)]) > 0))])
Null.JJ.r1 <- vegan::nullmodel(JJ[,11:ncol(JJ)], method="r1")
Null.JJ.c0 <- vegan::nullmodel(JJ[,11:ncol(JJ)], method="c0")
Null.JJ.r1 <- simulate(Null.JJ.r1, nsim = 9999)
Null.JJ.c0 <- simulate(Null.JJ.c0, nsim = 9999)
oecosimu(Null.JJ.r1, meandist)
oecosimu(Null.JJ.c0, meandist)

## Aug-Sep

AS<-dplyr::filter(sm, month2 == "Aug-Sep") 
AS <- cbind(AS[,1:10], AS[,(10 + which(colSums(AS[,11:ncol(AS)]) > 0))])
Null.AS.r1 <- vegan::nullmodel(AS[,11:ncol(AS)], method="r1")
Null.AS.c0 <- vegan::nullmodel(AS[,11:ncol(AS)], method="c0")
Null.AS.r1 <- simulate(Null.AS.r1, nsim = 9999)
Null.AS.c0 <- simulate(Null.AS.c0, nsim = 9999)
oecosimu(Null.AS.r1, meandist)
oecosimu(Null.AS.c0, meandist)

## Oct-Nov

ON<-dplyr::filter(sm, month2 == "Oct-Nov") 
ON <- cbind(ON[,1:10], ON[,(10 + which(colSums(ON[,11:ncol(ON)]) > 0))])
Null.ON.r1 <- vegan::nullmodel(ON[,11:ncol(ON)], method="r1")
Null.ON.c0 <- vegan::nullmodel(ON[,11:ncol(ON)], method="c0")
Null.ON.r1 <- simulate(Null.ON.r1, nsim = 9999)
Null.ON.c0 <- simulate(Null.ON.c0, nsim = 9999)
oecosimu(Null.ON.r1, meandist)
oecosimu(Null.ON.c0, meandist)

## Dec-Jan

DJ<-dplyr::filter(sm, month2 == "Dec-Jan") 
DJ <- cbind(DJ[,1:10], DJ[,(10 + which(colSums(DJ[,11:ncol(DJ)]) > 0))])
Null.DJ.r1 <- vegan::nullmodel(DJ[,11:ncol(DJ)], method="r1")
Null.DJ.c0 <- vegan::nullmodel(DJ[,11:ncol(DJ)], method="c0")
Null.DJ.r1 <- simulate(Null.DJ.r1, nsim = 9999)
Null.DJ.c0 <- simulate(Null.DJ.c0, nsim = 9999)
oecosimu(Null.DJ.r1, meandist)
oecosimu(Null.DJ.c0, meandist)


## Feb-Mar

FM<-dplyr::filter(sm, month2 == "Feb-Mar") 
FM <- cbind(FM[,1:10], FM[,(10 + which(colSums(FM[,11:ncol(FM)]) > 0))])
Null.FM.r1 <- vegan::nullmodel(FM[,11:ncol(FM)], method="r1")
Null.FM.c0 <- vegan::nullmodel(FM[,11:ncol(FM)], method="c0")
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
