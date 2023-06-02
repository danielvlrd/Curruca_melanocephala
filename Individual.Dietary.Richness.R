
### load packages

library(MASS)
library(readxl)
library(tidyverse)

### import data

data<-read_excel("SM.xlsx")

### calculate OTU and Family richness by sample

input<-data[!duplicated(data$sample.new), ]
input<-dplyr::select(input, sample.new, month2, age2, sex, ring)
input$richness.otu <- NA
input$richness.fam <- NA

for (i in 1:nrow(input)){
  
  temp<-data[which(data$sample.new == input$sample.new[i]),]
  input$richness.otu[i]<-length(unique(temp$final_id))
  rm(temp)
}

for (i in 1:nrow(input)){
  
  temp<-data[which(data$sample.new == input$sample.new[i]),]
  input$richness.fam[i]<-length(unique(temp$family))
  rm(temp)
}

### Modelling individual dietary richness for OTU

# fit poisson model

p.model.otu <- glm(richness.otu ~ month2 + age2 + sex, family = 'poisson', data = input)

# check dispersion

p_res <- resid(p.model.otu)
plot(fitted(p.model.otu), p_res, col='steelblue', pch=16,
     xlab='Predicted Offers', ylab='Standardized Residuals', main='Poisson')
abline(0,0)

mu <- predict(p.model.otu, type="response")
z <- ((input$richness.otu - mu)**2 - input$richness.otu)/ (mu * sqrt(2))
zscore <- lm(z ~ 1)
summary(zscore) # It is significant, so fitting with poisson gives over-dispersion

# fit negative binomial regression model

nb.model.otu <- glm.nb(richness.otu ~ month2 + age2 + sex, data = input)
nb_res <- resid(nb.model.otu)
plot(fitted(nb.model.otu), nb_res, col='steelblue', pch=16,
     xlab='Predicted Offers', ylab='Standardized Residuals', main='Negative Binomial')
abline(0,0)

# compare fit from poisson and negative binomial

pchisq(2 * (logLik(nb.model.otu) - logLik(p.model.otu)),
       df = 1, lower.tail = FALSE) # nb gives better fit

# check results of negative binomial model

summary(nb.model.otu)
anova(nb.model.otu)

### Modelling individual dietary richness for Family

# fit poisson model

p.model.fam <- glm(richness.fam ~ month2 + age2 + sex, family = 'poisson', data = input)

# check dispersion

p_res <- resid(p.model.fam)
plot(fitted(p.model.fam), p_res, col='steelblue', pch=16,
     xlab='Predicted Offers', ylab='Standardized Residuals', main='Poisson')
abline(0,0)

mu <- predict(p.model.fam, type="response")
z <- ((input$richness.fam - mu)**2 - input$richness.fam)/ (mu * sqrt(2))
zscore <- lm(z ~ 1)
summary(zscore) # It is significant, so fitting with poisson gives over-dispersion

# fit negative binomial regression model

nb.model.fam <- glm.nb(richness.fam ~ month2 + age2 + sex, data = input)
nb_res <- resid(nb.model.fam)
plot(fitted(nb.model.fam), nb_res, col='steelblue', pch=16,
     xlab='Predicted Offers', ylab='Standardized Residuals', main='Negative Binomial')
abline(0,0)

# compare fit from poisson and negative binomial

pchisq(2 * (logLik(nb.model.fam) - logLik(p.model.fam)), 
       df = 1, lower.tail = FALSE) # nb gives better fit

# check results of negative binomial model

summary(nb.model.fam)
anova(nb.model.fam)
