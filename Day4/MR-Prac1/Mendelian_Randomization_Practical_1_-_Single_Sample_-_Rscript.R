#Install dependencies
install.packages(c('ggplot2', 'scales', 'AER', 'devtools', 'dlyr', 
  'data.table', 'stringr', 'optparse'))

library(devtools)
install_github("MRCIEU/TwoSampleMR")

# Read in the dataset 
example <- read.table("data.txt", header=T)
attach(example)

# Look at the data 
head(example) 
summary(example)

# Units: SNP (0,1,2), CRP mmol/L, SBP mmHg, Income per $10,000, HDL mmol/L

# Run observational OLS (ordinary least squares) regression for BP & CRP
summary(lm(SBP~CRP)) 

# Plot the observational association between BP and CRP
plot(CRP,SBP)
abline(lm(SBP~CRP),col="red")

# Observational OLS regression of CRP on CRP SNP
summary(lm(CRP~rs3091244))

# Plot the relationship between CRP and rs3091244
plot(rs3091244, CRP)
abline(lm(CRP~rs3091244),col="red")

# Confounders
summary(lm(SBP~INCOME))
summary(lm(CRP~INCOME))
summary(lm(INCOME~rs3091244))

summary(lm(SBP~HDL))
summary(lm(CRP~HDL))
summary(lm(HDL~rs3091244))

# Run a covariate-adjusted model for the association between CRP & BP
summary(lm(SBP~CRP))
summary(lm(SBP~CRP+INCOME+HDL))

# OLS regression of CRP on CRP SNP
summary(lm(CRP~rs3091244))
# OLS regression of BP on CRP SNP
summary(lm(SBP~rs3091244))

# Observational OLS regression
summary(lm(SBP~CRP))

# Call the AER library to run TSLS (if the AER package has been installed)
library(AER)

# If AER has not been installed, run the command below first:
install.packages("AER")

# General format for TSLS command: 
# summary(ivreg(Outcome~Exposure | Instrument))

# TSLS regression
summary(ivreg(SBP~CRP | rs3091244))

# Regress the exposure (CRP) on the instrument (rs3091244)
First_Stage <- lm(CRP~rs3091244)

# Create predicted CRP values, from first-stage regression
Pred_CRP <- predict(First_Stage)

table(Pred_CRP)
plot(Pred_CRP, CRP)
abline(lm(CRP~ Pred_CRP), col="red")

# Second stage regression
Second_Stage <- lm(SBP~Pred_CRP)

# Look at the results:
summary(Second_Stage)

#Look at F-stat from the first-stage linear regression
summary(lm(CRP~rs3091244))

#Look at F-stat from 'diagnostics' by AER package
summary(ivreg(SBP~CRP | rs3091244), diagnostics=T)

