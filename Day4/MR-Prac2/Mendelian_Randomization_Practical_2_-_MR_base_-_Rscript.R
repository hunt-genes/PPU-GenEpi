# Install package (if these packages are already installed in your version of R then you can go down to the 'Load package' section below)
install.packages("devtools")
library(devtools)
install_github("MRCIEU/TwoSampleMR")
install.packages("plyr")
install.packages("ggplot2")

# Load package
library(TwoSampleMR)
library(plyr)
library(ggplot2)

# Set the working directory
setwd("C://Users/User/Desktop/PPU-GenEpi-main/Day4/")

# Load BMI (exposure) data
bmi_exp_data <- read_exposure_data(filename="BMIestimates-BMI_SNPs.txt", sep="\t", snp_col="SNP", beta_col="b", se_col="se", 
	pval_col="p", eaf_col="Freq1.Hapmap", effect_allele_col="A1", other_allele_col="A2", samplesize_col="N")
bmi_exp_data$exposure <- "Body mass index"

# Look at the data and check the dimensions to ensure all SNPs are present
head(bmi_exp_data)
dim(bmi_exp_data)

# Load CHD (outcome) data
chd_out_data <- read_outcome_data(filename="CHDestimates-BMI_SNPs.txt", sep="\t", snp_col="markername", beta_col="beta", 
	se_col="se_dgc", pval_col="p_dgc", eaf_col="effect_allele_freq", effect_allele_col="effect_allele", other_allele_col="noneffect_allele")    
chd_out_data$outcome <- "CHD"

# Look at the data and check the dimensions to ensure all SNPs are present
head(chd_out_data)
dim(chd_out_data)

# Harmonise the CHD and BMI datasets so that the effect alleles are the same.
# This syntax will flip the log odds ratio and effect alleles in the CARDIoGRAM dataset where the effect # alleles are different between  CARDIoGRAMplusC4D and GIANT.
dat <- harmonise_data(bmi_exp_data, chd_out_data)

# If you explore the dataset you'll notice that effect alleles and log odds ratios have been flipped in the # CHD dataset where the effect allele in the CHD dataset was different from the effect allele in the BMI # dataset
head(dat)

# Let's use the MR-Base R package to estimate the causal effect of BMI on CHD
res <- mr(dat,method_list=c("mr_ivw", "mr_simple_median", "mr_weighted_median", "mr_egger_regression"))

# Results from the MR-Base package using various methods including MR-Egger and weighted median sensitivity analyses 
res 

#estimate odds ratio and 95% confidence interval
exp(res$b[1])
exp(res$b[1]-1.96*res$se[1])
exp(res$b[1]+1.96*res$se[1])

# Test the intercept in MR-Egger to see if there is any evidence of directional pleiotropy
egg.int <- mr_pleiotropy_test(dat) # MR-Egger intercept test 
egg.int

# Let's see if there is any heterogeneity between the Wald ratio estimates
het <- mr_heterogeneity(dat, method_list="mr_ivw")
het

# Let's estimate the Wald ratio for each SNP to see if any particular SNP is having a big influence on the overall causal estimate
res_single <- mr_singlesnp(dat, all_method=c("mr_ivw", "mr_egger_regression", "mr_weighted_median"))
res_single

# Create a scatter plot of the SNP-CHD and SNP-BMI associations. 
p1 <- mr_scatter_plot(res, dat)
p1 
#save your plot using the png() function
png("scatter.png")
p1
dev.off()

# Create a funnel plot of the results. 
res_single <- mr_singlesnp(dat)
p2 <- mr_funnel_plot(res_single)
p2
#save your plot using the png() function
png("funnel.png")
p2
dev.off()

# 3. Create a forest plot of the results from the single SNP analysis in the previous section
p3 <- mr_forest_plot(res_single)
p3
#save your plot using the png() function
png("forest.png")
p3
dev.off()
