# DAY 2 - Quality Control and Testing Association for Genome-Wide Association Studies

Please fill the following form during the exercise today: 
  [https://nettskjema.no/a/284234](https://nettskjema.no/a/284234)

Thank you!


### FOCUS and LEARNING GOALS
> The aim for this session is to get familiar with the different data formats and 
> the quality control needed in preparation for running genetic analysis.
> The initial QC of a data‐set is very important for the quality of your down stream
> analysis and is the first step in any analysis. A take home message from this session
> is the difference in results generated from data before and after quality control.

Running a GWAS will often be the second step after 1) quality control (QC) of genotype+imputation data and 2) curation of the phenotype. This type of analysis is often referred to as “hypothesis generating” because it yields a long list of potentially disease-causing genetic loci which can be further studied via functional studies. 

**Suggested reading:**
* [Stephen Turner, et al. Quality Control Procedures for Genome Wide Association Studies, Curr Protoc Hum Genet. 2011](https://pubmed.ncbi.nlm.nih.gov/21234875/)
* [Yik Y. Teoa, Common statistical issues in genome‐wide association studies: a review on power, data quality control, genotype calling and population structure, Current Opinionin Lipidology 2008, 19:133–143](https://www.biorxiv.org/content/10.1101/583278v1)
* [Anderson et.al. Data quality control in genetic case‐control association studies. Nat Protoc. 2010 Sep; 5 (9): 1564‐73.](https://pubmed.ncbi.nlm.nih.gov/21085122/)

We will use Plink to run a GWAS in this practical:

We have already downloaded it on your PC today but you can also find it here if you need it one day
**PLINK:** [https://www.cog-genomics.org/plink/](https://www.cog-genomics.org/plink/)
PLINK is an open source and free toolset developed for the analysis of whole genome data.
It run a range of basic large scale analysis focused on genotype/phenotype data.
PLINK is used for a wide range of task, especially for quality controls, calculating relatedness between individuals and basic association analysis. 



## TASK: Association analysis and QC using PLINK
Today you will run a simple association analysis using PLINK. You should have installed PLINK on your computer. 
Because of time restraints we have made a small data-set that will run within reasonable time. 

The phenotype in today´s practical is dummy phenotype created for this practical but based on real data from the 1000-genomes project. The data-set consist of 1001 cases and 1001 controls and only chr 22.

#### TASK outline:
* identify and get to know your files
* run association analysis using PLINK
* plot and view your results
* run through the quality control of our data using PLINK
* re-run the association analysis using PLINK
* plot and view your results and evaluate differences

These are the files to use in this practical. We also need R/RStudio and a terminal window.
They contain the genotype data of individuals as well as the case control status:
  
* day2.bim
* day2.fam
* day2.bed

You can find useful information regarding these files here:
  [https://www.cog‐genomics.org/plink/1.9/formats#bed](https://www.cog‐genomics.org/plink/1.9/formats#bed)
   
*You can have a look at them in the terminal:*
*Open Ubuntu on your PC*
NB: This exercise will use Ubuntu as we are on Windows computers but will run with any terminal handling Linux commands.

We will first introduce some basic linux commands that will be nice to know (cd, less, head, tail, more) to navigate in a terminal.
    
###### In terminal:
You need to first make the terminal know where your files are/what directory to work from
and you can copy this and paste it into your terminal.  

For example:

```
cd /mnt/c/Users/User/Desktop/PPU-GenEpi-main/Day2/
```

We also need to tell the terminal where plink is:

```
plink=/mnt/c/Users/User/Desktop/PPU-GenEpi-main/PLINK/plink.exe
```

Then have a look at the top of the files:

```
less day2.bim | head
```
To exit `less` press `ESC`
```
less day2.fam | head
```

```
less day2.bed | head
```

This will display the the first ten lines in the terminal. If you use ```tail``` instead of ```head```,
you will the bottom of the file. By using ```more``` instead of ```head``` you can scroll down the file by 
hitting the space-button. If you wish to get out of the file, you press "q". 
These are examples of commands that are very useful when writing in a terminal with unix language.


When/If you try ```less day2.bed | head```, 
you will get the message 
> "day2.bed" may be a binary file. See it anyway?
- Just press "n" as it will make no sense for human eyes ;) 

*Then try:*
  
###### In terminal:
```
wc -l day2.fam
```

This will tell you how many lines there are in the day2.fam file. Try it on the .bim file as well

###### In terminal:
```
head -n1 day2.fam | wc -w
```

Will tell how many columns the day2.fam file has. Try it on the .bim file too. 


## Instructions/Scripts for running the analysis

#### PLINK:
The following command performs **association analysis in PLINK** using logistic regression:
  
###### In terminal:
```
$plink --bfile day2 --allow-no-sex --logistic beta --ci 0.95 --out output/results
```

This command uses the files with the name "day2*" as input, runs the analysis, and stores the output files in "output":
  
* results.assoc.logistic
* results.log
* results.nosex

Take a look at these files by using the commands enlisted above. The file "results.assoc.logistic" contains one line for each SNP. The entries
in each line are described by names in the header line of the file. The file contains the following columns. 

* **CHR**: Chromosome
* **SNP**: SNP identifier
* **BP**: Physical position (base-pair)
* **A1**: Tested allele (minor allele by default)
* **TEST**: Code for the test (see below)
* **NMISS**: Number of non-missing individuals included in analysis
* **BETA**: Regression coefficient
* **SE**: Standard error
* **L95**: Lower bound on confidence interval
* **U95**: Upper bound on confidence interval
* **STAT**: Coefficient t-statistic
* **P**: Asymptotic p-value for t-statistic

The other output files:
**results.log** shows the details regarding the analysis. Take time to read it as these
are the files which give you important information regarding the analysis run.
**results.nosex**: Take a look at this file. Are you able to find out what it is?

Every subsequent PLINK command will make these two other files together with the analyses output.

**Start an R session by writing "R" in the command line or open R/RStudio**

*The functions for R are also in the script R_script_Practical_Day2.R*

*We will set a working directory first in R as well to tell R where to find or save files.*  

###### An example in an R session:
```  
setwd("C://Users/User/Desktop/PPU-GenEpi-main/Day2/")
```

*The results*:, **results.assoc.logistic**, *can be read into R using the following command:*

###### In R session:
```  
assoc <- read.table("output/results.assoc.logistic", header=T, as.is=T)
```

Where the function **read.table** is used. Try **?read.table** to get information on the arguments used. When an argument is set to T, it means it is TRUE and is active. F denotes FALSE, and R ignores this argument. See the table above for more useful functions in R. 

The p-values (on –log10 scale) can be plotted across the chromosome:
  
**Note:** 
  * You need to have the packages "ggplot2" and "scales" installed

To install the R packages run this script:
```
install.packages(c('ggplot2', 'scales', 'AER', 'devtools', 'dlyr', 
  'data.table', 'stringr', 'optparse'))
```
```
library(devtools)
```
```
install_github("MRCIEU/TwoSampleMR")
```

* If you are not using Rstudio, you can view the plots in Finder or Outlook.
  
###### In R session:
```
library("ggplot2", "scales")
reverselog_trans <- function(base = exp(1)) {
  trans <- function(x) -log(x, base)
  inv <- function(x) base^(-x)
  scales::trans_new(paste0("reverselog-", format(base)), trans, inv,
            scales::log_breaks(base = base),
            domain = c(1e-100, Inf))
}

assoc_plot <- ggplot(assoc, aes(BP, P)) +
  geom_point() +  scale_y_continuous(trans=reverselog_trans(base=10),
                                     labels=scales::trans_format("log10", scales::math_format(10^.x))) +
  labs(title="-log10 p-values", x ="position", y = "-log10 p-value") +
  theme_bw()
assoc_plot

ggsave("output/assoc_plot.pdf", assoc_plot)
```


Now you have run association analysis by using PLINK on a dataset which have not been cleaned (QC´ed). 
Let us clean the dataset and see how it affects the association analysis. 

## Quality Control (QC)

We will continue to use the PLINK program to clean up the dataset. 
There are several different QC procedures that can be applied to a dataset. 
We will use PLINK to apply sample (person)-level filters on missing genotype 
call rate and heterozygosity, and to identify duplicate samples. 
Then.we will apply SNP-level filters on missing genotype call rate, allele 
frequency and HWE, and different genotype call rates between cases and controls. 
We will also examine quantile-quantile (Q-Q) plots of the p-values. 

Exit R session and return to the Terminal by writing **q()**.

### Sample Filters

**(a)  elevated missing data rates or outlying heterozygosity rate**
  
The first step involves calculating the sample statistics (missing genotype call rate and heterozygosity) of each sample. The following command calculates the 
sample statistic for missing rate:
  
###### In terminal:
```
$plink --bfile day2 --missing --allow-no-sex --out output/mis 
```
The sample statistics are written to

* mis.lmiss
* mis.imiss

The *fourth* column [N_MISS], in the file mis.imiss gives the number of missing SNPs  and the sixth column [F_MISS], gives the proportion of missing SNPs per individual.

*The following command **calculates the sample statistics for heterozygosity rate**:*
  
###### In terminal:
  ```  
$plink --bfile day2 --het --allow-no-sex --out output/het
```

The sample statistics are written to:
  
* het.het

The *third* column in the file **het.het** gives the observed number of 
homozygous genotypes [O(Hom)] and the *fifth* column gives the number of 
non-missing genotypes [N(NM)], per individual. 

The results can be visualized from within an R session by plotting a graph 
where the observed heterozygosity rate per individual is plotted on the x-axis 
and the proportion of missing SNPs per individuals is plotted on the y-axis:
  
### Evaluating QC results - missingness and heterozygosity - from Plink in R 
 
Reading the files:
  
###### In R session:
```
het <- read.table("output/het.het", header=T)

mis <- read.table("output/mis.imiss", header=T)
```


First, calculate the observed heterozygosity rate per individual using the formula (N(NM) - O(Hom))/N(NM)

Then, plot the proportion of missing genotypes and the heterozygosity rate

###### In R session:
```
mishet <- data.frame(FID=het$FID, IID=het$IID, het.rate=(het$N.NM - het$O.HOM)/het$N.NM, mis.rate=mis$F_MISS)

mishet_plot <- ggplot(mishet, aes(mis.rate, het.rate)) +
  geom_point() +
  labs(title="", x ="Proportion of missing genotypes", y = "Heterozygosity rate") +
  theme_bw()

mishet_plot

ggsave("output/mishet.pdf", mishet_plot)
```
**Note:** In this function we assume they are ordered in the same manner. 
This should always been checked before merging data from two or more 
data frames.

As you can see, there are samples with high levels of missing data and 
some samples with unusually high and low heterozygosity.
Based on this plot you can decide upon thresholds for removing 
individuals. For example, you might decide to remove all individuals 
with more than 1% missing data and a heterozygosity outside the range 
[0.23, 0.4].
Use the following R commands to visualize the thresholds and to save 
the R plot into a file:
  
### Evaluating which individuals to exclude based on set cut off values
  
###### In R session:
```
mishet_excl_plot <- ggplot(mishet, aes(mis.rate, het.rate)) +
  geom_point() +
  labs(title="", x ="Proportion of missing genotypes", y = "Heterozygosity rate") +
  theme_bw() +
  geom_vline(xintercept=0.01, linetype="dashed") +
  geom_hline(yintercept=0.23, linetype="dashed") +
  geom_hline(yintercept=0.4, linetype="dashed")

mishet_excl_plot

ggsave("output/mishet_excl_1.pdf", mishet_excl_plot)

```

*How many samples failed each of the filters? 
 - What is the total number of failures?*
 - These numbers often need to be reported in papers.

To answer these questions, you can extract the individuals you decided 
to exclude using the following R commands:
  
###### In R session:
```
#individuals with mis.rate > 0.01
fail_mis_qc <- mishet[mishet$mis.rate > 0.01,]

#individuals with het.rate <0.23 and individuals with  het.rate >0.4
fail_het_qc <- mishet[mishet$het.rate < 0.23 | mishet$het.rate > 0.4,]

```

You can then create files **fail_mis_qc.txt** and **fail_het_qc.txt** with the 
individuals that do not pass the thresholds:
  
##### In R session:
```{r }
write.table(fail_mis_qc, "output/fail_mis_qc.txt", row.names=F, col.names=T, quote=F)

write.table(fail_het_qc, "output/fail_het_qc.txt", row.names=F, col.names=T, quote=F)

```

*Take a look at these files to understand the format.*
  
  **(b)  duplicated or related individuals**
  
  Duplicated or related individuals can be identified using the 
identity-by-descent (IBD) report constructed from a reduced subset of 
frequent SNPs.

*First, create a subset of frequent snps (MAF > 0.35) while filtering
out SNPs that have > 5% missing rate and SNPs that don’t pass the HWE 
test at p-value =0.00000001:*
  
###### In Terminal:
 ```
$plink --bfile day2 --maf 0.35 --geno 0.05 --hwe 0.00000001 --allow-no-sex --make-bed --out output/frequent
```


*How many SNPs are we left with?*
```
head output/frequent.bim |head

wc -l output/frequent.bim
```

Reduce the subset of frequent SNPs by pruning so that no pair of SNPs
(within 50 base-pairs) has an r2 greater that a 0.2:
  
###### In Terminal:
```
$plink --bfile output/frequent --indep-pairwise 50 5 0.2 --allow-no-sex --out output/prunedsnplist

head output/prunedsnplist.prune.in | head
```

Generate the identity-by-descent (IBD) report from the reduced subset 
of frequent SNPs:
  
##### In Terminal:
```
$plink --bfile output/frequent --extract output/prunedsnplist.prune.in --genome --allow-no-sex --out  output/pruned

head output/pruned.genome |head
wc -l output/pruned.genome
```


Take a look at the file pruned.genome. 
* The columns [Z0], [Z1] and[Z2] contain the probability of sharing 0, 1 or 2 alleles IBD.  
* The column [PI_HAT] is a mean IBD per individual, i.e. (0\*genome$Z0 + 1\*genome$Z1 + 2\*genome$Z2)/2
This file is usually very large as it contains results for 
every pair of individuals. When using real data you may prefer to 
use the --Z-genome command instead, which outputs a gzipped version.

In R, you can visualize the IBD report and identify duplicated 
individuals by plotting the standard error of the IBD sharing vs.
the mean IBD sharing:
  
##### In R session:
```
#read in day2pruned.genome file
genome <- read.table("output/pruned.genome", header=T, as.is=T)

```
The pruned.genome file is big as it contains all possible pairs of 
individuals (2,003,001 in our case). Therefore, we look only at the 
pairs of individuals with PI_HAT > 0.1875, which corresponds to a half 
way between second and third degree relatives.

###### In R session:
```
genome <- genome[genome$PI_HAT > 0.1875,]

#compute Mean(IBD)
mean.ibd <- 0*genome$Z0 + 1*genome$Z1 + 2*genome$Z2

#compute Var(IBD)
var.ibd <- ((0 -mean.ibd)^2)*genome$Z0 +
  ((1 -mean.ibd)^2)*genome$Z1 +
  ((2 -mean.ibd)^2)*genome$Z2

#compute SE(IBD)
se.ibd <- sqrt(var.ibd)

ibd.stat <- data.frame(mean.ibd, se.ibd)

ibd.stat_plot <- ggplot(ibd.stat, aes(mean.ibd, se.ibd)) +
  geom_point() +
  labs(title="", x ="Mean IBD", y = "SE IBD") +
  theme_bw()

ibd.stat_plot

ggsave("output/ibd.stat.pdf", ibd.stat_plot)
```


People who on average share two alleles IBD are either monozygotic 
twins or duplicates. You can see two pairs of identical individuals 
(genome$Z2=1, i.e. they share all their alleles IBD).

###### In R session:
```
duplicate <- genome[mean.ibd == 2,]
duplicate
```

Save the ID of one in each pair into a file **fail_ibd_qc.txt** for subsequent removal:
  
###### In R session:
```{r }
fail_ibd_qc <- data.frame(FID=duplicate$FID2, IID=duplicate$IID2)

write.table(fail_ibd_qc, "output/fail_ibd_qc.txt", row.names=F, col.names=T, quote=F)
```

**NOTE:** In practice, you might want to look at the individual call 
rates stored in mis.imiss and output the IDs of the individual with 
the lowest call rate to **fail_ibd_qc.txt** for subsequent removal.

### Remove all individuals failing QC

In R, concatenate all the files listing individuals failing the 
previous QC steps into a single file **fail_qc.txt**:
  
###### In R session:
```{r }
fail_mis_qc <- read.table("output/fail_mis_qc.txt",header=T,as.is=T)

fail_het_qc <- read.table("output/fail_het_qc.txt", header=T, as.is=T)

fail_ibd_qc <- read.table("output/fail_ibd_qc.txt", header=T, as.is=T)

fail_qc <- data.frame(FID=c(fail_mis_qc$FID, fail_het_qc$FID,  fail_ibd_qc$FID), IID=c(fail_mis_qc$IID, fail_het_qc$IID,  fail_ibd_qc$IID))

fail_qc <- unique(fail_qc)

write.table(fail_qc, "output/fail_qc.txt", row.names=F, col.names=F, quote=F)
```

The file **fail_qc.txt** should now contain a list of unique individuals 
failing the previous QC steps.

Return to the terminal window and use plink to remove these from the 
dataset:
  
###### In Terminal:
```
$plink --bfile day2 --remove output/fail_qc.txt --allow-no-sex --make-bed --out output/qc.ind
```

**Note:** make sure to look at the screen output when these commands run.
Check the number of the removed individuals.

This command creates new files 
* qc.ind.bim
* qc.ind.fam
* qc.ind.bed

The new **ind.qc** dataset should contain 1939 individuals that passed 
the QC.

###  SNP Filters

In a similar way to the sample filters above, the first step involves 
calculating the SNP statistics (missing rate, allele frequency, 
                                p-value for the test of HWE and for the different call rates 
                                between cases and controls) for each SNP.


#### Excessive missing data rate

The following PLINK command computes the missing data rate:
  
###### In Terminal:
```
$plink --bfile output/qc.ind --missing --allow-no-sex --out output/qc.ind.mis

head output/qc.ind.mis.lmiss
```
The column **F_MISS** in the file **qc.ind.mis.lmiss** contains the 
fraction of missing genotypes per SNP. It can be visualised in R by 
plotting a histogram (here we choose to plot it on a log scale):
  
###### In R session:
```
#read in day2.qc.ind.mis.lmiss file
lmis <- read.table("output/qc.ind.mis.lmiss",header=T)

#plot histogram of the fraction of missing genotypes:
mis1_genot_plot <- ggplot(lmis, aes(F_MISS)) +
  geom_histogram(bins=9) +
  scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  labs(title="Fraction of missing data", x ="Fraction of missing genotypes", y = "Number of SNPs") +
  theme_bw()

mis1_genot_plot
ggsave("output/qc_ind_lmis_1.pdf", mis1_genot_plot)
```

Based on this plot you can decide upon thresholds for removing SNPs, 
e.g. you might want to remove SNPs with more than 4% missing data:
  
###### In R session:
```
mis_genot_plot_2 <- ggplot(lmis, aes(F_MISS)) +
  geom_histogram(bins=9) +
  scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  labs(title="Fraction of missing data", x ="Fraction of missing genotypes", y = "Number of SNPs") +
  theme_bw() +
  geom_vline(xintercept=0.04, linetype="dashed")

mis_genot_plot_2

ggsave("output/qc_ind_lmis_2.pdf", mis_genot_plot_2)
```

Next, examine the frequency of the minor alleles by computing them in 
PLINK and visualizing in R:
  
###### In terminal:
```
$plink --bfile output/qc.ind --freq  --allow-no-sex --out output/qc.ind.freq
```

###### In R session:
```
#read in frequencies file
freq <- read.table("output/qc.ind.freq.frq", header=T)

# plot minor alleles
maf_plot <- ggplot(freq, aes(MAF)) +
  geom_histogram(bins=10) +
  labs(title="Minor allele frequencies", x ="MAF", y = "Number of SNPs") +
  theme_bw() +
  geom_vline(xintercept=0.01, linetype="dashed")

maf_plot

ggsave("output/maf.pdf", maf_plot)
```

#### Different genotype call rates between cases and controls

The following PLINK command performs a test for differences between 
cases/controls missing call rates at each variant and outputs the 
results into a file **qc.ind.call.rate.missing**.

###### In Terminal:
```
$plink --bfile output/qc.ind --test-missing --allow-no-sex  --out output/qc.ind.call.rate

head output/qc.ind.call.rate.missing |head
wc -l output/qc.ind.call.rate.missing
```

Next, you create a file **fail_diffmiss_qc.txt** that contains SNPs with 
different call rates between cases and controls at p-value < 0.000001 
for subsequent removal:
  
###### In R session:
```
#read in qc.ind.call.rate.missing file
diffmiss <- read.table("output/qc.ind.call.rate.missing", header=T, as.is=T)

#save SNPs with p-value < 0.000001 into a file
diff.miss <- diffmiss[diffmiss$P<0.000001,]

head(diff.miss); dim(diff.miss)
write.table(diff.miss$SNP, "output/fail_diffmiss_qc.txt", row.names=F, col.names=F, quote=F)

```

#### Remove all SNPs failing QC

The following PLINK command removes the SNPs listed in the file 
**fail_difmiss_qc.txt**, and also SNPs with MAF <0.01, 
SNPs with >4% missing genotypes, and SNPs that depart from the HWE 
at p-value <0.000001 (in controls):
  
###### In Terminal:
```
$plink --bfile output/qc.ind  --exclude output/fail_diffmiss_qc.txt  --maf 0.01 --hwe 1e-6  --geno 0.04  --make-bed --allow-no-sex  --out output/practical2.qc
```


The new dataset **day2.qc** contains 5016 SNPs that passed the QC.
1939 persons, 969 cases, 970 controls

## Testing for association after QC

Now the filtering has been carried out you can re-test for association 
to see if it has made any difference to the results using PLINK:
  
###### In Terminal:
```
$plink --bfile output/practical2.qc --allow-no-sex --logistic beta --ci 0.95 --out  output/day2.qc.log
```

###### In R session:
```
#read in the association test results (before QC)
assoc <- read.table("output/results.assoc.logistic", header=T, as.is=T)

#read in the association test results (after QC)
assoc.qc <- read.table("output/day2.qc.log.assoc.logistic", header=T, as.is=T)


QC_before_plot <- ggplot(assoc, aes(BP, P)) +
  geom_point() +  scale_y_continuous(trans=reverselog_trans(base=10),
                                     labels=scales::trans_format("log10", scales::math_format(10^.x))) +
  labs(title="-log10 p-values (before QC)", x ="position", y = "-log10 p-value") +
  theme_bw()
QC_before_plot

QC_after_plot <- ggplot(assoc.qc, aes(BP, P)) +
  geom_point() +  scale_y_continuous(trans=reverselog_trans(base=10),
                                     labels=scales::trans_format("log10", scales::math_format(10^.x))) +
  labs(title="-log10 p-values (after QC)", x ="position", y = "-log10 p-value") +
  theme_bw()
QC_after_plot

# Saving plots
ggsave("output/comparison_beforeQC.pdf", QC_before_plot)
ggsave("output/comparison_afterQC.pdf", QC_after_plot)
```

It seems as if most of the highly significant results have disappeared 
after QC. There is one (possibly slightly suspicious?) signal with a 
p-value around 1e-14 and a few better-supported signals with p-values 
around 1e-06 or 1e-07.


#### Q-Q plots of p-values

You can examine the quantile-quantile (Q-Q) plot of the p-values before 
and after the QC by plotting the observed p-values vs. the expected 
p-values:
  
  
###### In R session:
  
```
#observed –log10 p-values (before QC)
p.obs <- -log10(sort(assoc$P,decreasing=F))

#expected –log10 p-values (before QC)
p.exp <- -log10( 1:length(p.obs)/length(p.obs) )

qq.stat <- data.frame(p.obs, p.exp)

qq_plot_beforeQC <- ggplot(qq.stat, aes(p.exp, p.obs)) +
  geom_point() +
  labs(title="QQ plot (before QC)", x ="expected -log10 p-value", y = "observed -log10 p-value") +
  theme_bw() +
  geom_abline(intercept = 0, slope = 1, color="grey", linetype="dashed")

qq_plot_beforeQC

#save plot
ggsave("output/qq_plot_beforeQC.pdf", qq_plot_beforeQC)

#observed –log10 p-values (after QC)
p.obs <- -log10(sort(assoc.qc$P,decreasing=F))

#expected –log10 p-values (after QC)
p.exp <- -log10( 1:length(p.obs)/length(p.obs) )

qq.stat <- data.frame(p.obs, p.exp)

qq_plot_afterQC <- ggplot(qq.stat, aes(p.exp, p.obs)) +
  geom_point() +
  labs(title="QQ plot (after QC)", x ="expected -log10 p-value", y = "observed -log10 p-value") +
  theme_bw() +
  geom_abline(intercept = 0, slope = 1, color="grey", linetype="dashed")

qq_plot_afterQC

#save plot
ggsave("output/qq_plot_after.pdf", qq_plot_afterQC)

```

#### Identifying regions of association

Now the dataset has been cleaned you can identify a set of SNPs for 
follow-up replication. To do this you will need to pick a threshold 
for the p-values. For example, SNPs with a -log10 (p-value) > 4 can be 
identified by typing in R:
  
###### In R session:

```
assoc.qc[-log10(assoc.qc$P)>4,]
```


###### Final note:

Remember to fill the following form during the exercise today: 
  [https://nettskjema.no/a/284234](https://nettskjema.no/a/284234)

Thank you!
