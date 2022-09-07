# DAY 3 - Meta-analysis for quantitative traits using METAL

## FOCUS and LEARNING GOALS

>  The aim for this session is to get familiar with running a genome wide association study meta-analysis which  enables researchers to  gather data from many studies and analyse them together. 

There are several motivations for meta-analysis. One is the ability to increase power to detect small effect sizes or rare variant effects by increasing the study sample size. Many methods for meta analysis rely on using summary statistics therefore rendering the need to share individual level data unnecessary. This makes it easier to share data for meta analysis as the summary statistics are not deemed sensitive information and are typically made publicly available when papers are published in peer-reviewed journals. Finally, meta-analysis across genetic ancestries is the most statistically robust approach rather than pooling all ancestries together in one GWAS as it results in little or no loss of efficiency (as compared to analysis of combined data-sets) and reduces population stratification.
Some [slides](SMED8020_2022_MetaAnalysis.pdf) with extra informaiton may be helpful.  

**Suggested reading:**

* [Willer, C. J., Li, Y. & Abecasis, G. R. METAL: fast and efficient meta-analysis of genomewide association scans. Bioinformatics 26, 2190–2191 (2010).](https://academic.oup.com/bioinformatics/article/26/17/2190/198154)  
* [Nielsen, J. et al. Loss-of-function genomic variants with impact on liver-related blood traits highlight potential therapeutic targets for cardiovascular disease. Biorxiv. (2019) ](https://www.biorxiv.org/content/10.1101/597377v1)  
* [Evangelos Evangelou & John P. A. Ioannidis. Meta-analysis methods for genome-wide association studies and beyond. Nature Reviews Genetics. 2013](https://www.nature.com/articles/nrg3472)

[METAL](http://csg.sph.umich.edu/abecasis/metal/) was developed at the University of Michigan as a tool for meta-analysis of  genome-wide association analysis. For running METAL you can use either test statistics and standard error or p-values. For more info on METAL see the web-links in this document and the suggested readings paper. 

## TASK: Running a trans-ethnic meta-analysis using METAL  

Today you will run a meta-analysis to combine three studies using METAL. METAL has been pre-installed on our lab. Because of time restraints we have made a small data-set that will run within reasonable time. The data will not generate significant results. A separate set of files will therefore be used for plotting results. 
    
The phenotype in today's practical is low density lipoprotein (LDL) cholesterol and we will be using data from three large studies: HUNT, Biobank Japan and Global Lipids Genetics Consortium (GLGC).  

### Task outline  

0. Install METAL
1. Gather summary statistics from GWAS for low density lipoprotein (LDL) cholesterol in three separate studies (HUNT, Global Lipids Genetics Consortium, and UK Biobank) and check details for the files.
2. Run a meta-analysis using METAL
3. Have a coffee or a biobreak or ask questions to the lecturers. Questions for consideration are in ****bold****. 
4. View meta-analysis results

### Important points to consider  

When running a meta analysis there are many issues that need to be addressed.
* the availability of summary statistics
* phenotype related questions such as:
  * what phenotypes are available?
  * what are the phenotypes based on (self-reported, Electronic Health Records, physician curated)?
  * how are the phenotypes constructed?
  * are they comparable to your defined phenotype?
* genotype related questions such as: 
  * differences in genotyping and imputation, some markers will be study specific
  * genome build
  * flipped markers
  * population stratification

## Instructions  
#### 0. [Install METAL](http://csg.sph.umich.edu/abecasis/Metal/download/) using the pre-compiled binaries for your operating system. See [Practial Day 1](https://github.com/hunt-genes/SMED8020/tree/main/Day1) for details. 

#### 1. Organizing summary statistics  

Usually you would download publically available summary statistics from the internet to your local machine. For convience for this practical, the data can be downloaded from [here](https://ntnu.box.com/s/rvytm8ymd8iple8negy8ix8x5vp7qs9a). You will need about 1.7 GB.

* The original summary statistics from Biobank Japan (BBJ) of LDL cholesterol in N=72,866 can be found [here](https://humandbs.biosciencedbc.jp/files/hum0014/hum0014_README_QTL_GWAS.html)  
`BBJ-LDL-preMeta.txt`  
The columns are CHR     POS38   SNPID   Allele1 Allele2 AC_Allele2      AF_Allele2      N       BETA    SE      p.value log10P

* The original summary statistics of joint analysis of metabochip and GWAS data for LDL cholesterol in N=89,138 from the Global Lipids Genetics Consortium (GLGC) can be found [here](http://csg.sph.umich.edu/willer/public/lipids2013/)  
`GLGC-LDL-preMeta.txt`  
The columns are SNP_hg18        SNP_hg19        rsid    A1      A2      beta    se      N       P-value Freq.A1.1000G.EUR

* The summary statistics of LDL cholesterol from the HUNT study in N=67,429.   
`HUNT-LDL-preMeta.txt`  
The columns are CHR     POS38   SNPID   Allele1 Allele2 AC_Allele2      AF_Allele2      N       BETA    SE  p.value 

#### 2. Check your summary statistics to make sure they're ready for meta-analysis.

2.1 Which genome build is used?

The human reference genome has been updated over the years and variants are given different coordinates in different versions. 
The latest human reference genome GRCh38 was released from the Genome Reference Consortium on 17 December 2013.  
The previous human reference genome (GRCh37) was the nineteenth version (hg19).  
The version before this was NCBI Build 36.1	released March 2006	(hg18). 
You can see more [here](https://genome.ucsc.edu/FAQ/FAQreleases.html#release1). hg19 is still widely used and people are slowly converting to hg38.  
****From the summary statistic headers, can you tell what reference genome versions are used for each study?****  

It looks like BBJ and HUNT have SNP coordinates from hg38, but GLGC has summary statistics from hg18 and hg19. 
We must use [UCSC listOver](https://genome.ucsc.edu/cgi-bin/hgLiftOver) to convert the hg19 coordinates to hg38 before meta-analysis. We can use liftOver on the command line or via the web. To avoid extensive file manipulation on your part, we already used this .bed file to make a new version of the GLGC results: `GLGC-LDL-hg38-preMeta.txt`. This file also has a header that is consistent with the other two files. You will use this in the meta-analysis. The instructions for using liftOver are below.

Create a .bed file file from GLGC-LDL-preMeta.txt using linux tools `awk` and `sed` in the following the command:
`awk ' NR > 1 {print $2"\t"$3"\t"$4"\t"$5}' GLGC-LDL-preMeta.txt | sed 's/:/\t/g' | awk '{print $1"\t"$2-1"\t"$2"\t"$1":"$2"\t"$4}' > GLCG.hg19.bed`

**Web option:**
Upload the .bed file you made [here](http://genome.ucsc.edu/cgi-bin/hgLiftOver)

**Command line option:**
[Download liftOver](https://hgdownload.soe.ucsc.edu/admin/exe/)
You can use `wget` like so: `wget https://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/liftOver`
Turn on the executable bit `chmod +x ./filePath/utility_name`
Now `./filePath/utility_name` is executable.

[Download the map.chain](https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/) for hg19 to hg38  

The liftover command requires 4 parameters in this order: 
1) oldFile (in .bed format) 
2) map.chain 
3) newFile (just the name) 
4) unMapped
Execute this command:
`liftOver GLCG.hg19.bed hg19ToHg38.over.chain GLGC.h38.bed GLGC.hg38.unmapped`

Look in GLGC.hg38.unmapped. ****Were there some markers that did not get converted from hg19 to hg38? Why do you think that is?****  

The code to create the file with compatible header is here:
`join -1 4 -2 2 <(sort -k 4 GLGC.h38.bed) <(sort -k 2 GLGC-LDL-preMeta.txt) | awk -v OFS='\t' '{$5=toupper($5);$9=toupper($9)}1' | awk '{print $0"\t"substr($2, 4)"\t"$2":"$4":"$9":"$5}'  | sed  '1i\CHRPOS\tchr\tstart\tPOS38\tAllele2\tCHRPOS37\trsid\ta2\tAllele1\tBETA\tSE\tN\tp.value\tAF_Allele2\tCHR\tSNPID' > GLGC-LDL-hg38-preMeta.txt`

2.2 Check the file formats and headers

What is the header?
`head -n 1 file`

****Are your SNPIDs across the files formatted in the same way?****

2.3 How many variants will we be meta-analyzing?
****How many variants are in each of the files?****  
```
wc -l BBJ-LDL-preMeta.txt
wc -l HUNT-LDL-preMeta.txt
wc -l GLGC-LDL-hg38-preMeta.txt
```
The HUNT summary statistics originally had many variants because imputation was done with the TOPMed imputation panel, which allows for higher resolution imputation due to the large amount of sequencing samples which make up the reference panel. We have subsetted the HUNT file to include only variants seen in GLGC or BBJ. This makes the file a more manageable size, and we only will perform meta-analysis on variants tested in 2 or more studies.

****What imputation panel was used for GLGC?**** HINT: Check the methods of the [paper](https://www.nature.com/articles/ng.2797).

****How many genome wide significant results are in each of the input files?****  
```
awk '$11 < 5e-8 {print 0}' HUNT-LDL-preMeta.txt | wc -l
awk '$11 < 5e-8 {print 0}' BBJ-LDL-preMeta.txt | wc -l
awk '$10 < 5e-8 {print 0}' GLGC-LDL-hg38-preMeta.txt | wc -l
```

#### 3. Running METAL   

The [Wiki page for METAL](https://genome.sph.umich.edu/wiki/METAL_Documentation#Brief_Description)  may be useful.

Input files: 
* A text file for each study with results, summarized as a table. NB column separators must be specified. 
* A column with marker name, which should be consistent across studies 
* A column indicating the tested allele 
* A column indicating the other allele 
* If you are carrying out a sample size weighted analysis (based on p-values), you will also need: 
  * A column indicating the direction of effect for the tested allele 
  * A column indicating the corresponding p-value 
  * An optional column indicating the sample size (if the sample size varies by marker) 
* If you are carrying out a meta-analysis based on standard errors, you will need: 
  * A column indicating the estimated effect size for each marker 
  * A column indicating the standard error of this effect size estimate 
  * The header for each of these columns must be specified so that METAL knows how to interpret the data. 
 
A shell wrapper script will be used to create the config file needed to run METAL. `LDL_metal.sh` has been created for you. You can run it with the following commands:    
3.1. Create a config file with the bash script `LDL_METAL.sh` by filling in the appropriate arguments instead of "file1",  "file2",  "file3" and using "LDL_METAL" as your output prefix.
`bash LDL_METAL.sh  file1 file2 file3 LDL_METAL > LDL_METAL.conf`   

e.g. `bash LDL_metal.sh HUNT-LDL-preMeta.txt GLGC-LDL-hg38-preMeta.txt BBJ-LDL-preMeta.txt LDL_METAL_META > LDL_METAL.conf`    

3.2. Run metal with the config file (this should take less than 20 minutes)
`metal LDL_METAL.conf > LDL_METAL.log`  
Note: If you would like to time your analysis you can use the time program.  
`/usr/bin/time -o test_time -v metal LDL_METAL.conf`

****What type of meta-analysis did you run (fixed or random effects? sample size or inverse variance based?) What is the difference?****  
****Did you use genomic control? In what situations is it useful to use genomic control****  
****What does it mean to set the minimum weight to 10,000?****   
****What is the difference between "ANALYZE" and "ANALYZE HETEROGENEITY"?****  
****How might you create the config file if your summary statistics files had different header labels?****  

#### 4. View the meta-analysis results

Some informative outptut was printing to stdout as METAL was running. ****What was the smallest p-value and how many markers was the meta-analysis completed for?****  
There will be a .tbl and .tbl.info file created from the meta-analysis. You can use `less` to view the files.
****Will we use the same genome wide significance threshold for the meta-analysis as we used for the GWAS? Why or why not?****  
****How many genome wide significant results are there now?**** HINT: Use code like in #2 but replace `$10` with the column number with the p-value and use the file name for your meta-analysis results.

5. Note: We pre-processed the files so you don't have to subset the results to markers in >1 study, but you might need this information in the future if you have not pre-processed your input files.

METAL will perform a meta-analysis even on markers which are only present in one of the sub-studies. We are only interested in markers present in more than one study. 
The column labelled "direction" shows '?', '+', or '-' to indicate missingness, positive direction of effect, or negative direction of effect, respectively.  
One can use the `subset_meta_analysis.r` Rscript to exclude markers with more than one '?'.  
Execute the following command to subset the results. This will take < 5 minutes.  
`Rscript subset_meta_analysis.r --input LDL_METAL_META1.tbl --output LDL_METAL_MultiStudy.txt`

#### 6. Plot the meta-analysis results

To visually inspect your results for significant findings you can make a QQ-plot like Day 2's practical. ****How does the inflation appear to you?****  
`Rscript QQplot.r --input LDL_METAL_META1.tbl --pvalue P-value --af Freq1 --prefix LDL_METAL_MultiStudy --break.top 120`  
The file should exist in whatever the default directory your R is writing into. You can find this with `pwd`.  


****What is the lambda value for the smallest minor allele frequency (MAF) bin?****  
`cat *_lambda.txt`. 


