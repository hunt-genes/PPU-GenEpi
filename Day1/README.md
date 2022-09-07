# DAY 1 - Gitlab orientation and laptop setup for practicals

### FOCUS and LEARNING GOALS
> The aim for this session is for you to 1) become familar with gitlab and how to navigate the practicals, 2) 
> download and install the required programs, and 3) test the scripts used in the practicals.

### Required downloads   

R version 4.0.0 https://cran.uib.no/   
RStudio https://www.rstudio.com/products/rstudio/download/   

If you are using windows, please install bash to run the suggested commands in the terminal:
https://itsfoss.com/install-bash-on-windows/   

You can use Plink from the zip archive here: https://github.com/hunt-genes/SMED8020/tree/main/Day2 
or from the originating site for PLINK: (PLINK 1.9) https://www.cog-genomics.org/plink/ 

### Specific installations and tests for practicals

#### Practical 2 - Setup for PLINK QC: 

Download the zip archive in https://github.com/hunt-genes/SMED8020/tree/main/Day2

Store the files on your local computer and direct the terminal to them by writing ```cd filepath```
This is how you make the terminal know where your files are/what directory to work from
and you can copy this and paste it into your terminal.
For example:
```cd /mnt/c/Users/name00/Desktop/SMED8020/2021/Day2/```

The suggested structure would be a folder that has Plink, the three data day2 files and a subfolder called output.
Note that the Plink.exe file needs to have the same file path as the data

If you are new to linux commands, you could make yourself familiar with basic commands in sites like: https://www.hostinger.com/tutorials/linux-commands

###### In R session:
Libraries needed for R are ggplot2 and scales   

```
install.packages("ggplot2", "scales")
```
#### Practical 3 - GWAS Meta-analysis

**Option 1:** Download the precompiled binary for METAL that matches your operating system from [University of Michigan Center for Statistical Genetics](http://csg.sph.umich.edu/abecasis/metal/download/).
Note: This doesn't seem to work on newer Mac OS. 

Store the exe file on your local computer. 

Mac: The .tar.gz will likely be automatically unzipped on a Mac into a folder called `generic-metal`.   
PC: You may need to use WinZip to unzip the file on a PC.  
Linux: Via the linux command line you can unzip with `tar -xf <path/to/.tar.gz`.  


**Option 2:** Download and build a new version of [METAL from source code hosted on GitHub](https://github.com/statgen/METAL)

You can download the GitHub repo as a zip file or use `git clone https://github.com/statgen/METAL.git` (this requires `git` to be installed).

You also need to install of [CMake](https://cmake.org/install/) with downloads of the pre-compiled binaries for your OS [here](https://cmake.org/download/)

```
mkdir build && cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
make
make test
make install
```

Depending on your working directory will use the path to `generic-metal/metal` as the command on Day 3. We will download the necessary data on Day 3.
You can make sure it works by running `path/to/metal --help`

Download the data you wil use [here](https://ntnu.box.com/s/rvytm8ymd8iple8negy8ix8x5vp7qs9a). You will need about 1.7 GB free on your machine.

#### Practical 5 - Mendelian randomization: https://github.com/hunt-genes/SMED8020/tree/main/SMED8020-2022-master

Test your setup   
-Download the files for practical 5 from github: https://github.com/hunt-genes/SMED8020/tree/main/SMED8020-2022-master 
-Open Rstudio or R   
-Open  Mendelian_Randomization_Practical_1_-_Single_Sample_-_Rscript.R   
-Edit input directories   
-Run   
-Open  Mendelian_Randomization_Practical_2_-_MR_base_-_Rscript.R   
-Edit input directories   
-Run   
