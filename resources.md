# resources
List of resources to share for learning

# Learning R
* [swirl](https://swirlstats.com/): R package that turns your R consule into a learning environment, for beginners
* [Riffomonas](https://riffomonas.org/): University of Michigan professor makes tutorials and YouTube videos 
* [CodeAcademy](https://www.codecademy.com/learn/learn-r): Modules for R *note this is paid*. 
* [Cookbook for R](http://www.cookbook-r.com/): Solutions to common tasks with code examples
* [DataCamp](https://www.datacamp.com/): Need to register but introductory courses and chapters of modules are free
* [R Tutorial by RBloggers](https://www.r-bloggers.com/2015/12/how-to-learn-r-2/): Learning path that brings together many resources
* [More R resources for free](https://www.r-bloggers.com/2016/02/learning-r-for-free-free-online-resources/)

# Learning Python
* [Girls Who Code at UM DCMB](https://github.com/GWC-DCMB/GWC-DCMB/blob/master/get-started.md): Python notebooks you can use to learn data science for Python
* [Google's Python class](https://developers.google.com/edu/python): For people with a little programming experience who want to learn Python; has videos and code exercises
* [Python documentation](https://docs.python.org/3/)
* [LearnPython](https://www.learnpython.org/): Free interactive Python tutorial

# Stats
* [StatQuest](https://www.youtube.com/channel/UCtYLUTtgS3k1Fg4y5tAhLbw)
* [Linear-regression](https://mlu-explain.github.io/linear-regression)

# Genomics Tools 
* [gnomAD browser](https://gnomad.broadinstitute.org/): Aggregating exome and genome sequencing data
* [PLINK](https://www.cog-genomics.org/plink/): Tool for storing and analyzing genomic data
* [UCSC Genome Browser](https://genome.ucsc.edu/): Integrate datasets to understand the genome
* [OpenTargets](https://www.opentargets.org/): Using human genetics for systematic drug discovery
* [GWAS Catalog](https://www.ebi.ac.uk/gwas/): Browse published genome wide association summary statistics
* [PGS Catalog](https://www.pgscatalog.org/): Download markers and weights from published polygenic risk scores
* [HaploReg](https://pubs.broadinstitute.org/mammals/haploreg/haploreg_v4.php): Tool for exploring annotations of the noncoding genome at variants on haplotype blocks
* [TOPMed BRAVO](https://bravo.sph.umich.edu/freeze8/hg38/): Browse genetic variants from TOPMed

# Biomedical research tools
* [BioConda](https://bioconda.github.io/): Python software packages for biomedical research 
* [Type 2 diabetes Knowledge Portal](https://t2d.hugeamp.org/): Data and tools for T2D research
* [Cardiovascular disease knowledge portal](https://cvd.hugeamp.org/): Data and tools for CVD research

# UNIX Cheat Sheet
* `ssh user@server`: Start SSH sessiom (exit with exit)
* `ls`: List files in directory
* `ls -lah <directory>`: Show content of directory for humans, including hidden files
* `cd <directory>`: Change directory
* `mkdir <directory>`: Create directory
* `pwd`: Print working directory (the directory you are in)
* `rm -rf <directory>`: Delete full directory
* `rm file`: Delete file
*  `wc -l`: Number of lines in a file
*  `zless <file.gz> | wc -l`: Count lines in a gzipped file, can also use zcat and possibly less (if computer knows to use less with zipped files)
*  `less <file>`: Display the file; to exit press q, to move down, use space bar
*  `head <file>`: Display the first 10 line of the file to stdout
*  `mv <file1> <file2>`: Move or rename a file
*  `cp <file1> <file2>`: Copy a file to another location
*  `cd ../`: Change directory up to parent directory
*  `cat <file>`: prints the contents of files and can be use to concatenate multiple files together
*  `gzip`: Zip a file with gzip compression format
*  `gunzip -c`: Write to std out and keep original files when you decompress a .gz file
*  `unzip`: Decompress  a .zip file
*  `apt install <package>`: Install a package 
*  `echo "your text here"`: displaying lines of text passed 
*  `<command> | <command>`: The pipe signal passes the output from command 1 into command 2
*  `grep "string" <file> `: Searches for a string in a file
*  `sort <file>`: Sorts a file alphabetically or numerically, can use -k to specify specific columns
*  `cut -f 1 <file>`: Cut the first column from a file
*  `sort <file> | uniq`: Find the unique elements in a list, need to sort first
*  `sort <file> | uniq -c`: Count how many times you see each unique element
* Control-c : Sends interrupt signal, can be used to get to a new line
* Control-k: Delete everything on the line
* Control-a: Jump to the front of the line
* Control-e: Jump to the end of the line
* `*`: wild card; you can use like `ls *.gz` to list anything with the `.gz`. suffix 
* `>`: write to a file; you can use like `echo "my text" > my_file` 
* `>>`: append to a file; you an use like `echo "my text" > my_file` 
* If you hit the arrow key upwards you can access previous commands
* [Unix tutorial](https://www.coursera.org/learn/unix) on Coursera

# Awk
Awk is a language you can use on the command line. [User's Guide](https://www.gnu.org/software/gawk/manual/gawk.html)
* `awk 'NR==1 {print}'` Print the first row (NR==1)
* `awk 'NR==1 || ($1 > 1) {print}'` Print the first row OR `||` lines with the a first column greater than 1
* `awk '($1 > 1) && ($2 > 3) {print}'` Print rows with the first column (`$1`) greate than 1 AND `&&` rows with the second column (`$2`) greater than 3
* `awk '!($1 > 1) {print}'` Print rows with the first column NOT (`!`) greater than 1

# Command line OS
* [Unix vs Linux](https://www.geeksforgeeks.org/linux-vs-unix/): Unix is a operating system, Linux is a Unix-like operating system that is free and opensource 
* Ubuntu is an operating system based on the Debian Linux distribution and is free and open source 
* [sh vs bash](https://www.geeksforgeeks.org/difference-between-sh-and-bash/): bash is a shell of the Unix operating system, like sh but with more features, sh is a shell of the Unix operating system 
* Drag and drop a folder into the terminal to have the path show up
