# Pincho
A Modular Approach to High Quality de novo Transcriptomics

# Table of Contents
(1) Introduction<br />
(2) What's Inside?<br />
(3) Installation<br />
(4) Running Pincho<br />
(5) Troubleshooting<br />
(6) Subprograms

# Introduction
Pincho is:<br />
- a modular methods workflow for transcriptomics.<br />
- designed for short read cleaning, assembly, validation, annotation and expression analysis.<br />
- currently designed for Ubuntu 20.04 Focal Fossa LTS for workstations with at least 12 threads, 64GB ram.
  
# What's Inside?
Pincho Master folder consists of:<br />
- bin folder with all subprograms<br />
- lib folder with some premade annotation databases (kegg, trembl_amphibians, uniprot)<br />
- COPYING.txt our license agreement<br />
- png files for our pincho icon and UI<br />
- pincho.conf or a printed methods list of the last run<br />
- pincho_conf.py the config script used to launch pincho UI<br />
- pincho.py the main workflow of pincho<br />
- P_installer.py an automated installer script to be launched in its residing location
  
# Installation
- [1] Download Pincho_Master from this github page<br />
- [2] Open Terminal<br />
- [3] $cd Pincho_Master<br />
- [4] $P_installer.py<br />
- [5] be sure to answer 'y' when prompted by the automated installer<br />
- [6] installation complete!<br />
<br />
[Alt] as an alternative to launching the installer you can open the installer file as a text file and manually run the commands in terminal

# Running Pincho
Pincho is designed to iterate through a large quantity of folders in a given directory, thus increasing its performance as it can run for months without user input. Because of this mechanic we impose a strict regimen of directory structure prior to launching Pincho.<br />
<br />
[1] Working folders, or folders containing data to be processed by Pincho must conform to one of the following formats:<br />
- [a] folder containing solely raw reads as fastq.gz or fasta.gz. Pincho will interpret 1 file as single-end data, 2 files as paired-end data. If there are more than 2 fastq.gz or fasta.gz files then Pincho will NOT perform ideally and it will grab the first 2 fasta.gz/fastq.gz files detected.<br />
- [b] if you are not starting from raw NGS data and are illiciting a portion of Pincho then the above rule does not necessarily need to be followed unless the activity requested requires raw data.<br />
  <br />
[2] Working folders, and ONLY working folders, must be in the same parent directory. Pincho will process the folders in the parent directory, thus it is important that files that are not intended to be processed, not reside in that directory.<br />
<br />
Example of ideal directory structure for running Pincho:<br />
<br />
-----Desktop-----NGS_Project_1-----Species_A_Tissue_Type_A-----read_1.fq.gz, etc.<br />
__________________________________-----Species_A_Tissue_Type_B-----read_1.fq.gz, etc.<br />
_________________________________________________________________-----etc.<br />
<br />                              
Once your directory is in the proper format:<br />
- [1] Open terminal<br />
- [2] Change working directory to the parent folder of your working folders (i.e. NGS_Project_1)<br />
- [3] $Pincho_conf.py<br />
- [4_alt] If you wish to test Pincho, then press the 'Test' button on the UI and click 'Start.'<br />
- [4] Familiarize yourself with the parameters, then when ready press 'Clear' and create your pipeline.<br />
- [5] Press 'Start' to start your created pipeline

# Troubleshooting
The most common source of stress for the program will be the thread count and memory count. Please make sure to have at least 12 threads and 64GB of ram on the systema and please also scale in ratio (3:16). Having more threads than ram in this ratio will sometimes cause subprograms to fail.<br />
<br />
For any other problem please reach out to us directly.

# Subprograms
As Pincho is a methods workflow, subprograms used are stated and given credit below:<br />
- abyss
- augustus
- bbmap
- binpacker
- blat
- boost
- busco
- bwa
- cd-hit
- circos
- circos-tools
- diamond
- fastx_toolkit
- hisat2
- hmmer
- idba
- jellyfish
- jre
- kallisto
- map2mitos
- megahit
- mira
- mitobim
- mitoz
- ncbi-blast
- oases
- orp-transrate
- paladin
- python
- quorum
- rcorrector
- rsem
- salmon
- samtools
- seqkit
- shannon
- soapdenovo-trans
- spades
- sratoolkit
- transabyss
- transcriptomeAssemblyTools
- translig
- transrate
- trimmomatic
- trinity
- velvet

Packages:
- curl
- python
- python3
- pip
- pip3
- numpy
- gcc
- make
- cmake
- libbz2-dev
- libncurses5-dev
- zlib1g-dev
- libncursesw5-dev
- liblzma-dev
- libcurl4-openssl-dev
- libssl-dev
- libsparsehash-dev
- shmlast
- bowtie2
- samtools
- pandas
- biopython==1.77
- tqdm
- default-jre
- git
- dos2unix
- cd-hit
- pkgconf
- sh
- cvxopt
- metis
- r-base
- python-igraph==0.8.3
- csvtool
- libboost-iostreams-dev
- libgsl-dev
- libboost-graph-dev
- libboost-all-dev
- libsuitesparse-dev
- liblpsolve55-dev
- libsqlite3-dev
- libmysql++-dev
- libbamtools-dev
- python-tk
- python3-tk
- pigz

What is not currently in use is intended for use in future versions of Pincho.
