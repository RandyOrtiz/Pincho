# Pincho v0.1
A Modular Approach to High Quality de novo Transcriptomics

# Download
Available to download from our github releases.

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
  
# Installation - Normal
- [1] Download Pincho_Master from releases<br />
- [2] Open Terminal<br />
- [3] $cd Pincho_Master<br />
- [4] $P_installer.py<br />
- [5] be sure to answer 'y' when prompted by the automated installer<br />
- [6] installation complete!<br />
<br />
[Alt] as an alternative to launching the installer you can open the installer file as a text file and manually run the commands in terminal

# Installation - Light
- [1] Download Pincho_Master from releases<br />
- [2] Open Terminal<br />
- [3] $cd Pincho_Master/bin <br />
- [4] $P_application_fetcher.sh <br />
- [5] $cd Pincho_Master <br />
- [6] $P_installer.py<br />
- [7] be sure to answer 'y' when prompted by the automated installer<br />
- [8] installation complete!<br />
<br />
[Alt] as an alternative to launching the installer you can open the installer file as a text file and manually run the commands in terminal

# Running Pincho -Single Run
- [1] Open terminal<br />
- [3] $Pincho_conf.py<br />
- [4_alt] If you wish to test Pincho, then press the 'Test' button on the UI adjust the memory allotment and click 'Start.'<br />
- [4] Familiarize yourself with the parameters, then when ready press 'Clear' and create your pipeline.<br />
- [5] Press 'Start' to start your created pipeline

# Bulk Runs
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
- [1] Open terminal<br />
- [2] Change working directory to the parent folder of your working folders (i.e. NGS_Project_1)<br />
- [3] $Pincho_conf.py<br />
- [4] Familiarize yourself with the parameters, then when ready press 'Clear' and create your pipeline.<br />
- [5] Press 'Start' to start your created pipeline

# Troubleshooting
The most common source of stress for the program will be the thread count and memory count. Please make sure to have at least 12 threads and 64GB of ram on the systema and please also scale in ratio (3:16). Having more threads than ram in this ratio will sometimes cause subprograms to fail.<br />
<br />
For any other problem please reach out to us directly.

# Subprograms
As Pincho is a methods workflow, subprograms used are stated and given credit below:<br />
- ABySS version 2.2.4 https://github.com/bcgsc/abyss/releases/tag/2.2.4
- Augustus version 3.3.3 https://github.com/Gaius-Augustus/Augustus/releases/tag/v3.3.3
- BBMap version 38.86 https://sourceforge.net/projects/bbmap/files/
- BinPacker version 1.0 https://sourceforge.net/projects/transcriptomeassembly/files/
- BLAT version 35 https://genome-test.gi.ucsc.edu/~kent/src/
- boost version 1.72.0 https://www.boost.org/users/download/
- BUSCO version 4.0.1 https://gitlab.com/ezlab/busco/-/releases/4.0.1
- BWA version 0.7.17 https://github.com/lh3/bwa/releases/tag/v0.7.17
- CD-HIT version 4.8.1 https://github.com/weizhongli/cdhit/releases/tag/V4.8.1
- Circos version 0.69.9 http://circos.ca/software/download/
- Circos-tools version 0.22 http://mkweb.bcgsc.ca/dev/circos/software/download/utilities/
- DIAMOND version 2.0.9 https://github.com/bbuchfink/diamond/releases/tag/v2.0.9
- FASTX-Toolkit version 0.0.13 http://hannonlab.cshl.edu/fastx_toolkit/download.html
- HISAT2 version 2.1.0 http://daehwankimlab.github.io/hisat2/download/#version-hisat2-210
- HMMER version 3.2.1 http://eddylab.org/software/hmmer/
- IDBA version 1.1.3 https://github.com/loneknightpy/idba/releases/tag/1.1.3
- Jellyfish version 2.3.0 https://github.com/gmarcais/Jellyfish/releases/tag/v2.3.0
- JRE version 1.8.0_231 https://www.java.com/en/download/manual.jsp
- kallisto version 0.46.1 https://pachterlab.github.io/kallisto/download
- MEGAHIT version 1.2.9 https://github.com/voutcn/megahit/releases/tag/v1.2.9
- NCBI-BLAST version 2.3.0+ https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.3.0/
- NCBI-BLAST version 2.10.0+ https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.10.0/
- Oases version 0.2.09 https://github.com/dzerbino/oases/releases/tag/0.2.09
- ORP-transrate version 1.0.3 https://github.com/macmanes-lab/Oyster_River_Protocol/tree/master/software
- PALADIN version 1.4.6 https://github.com/ToniWestbrook/paladin/releases/tag/v1.4.6
- Python version 3.8.0 https://www.python.org/downloads/
- Quorum version 1.1.1 https://github.com/gmarcais/Quorum/releases/tag/v1.1.1
- Rcorrector version 1.0.4 https://github.com/mourisl/Rcorrector/releases/tag/v1.0.4
- RSEM version 1.3.1 https://github.com/deweylab/RSEM/releases/tag/v1.3.1
- Salmon version 1.0.0 https://github.com/COMBINE-lab/salmon/releases/tag/v1.0.0
- Samtools version 1.10 https://github.com/samtools/samtools/releases/tag/1.10
- SeqKit version 0.16.0 https://bioinf.shenwei.me/seqkit/download/
- Shannon_Cpp version 0.4.0 https://github.com/bx3/shannon_cpp/releases/tag/v0.4.0
- SOAPdenovo-Trans version 1.0.4 https://github.com/aquaskyline/SOAPdenovo-Trans/releases/tag/1.0.4
- SPAdes version 3.14.1 https://github.com/ablab/spades/releases/tag/v3.14.1
- SRAtoolkit version 2.11.0 https://github.com/ncbi/sratoolkit
- transabyss version 2.0.1 https://github.com/bcgsc/transabyss/releases/tag/2.0.1
- TranscriptomeAssemblyTools https://github.com/harvardinformatics/TranscriptomeAssemblyTools
- TransLiG version 1.3 https://sourceforge.net/projects/transcriptomeassembly/files/TransLiG/
- transrate version 1.0.3 https://github.com/blahah/transrate/releases/tag/v1.0.3
- Trimmomatic version 0.39 http://www.usadellab.org/cms/?page=trimmomatic
- Trinity version 2.11.0 https://github.com/trinityrnaseq/trinityrnaseq/releases/tag/v2.11.0
- Velvet version 1.2.10 https://www.ebi.ac.uk/~zerbino/velvet/


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
