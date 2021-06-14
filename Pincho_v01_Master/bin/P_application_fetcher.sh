tar -xf orp-transrate.tar.xz
tar -xf transrate-1.0.3-linux-x86_64.tar.xz
tar -xf velvet_1.2.10.tar.xz
tar -xf busco.tar.xz
rm orp-transrate.tar.xz transrate-1.0.3-linux-x86_64.tar.xz velvet_1.2.10.tar.xz busco.tar.xz
cd blat
chmod +x blat
cd ..

cd seqkit_linux_amd64
chmod +x seqkit
cd ..

wget http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.39.zip
unzip Trimmomatic-0.39.zip
rm Trimmomatic-0.39.zip

wget https://github.com/mourisl/Rcorrector/archive/refs/tags/v1.0.4.tar.gz
tar -xvzf v1.0.4.tar.gz
rm v1.0.4.tar.gz
cd Rcorrector-1.0.4
make
chmod +x run_rcorrector.pl
cd ..

wget https://github.com/bcgsc/abyss/releases/download/2.2.4/abyss-2.2.4.tar.gz
tar -xvzf abyss-2.2.4.tar.gz
rm abyss-2.2.4.tar.gz
cd abyss-2.2.4
./autogen.sh
./configure
make
sudo make install
cd ..

wget https://sourceforge.net/projects/bbmap/files/BBMap_38.86.tar.gz/download
tar -xvzf download
rm download
mkdir BBMap_38.86
mv -v bbmap* BBMap_38.86

wget https://sourceforge.net/projects/transcriptomeassembly/files/BinPacker_1.0.tar.gz/download
tar -xvzf download
rm download
cd BinPacker_1.0
./configure
make
cd ..

wget https://sourceforge.net/projects/transcriptomeassembly/files/TransLiG/TransLiG_1.3.tar.gz/download
tar -xvzf download
rm download
cd TransLiG_1.3
./configure
make
cd ..

wget https://github.com/loneknightpy/idba/releases/download/1.1.3/idba-1.1.3.tar.gz
tar -xvzf idba-1.1.3.tar.gz
rm idba-1.1.3.tar.gz
cd idba-1.1.3
./configure
make
cd ..

wget https://github.com/voutcn/megahit/releases/download/v1.2.9/MEGAHIT-1.2.9-Linux-x86_64-static.tar.gz
tar -xvzf MEGAHIT-1.2.9-Linux-x86_64-static.tar.gz
rm MEGAHIT-1.2.9-Linux-x86_64-static.tar.gz

wget https://github.com/dzerbino/oases/archive/refs/tags/0.2.09.tar.gz
tar -xvzf 0.2.09.tar.gz
rm 0.2.09.tar.gz
mv oases-0.2.09 oases
mkdir oases-0.2.09
mv oases* oases-0.2.09
cp -a velvet_1.2.10/. oases-0.2.09/oases/velvet
cd oases-0.2.09/oases
make 'MAXKMERLENGTH=149' 'CATEGORIES=5'
cd ..

wget https://github.com/bx3/shannon_cpp/releases/download/v0.4.0/shannon_cpp-Linux64-v0.4.0.tar.gz
tar -xvzf shannon_cpp-Linux64-v0.4.0.tar.gz
rm shannon_cpp-Linux64-v0.4.0.tar.gz

wget https://github.com/ablab/spades/releases/download/v3.14.1/SPAdes-3.14.1-Linux.tar.gz
tar -xvzf SPAdes-3.14.1-Linux.tar.gz
rm SPAdes-3.14.1-Linux.tar.gz

wget https://github.com/bcgsc/transabyss/releases/download/2.0.1/transabyss-2.0.1.zip
unzip transabyss-2.0.1.zip
rm transabyss-2.0.1.zip

wget https://github.com/trinityrnaseq/trinityrnaseq/releases/download/v2.11.0/trinityrnaseq-v2.11.0.FULL.tar.gz
tar -xvzf trinityrnaseq-v2.11.0.FULL.tar.gz
rm trinityrnaseq-v2.11.0.FULL.tar.gz
mkdir trinityrnaseq-v2.11.0.FULL
mv -v trinityrnaseq-v2.11.0* trinityrnaseq-v2.11.0.FULL

wget https://github.com/weizhongli/cdhit/releases/download/V4.8.1/cd-hit-v4.8.1-2019-0228.tar.gz
tar -xvzf cd-hit-v4.8.1-2019-0228.tar.gz
rm cd-hit-v4.8.1-2019-0228.tar.gz

wget https://cloud.biohpc.swmed.edu/index.php/s/hisat2-210-Linux_x86_64/download
unzip download
rm download
mkdir hisat2-2.1.0-Linux_x86_64
mv hisat2-2.1.0* hisat2-2.1.0-Linux_x86_64

wget https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.3.0/ncbi-blast-2.3.0+-x64-win64.tar.gz
tar -xvzf ncbi-blast-2.3.0+-x64-win64.tar.gz
rm ncbi-blast-2.3.0+-x64-win64.tar.gz
mkdir ncbi-blast-2.3.0+-x64-linux
mv ncbi-blast-2.3.0+* ncbi-blast-2.3.0+-x64-linux

wget https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.10.0/ncbi-blast-2.10.0+-x64-linux.tar.gz
tar -xvzf ncbi-blast-2.10.0+-x64-linux.tar.gz
rm ncbi-blast-2.10.0+-x64-linux.tar.gz
mkdir ncbi-blast-2.10.0+-x64-linux
mv ncbi-blast-2.10.0+* ncbi-blast-2.10.0+-x64-linux

wget https://github.com/pachterlab/kallisto/releases/download/v0.46.1/kallisto_linux-v0.46.1.tar.gz
tar -xvzf kallisto_linux-v0.46.1.tar.gz
rm kallisto_linux-v0.46.1.tar.gz
mkdir kallisto_linux-v0.46.1
mv kallisto* kallisto_linux-v0.46.1

wget https://github.com/deweylab/RSEM/archive/refs/tags/v1.3.1.tar.gz
tar -xvzf v1.3.1.tar.gz
rm v1.3.1.tar.gz

wget https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/2.11.0/sratoolkit.2.11.0-ubuntu64.tar.gz
tar -xvzf sratoolkit.2.11.0-ubuntu64.tar.gz
rm sratoolkit.2.11.0-ubuntu64.tar.gz

wget https://github.com/samtools/samtools/releases/download/1.10/samtools-1.10.tar.bz2
tar -xf samtools-1.10.tar.bz2
rm samtools-1.10.tar.bz2

wget https://github.com/aquaskyline/SOAPdenovo-Trans/archive/refs/tags/1.0.4.tar.gz
tar -xvzf 1.0.4.tar.gz
rm 1.0.4.tar.gz

wget https://github.com/gmarcais/Quorum/releases/download/v1.1.1/quorum-1.1.1.tar.gz
tar -xvzf quorum-1.1.1.tar.gz
rm quorum-1.1.1.tar.gz

wget https://boostorg.jfrog.io/artifactory/main/release/1.72.0/source/boost_1_72_0.tar.gz
tar -xvzf boost_1_72_0.tar.gz
rm boost_1_72_0.tar.
cd boost_1_72_0
./bootstrap.sh
./b2
cd ..

wget https://github.com/lh3/bwa/releases/download/v0.7.17/bwa-0.7.17.tar.bz2
tar -xf bwa-0.7.17.tar.bz2
rm bwa-0.7.17.tar.bz2

wget http://eddylab.org/software/hmmer/hmmer-3.2.1.tar.gz
tar -xvzf hmmer-3.2.1.tar.gz
rm hmmer-3.2.1.tar.gz
mkdir hmmer
mv hmmer-3.2.1* hmmer
cd hmmer/hmmer-3.2.1
./configure
make
cd ..
cd ..

wget https://github.com/gmarcais/Jellyfish/releases/download/v2.3.0/jellyfish-2.3.0.tar.gz
tar -xvzf jellyfish-2.3.0.tar.gz
rm jellyfish-2.3.0.tar.gz

wget https://github.com/Gaius-Augustus/Augustus/releases/download/v3.4.0/augustus-3.4.0.tar.gz
tar -xvzf augustus-3.4.0.tar.gz
rm augustus-3.4.0.tar.gz 
mv augustus-3.4.0 Augustus-master

wget https://github.com/COMBINE-lab/salmon/releases/download/v1.5.1/salmon-1.5.1_linux_x86_64.tar.gz
tar -xvzf salmon-1.5.1_linux_x86_64.tar.gz
rm salmon-1.5.1_linux_x86_64.tar.gz

cd ..
cd lib
tar -xf eukaryota_odb10.2019-11-20.tar.xz
rm eukaryota_odb10.2019-11-20.tar.xz

cd kegg_libraries
tar -xf npr_kegg.tar.xz
tar -xf xla_kegg.tar.xz
rm npr_kegg.tar.xz xla_kegg.tar.xz
cd xla_npr_kegg
cat xla_npr_KEGG_nucleotide_1.fasta xla_npr_KEGG_nucleotide_2.fasta xla_npr_KEGG_nucleotide_3.fasta xla_npr_KEGG_nucleotide_4.fasta > xla_npr_KEGG_nucleotide.fasta
rm xla_npr_KEGG_nucleotide_1.fasta xla_npr_KEGG_nucleotide_2.fasta xla_npr_KEGG_nucleotide_3.fasta xla_npr_KEGG_nucleotide_4.fasta
cd ..
cd ..

cd trembl_amphi
cat trembl_amphibia_1.fasta trembl_amphibia_2.fasta trembl_amphibia_3.fasta trembl_amphibia_4.fasta trembl_amphibia_5.fasta trembl_amphibia_6.fasta trembl_amphibia_7.fasta trembl_amphibia_8.fasta > trembl_amphibia.fasta
rm trembl_amphibia_1.fasta trembl_amphibia_2.fasta trembl_amphibia_3.fasta trembl_amphibia_4.fasta trembl_amphibia_5.fasta trembl_amphibia_6.fasta trembl_amphibia_7.fasta trembl_amphibia_8.fasta 
cat trembl_amphibia_1.phr trembl_amphibia_2.phr trembl_amphibia_3.phr > trembl_amphibia.phr
rm trembl_amphibia_1.phr trembl_amphibia_2.phr trembl_amphibia_3.phr
cat trembl_amphibia_1.psq trembl_amphibia_2.psq trembl_amphibia_3.psq trembl_amphibia_4.psq trembl_amphibia_5.psq trembl_amphibia_6.psq trembl_amphibia_7.psq trembl_amphibia_8.psq > trembl_amphibia.psq
rm trembl_amphibia_1.psq trembl_amphibia_2.psq trembl_amphibia_3.psq trembl_amphibia_4.psq trembl_amphibia_5.psq trembl_amphibia_6.psq trembl_amphibia_7.psq trembl_amphibia_8.psq
cd ..

cd uniprot_db
cat uniprot_sprot_1.fasta uniprot_sprot_2.fasta uniprot_sprot_3.fasta uniprot_sprot_4.fasta uniprot_sprot_5.fasta uniprot_sprot_6.fasta uniprot_sprot_7.fasta uniprot_sprot_8.fasta uniprot_sprot_9.fasta uniprot_sprot_10.fasta uniprot_sprot_11.fasta > uniprot_sprot.fasta
rm uniprot_sprot_1.fasta uniprot_sprot_2.fasta uniprot_sprot_3.fasta uniprot_sprot_4.fasta uniprot_sprot_5.fasta uniprot_sprot_6.fasta uniprot_sprot_7.fasta uniprot_sprot_8.fasta uniprot_sprot_9.fasta uniprot_sprot_10.fasta uniprot_sprot_11.fasta
cat uniprot_db_1.phr uniprot_db_2.phr uniprot_db_3.phr uniprot_db_4.phr uniprot_db_5.phr uniprot_db_6.phr > uniprot_db.phr
rm uniprot_db_1.phr uniprot_db_2.phr uniprot_db_3.phr uniprot_db_4.phr uniprot_db_5.phr uniprot_db_6.phr
cat uniprot_db_1.psq uniprot_db_2.psq uniprot_db_3.psq uniprot_db_4.psq uniprot_db_5.psq uniprot_db_6.psq uniprot_db_7.psq uniprot_db_8.psq uniprot_db_9.psq > uniprot_db.psq
rm uniprot_db_1.psq uniprot_db_2.psq uniprot_db_3.psq uniprot_db_4.psq uniprot_db_5.psq uniprot_db_6.psq uniprot_db_7.psq uniprot_db_8.psq uniprot_db_9.psq
tar -xf taxdb.tar.xz
rm taxdb.tar.xz
cd ..

cp uniprot_db/taxdb* trembl_amphi
cp uniprot_db/taxdb* kegg_libraries/xla_npr_kegg
