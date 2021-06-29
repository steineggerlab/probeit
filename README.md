# Probeit: probes generator for pathogen detection and genotyping
Probeit is a software to generate probes which are for pathogen detection and genotyping. Probeit is open source GPL-licensed software implemented in python for Linux and MacOS.

### Publication

### Installation
```sh
# to get source code of probeit 
# move to the directory where you want to install probeit.
git clone https://github.com/steineggerlab/probeit.git
cd probeit

# to create probeit environment
conda create -n probeit -c conda-forge -c anaconda -c bioconda pandas entrez-direct fire primer3-py bedtools  mmseqs2 seqkit genmap primer3 biopython

# to initiate setcover
cd set cover
make
cd ..

# Before run probeit 
conda activate probeit  
# or ...
source activate probeit
```
### Getting started
We provide 2 types of workflows: posnegset for detection and snp for genotyping.
##### Probes generating for detection
For generating probes for detection Probeit posnegset is available. The posnegset  workflow finds sequences that are included in the positive genome but not included in the negative genome.
```sh
python probeit.py posnegset -p [positive genome] -n [negative genome] -o [output directory]
```
###### [Essential]
**Positive genome** is the name of the genome file to be included in the probes, and it should be in fasta format.
**negative genome** is the name of the genome file **NOT** to be included in the probes, and it should be in fasta format.
**Output directory** means the name of directory for output files and Probeit posnegset workflow will make a new directory with that name. 
###### [Optional]
- **--never-cluster-positive** option for use this option if you don't need to cluster positive genome. [NONE]       
- **--probe-len1** option for length of probe1 (ligation probe) (default 40)[INT]       
- **--probe-len2** option for length of probe2 (capture probe) (default 20)[INT]      
- **--probe-error1** option for error allowed in probe1 (ligation probe) (default 0)[INT]            
- **--probe-error2** option for error allowed in probe2 (capture probe) (default 1)[INT]                       
- **--minimizing-covered1** option for how many times should each positive seq be covered by probe1 while minimizing probe1 (default 1)[INT]                                
- **--minimizing-covered2** option for how many times should each probe1 be covered by probe2 while minimizing probe2 (default 1)[INT]                           
- **--minimizing-repeats1** option for randomly N iterations while minimizing probe1 (ligation probe) (default 1)[INT]                       
- **--minimizing-repeats2** option for randomly N iterations while minimizing probe2 (capture probe) (default 10)[INT]                                
- **--minimizing-earlystop-criteria1** option for proportion of positive seq covered by probe1 to earlysotp minimizing. probe1 (default 0.9)[FLOAT]                                  
- **--minimizing-earlystop-criteria2** option for proportion of probe1 covered by probe2 to earlysotp minimizing probe2 (default 0.99)[FLOAT]                                              

##### Probes generating for genotyping. 
For genotyping Probeit snp is available. The snp workflow extracts sequences containing a snp from a strain genome.
```sh
python probeit.py snp  -r [reference genome] -s [strain genome]  -p [position list] -m [mutation list]  -o [output directory] -a [reference annotation]
```
###### [Essential]
**Reference genome** means the name of the. wildtype genome file and it must be in fasta format.
**Strain genome** means the name of the strain genome file and it must be in fasta format.
The positions in the **position list** indicate the positions of snp in the ligation probe. This must be a list of integers. ex) ”10,21,31”
**Mutation list** means the snp list of the strain. snps can be written at the amino acid or nucleotide level. ex) ”aa:orf1ab:L4715F,nt:A21716G”
**Output directory** means the name of directory for output files and Probeit snp workflow will make a new directory with that name. 
**Reference annotation** means the wildtype genome annotation file and it must be in GFF format, and it is only needed when mutations are written at the amino acid level. It is only needed when you input SNPs in amino acid level.
###### [Optional]
-  **--probe-len1** option for length of probe1 (ligation probe) (default 40)[INT]                           
- **--probe-len2** option for length of probe2 (capture probe) (default 20)[INT]                           
- **--probe-error2** option for error allowed in probe2 (capture probe) (default 1)[INT]                                                     
- **--minimizing-covered2** option forhow many times should each probe1 be covered by probe2 while minimizing probe2 (default 1)[INT]                                                     
- **--minimizing-repeats2** option for randomly N iterations while minimizing probe2 (capture probe) (default 10)[INT]                                            
- **--minimizing-earlystop-criteria2** option for proportion of probe1 covered by probe2 to earlysotp minimizing probe2 (default 0.99)[FLOAT]                                                 

