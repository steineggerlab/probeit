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
**negative genome**  is the name of the genome file **NOT** to be included in the probes, and it should be in fasta format.
**Output directory** means the name of directory for output files and Probeit posnegset workflow will make a new directory with that name. 

##### Probes generating for genotyping 
For genotyping Probeit snp is available. The snp workflow extracts sequences containing a snp from a strain genome.
```sh
python probeit.py snp  -r [reference genome] -s [strain genome]  -p [position list] -m [mutation list]  -o [output directory] -a [reference annotation]
```
###### [Essential]
**Reference genome** means the name of the wildtype genome file and it must be in fasta format.
**Strain genome** means the name of the strain genome file and it must be in fasta format.
The positions in the **position list** indicate the positions of snp in the ligation probe. This must be a list of integers. ex) ”10,21,31”
**Mutation list** means the snp list of the strain. snps can be written at the amino acid or nucleotide level. ex) ”aa:orf1ab:L4715F,nt:A21716G”
**Output directory** means the name of directory for output files and Probeit snp workflow will make a new directory with that name. 
###### [Optional]
**Reference annotation** means the wildtype genome annotation file and it must be in GFF format, and it is only needed when mutations are written at the amino acid level.

