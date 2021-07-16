# Probeit: probes generator for pathogen detection and genotyping
Probeit is a software to generate probes which are for pathogen detection and genotyping. Probeit is open source GPL-licensed software implemented in python for Linux and MacOS.



### Publication


### Get started
###### Linux
```
# install probeit
conda install -c https://161832-42372094-gh.circle-artifacts.com/0/tmp/artifacts/packages probeit
# use posnetset workflow
python probeit.py posnegset -p sample/positive.fasta -n sample/negative.fasta -o posnegset_output
# to use snp workflow with amino acid SNPs
python probeit.py snp  -r sample/ref.fasta -s sample/str.fasta  -p "10,11,19,20,21,22" -m "aa:orf1ab:L4715F,aa:S:Q52R,aa:S:E484K,aa:S:Q677H,aa:S:F888L,aa:E:L21F,aa:M:I82T"  -o 
# to use snp workflow with nucleotide SNPs
snp_aa_output -a sample/ref.gff
python probeit.py snp  -r sample/ref.fasta -s sample/str.fasta  -p "10,11,19,20,21,22" -m "nt:A21716G,nt:G23011A,nt:G23592C,nt:T24223C,nt:C26304T,nt:T26766C"  -o snp_nt_output
```
###### MacOS
```
# install probeit
conda install -c https://161837-42372094-gh.circle-artifacts.com/0/tmp/artifacts/packages probeit
# use posnetset workflow
python probeit.py posnegset -p sample/positive.fasta -n sample/negative.fasta -o posnegset_output
# to use snp workflow with amino acid SNPs
python probeit.py snp  -r sample/ref.fasta -s sample/str.fasta  -p "10,11,19,20,21,22" -m "aa:orf1ab:L4715F,aa:S:Q52R,aa:S:E484K,aa:S:Q677H,aa:S:F888L,aa:E:L21F,aa:M:I82T"  -o 
# to use snp workflow with nucleotide SNPs
snp_aa_output -a sample/ref.gff
python probeit.py snp  -r sample/ref.fasta -s sample/str.fasta  -p "10,11,19,20,21,22" -m "nt:A21716G,nt:G23011A,nt:G23592C,nt:T24223C,nt:C26304T,nt:T26766C"  -o snp_nt_output
```


### Installation
#### Installation using conda
###### Linux
```
conda install -c https://161832-42372094-gh.circle-artifacts.com/0/tmp/artifacts/packages probeit
```
###### MacOS
```
conda install -c https://161837-42372094-gh.circle-artifacts.com/0/tmp/artifacts/packages probeit
```
#### Installation from source
###### To get source code of Probeit
```
# move to the directory where you want to install probeit.
git clone https://github.com/steineggerlab/probeit.git
cd probeit
```
###### To create probeit environment
```
conda create -n probeit -c conda-forge -c anaconda -c bioconda pandas entrez-direct fire primer3-py bedtools  mmseqs2 seqkit genmap primer3 biopython
```
###### To initiate setcover
```
cd setcover
make
cd ..
```
##### Before run probeit 
```
conda activate probeit  
```


### How to use
We provide 2 types of workflows: posnegset for detection and snp for genotyping.
#### **Probes generating for detection**
For generating probes for detection Probeit posnegset is available. The posnegset  workflow finds sequences that are included in the positive genome but not included in the negative genome.
##### Easy User's Guide
You can use some data in **probeit/sample** to test Probeit posnegset
```
python probeit.py posnegset -p positive_genome -n negative_genome -o output_dir [additional opts]
```
##### Required Options
###### **-p/--positive** FASTA file : Option for the genome which **MUST** be covered by the probes.
###### **-n/--negative** FASTA file: Option for the genome which **MUST NOT** to be covered by the probes. 
###### **-o/--output** Dir: Option for the output directory. Because Probeit posnegset workflow makes a new directory with given name, you don't have to create a directory. 

#### **Probes generating for genotyping** 
For genotyping Probeit snp is available. The snp workflow extracts sequences containing a snp from a strain genome.
##### Easy User's Guide
You can use some data in **probeit/sample** to test Probeit snp
```
python probeit.py snp  -r reference_genome -s strain_genome  -p podition_list -m snp_list  -o output_dir -a reference_annotation [additional opts]
```

##### Required Options
###### **-r/--reference** FASTA file: Option for the wildtype(reference) genome. 
###### **-s/--strain** FASTA fileOption for the strain genome.  format** as a parameter.
###### **-p/--positions** COMMA SEPARATED INT ARRAY: Option for The : position list. Positons in the position list indicate the positions of SNPs in the 1st Probes.  
###### **-m/--mutations** COMMA SEPARATED SNP ARRAY: Option for SNP list of the strain. Both amino acid differences and nucleotide differences are allowed. 
###### **-o/--output** DIR: Option for the output directory.Because Probeit snp workflow makes a new directory with given name, you don't have to create a directory. 
###### **-a/--annotation** GFF file : Option for wild-type genome annotation. Only required when using amino acid differences in the -m option.
