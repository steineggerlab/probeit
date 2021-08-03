# **Probeit**: probes generator for pathogen detection and genotyping
**Probeit** is a software to generate probes which are for pathogen detection and genotyping. **Probeit** is open source **GPL-licensed** software implemented in python for **Linux** and **MacOS**.
## About
Probeit is a probe set designer for pathogen detection and genotyping.
Probeit provides 2 types of workflows: **posnegset** and **snp**. 
* **posnegset**: probe desinger for pathogen detection 
* **snp**: probe desinger for pathogen genotyping.

## Publication

## Get started
#### install **Probeit** with conda

```
conda install -c bioconda probeit
```

#### use **posnegset** workflow
* It will generate a probe set with sequences included in the positive genome but not in the negative genome
* It demands a positive genome(fasta) and a negative genome(fasta)

```
# to get sample files
wget https://raw.githubusercontent.com/steineggerlab/probeit/master/sample/negative.fa
wget https://raw.githubusercontent.com/steineggerlab/probeit/master/sample/positive.fa
wget https://raw.githubusercontent.com/steineggerlab/probeit/master/sample/posnegset1.fa
wget https://raw.githubusercontent.com/steineggerlab/probeit/master/sample/posnegset2.fa
# to run Probeit
probeit.py posnegset -p sample/positive.fasta -n negative.fasta -o posnegset_output
# to compare result files 
diff posnegset_output/sorted1.fa posnegset1.fa
diff posnegset_output/sorted2.fa posnegset2.fa
```
#### to use **snp** workflow with **amino acid** SNPs
* It will generate a probe set which detect input amino acid SNPs from strain genome.
* It demands a reference genome(fasta), a strain genome(fasta), positions, SNPs and a reference annotation(gff).

```
# to get sample files
wget https://raw.githubusercontent.com/steineggerlab/probeit/master/sample/str.fa
wget https://raw.githubusercontent.com/steineggerlab/probeit/master/sample/ref.fa
wget https://raw.githubusercontent.com/steineggerlab/probeit/master/sample/ref.gff
wget https://raw.githubusercontent.com/steineggerlab/probeit/master/sample/snp_aa1.fa
wget https://raw.githubusercontent.com/steineggerlab/probeit/master/sample/snp_aa2.fa
# to run Probeit
probeit.py snp  -r ref.fasta -s str.fasta  -p "10,11,19,20,21,22" -m "aa:orf1ab:L4715F,aa:S:Q52R,aa:S:E484K,aa:S:Q677H,aa:S:F888L,aa:E:L21F,aa:M:I82T"  -o snp_aa_output -a ref.gff
# to compare result files
diff snp_aa_output/sorted1.fa snp_aa1.fa
diff snp_aa_output/sorted2.fa snp_aa2.fa
```
#### to use **snp** workflow with **nucleotide** SNPs
* It will generate a probe set which detect input necleotide SNPs from strain genome.
* It demands a reference genome(fasta), a strain genome(fasta), positions and SNPs.

```
# to get sample files
wget https://raw.githubusercontent.com/steineggerlab/probeit/master/sample/str.fa
wget https://raw.githubusercontent.com/steineggerlab/probeit/master/sample/ref.fa
wget https://raw.githubusercontent.com/steineggerlab/probeit/master/sample/snp_nt1.fa
wget https://raw.githubusercontent.com/steineggerlab/probeit/master/sample/snp_nt2.fa
# to run Probeit
probeit.py snp  -r ref.fasta -s str.fasta  -p "10,11,19,20,21,22" -m "nt:A21716G,nt:G23011A,nt:G23592C,nt:T24223C,nt:C26304T,nt:T26766C"  -o snp_nt_output
# to compare result files
diff snp_nt_output/sorted1.fa snp_nt1.fa
diff snp_nt_output/sorted2.fa snp_nt2.fa
```

## result
**Probeit** produces two result files: **sorted1.fa** and **sorted2.fa**. **Probeit** desings two kinds of probes: **probe1** and **probe2**. **Probe1** is usually 40 nt long. **Probe1** covers a pattern of certain pathogen or strain so it can be used as a ligation probe in cRASL-seq. **Probe2** is usually 20nt long and covers **Probe1**. Usually, **probe2** is not more than 200 nt away from **probe1**, but does not overlap with **probe1**. When used for cRASL-seq, **probe2** is used as a capture probe.
* **sorted1.fa** is a fasta fromat file and it contains **probe1**.
* **sorted2.fa** is a fasta fromat file and it contains **probe2**.
 
## **posnegset** 
For generating probes for detection **Probeit posnegset** is available. The posnegset  workflow finds sequences that are included in the positive genome but not included in the negative genome.
#### Easy User's Guide
```
python probeit.py posnegset -p positive_genome -n negative_genome -o output_dir [additional opts]
```
#### Required Options
###### **-p/--positive** FASTA file : Option for the genome which **MUST** be covered by the probes.
###### **-n/--negative** FASTA file: Option for the genome which **MUST NOT** to be covered by the probes. 
###### **-o/--output** Dir: Option for the output directory. Because Probeit posnegset workflow makes a new directory with given name, you don't have to create a directory. 

## **snp** 
For genotyping **Probeit snp** is available. The snp workflow extracts sequences containing a snp from a strain genome.
#### Easy User's Guide
```
python probeit.py snp  -r reference_genome -s strain_genome  -p podition_list -m snp_list  -o output_dir -a reference_annotation [additional opts]
```

#### Required Options
###### **-r/--reference** FASTA file: Option for the wildtype(reference) genome. 
###### **-s/--strain** FASTA fileOption for the strain genome.  format** as a parameter.
###### **-p/--positions** COMMA SEPARATED INT ARRAY: Option for The : position list. Positons in the position list indicate the positions of SNPs in the 1st Probes.  
###### **-m/--mutations** COMMA SEPARATED SNP ARRAY: Option for SNP list of the strain. Both amino acid differences and nucleotide differences are allowed. 
###### **-o/--output** DIR: Option for the output directory.Because Probeit snp workflow makes a new directory with given name, you don't have to create a directory. 
###### **-a/--annotation** GFF file : Option for wild-type genome annotation. Only required when using amino acid differences in the -m option.

## Installation from source
#### Dependency

* pandas
* primer3-py
* bedtools
* mmseqs
* seqkit
* genmap
* biopython

#### To get source code of Probeit
```
git clone https://github.com/steineggerlab/probeit.git
```
#### To initiate **setcover**
* **Setcover** is a tool for minimizing probe sets.
* You must compile **setcover** before using **Probeit**.

```
cd probeit/setcover
cd setcover
make
cd ..
```

