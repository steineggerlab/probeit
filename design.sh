#!/bin/bash -ex

PARAMCNT=0
PROBLEN=40
ERRORINPROB=0

while [[ $# -gt 0 ]]
do
key="$1"
case $key in
    -t|--taxid)
    TAXID="$2"
    ((PARAMCNT++))
    shift # past argument
    shift # past value
    ;;
    -l|--len)
    PROBLEN="$2"
    shift # past argument
    shift # past value
    ;;
    -o|--output)
    OUTPUT="$2"
    ((PARAMCNT++))
    shift # past argument
    shift # past value
    ;;
    -d|--db)
    MMDB="$(greadlink -f $2)"
    ((PARAMCNT++))
    shift # past argument
    shift # past value
    ;;
    -e)
    ERRORINPROB="$2"
    shift # past argument
    shift # past value
    ;;
    *)    # unknown option
    POSITIONAL+=("$1") # save it in an array for later
    shift # past argument
    ;;
esac
done
if [[ "$PARAMCNT" -lt 3 ]]; then
    echo "Usage: design.sh -t taxid -l prob_len -o output_folder"

    echo "REQUIRED"
    echo " -t|--taxid NCBI taxonomical identifier"
    echo " -o|--output result output folder"
    echo " -d|--db path to db"
    echo "OPTIONAL"
    echo " -l|--len probe length (default 40)"
    echo " -e compute k-mer conservation with N errors (default 0)"
    exit
fi

set -- "${POSITIONAL[@]}" # restore positional parameters

mkdir "$OUTPUT" && cd "$OUTPUT"
mkdir input && cd input

QUERY="txid${TAXID}[Organism:exp] AND refseq[filter]"
echo "Query NCBI ${QUERY}"
# download ref. seq
FILES=$(esearch -db "nucleotide" -query "${QUERY}" | efetch -format fasta | awk '/^>/{s=++d".fasta"} {print > s}END{print "FILES="d}'|awk -F'=' '/FILES=/{print $2}')
QUERY="txid${TAXID}[Organism:exp] and \"complete genome\"[title]"
echo "Query NCBI ${QUERY}"

# download genomes
esearch -db nuccore -query "${QUERY}" | efetch -format fasta | awk -v d=$FILES '/^>/{d++;s=d".fasta"} {print > s}'
cd ..

# create index for leading sequence
seqkit faidx -f input/1.fasta
REFSEQID=$(awk -F'\t' '{print $1}' input/1.fasta.seqkit.fai)
REFSEQID2=$(awk '{print $1}' input/1.fasta.seqkit.fai)
SEQLEN=$(awk -F'\t' '{print $2}' input/1.fasta.seqkit.fai)
awk -v ref="$REFSEQID" -v len="$SEQLEN" 'BEGIN{print ref"\t"1"\t"len; exit}' > input/1.fasta.bed

# compute mappability
genmap index -FD input -I index
mkdir mapping
genmap map -S input/1.fasta.bed -E $ERRORINPROB --csv -K $PROBLEN -t -b --frequency-large  -I index -O mapping

# reformat genmap output
awk -F';' 'BEGIN{getline; total= NF-1} {split($1,a,","); pos=a[2]; cons="1"; cnt=0; for(i = 2; i <= NF; i++) { if(length($i) != 0) { cnt++; cons=cons";"i} } print pos"\t"cnt/total"\t"cons}' "mapping/1.genmap.csv" > "mapping/kmers.conserv"

#awk -F' ' '!/^>/{for(i=1; i <= NF; i++){print i"\t"$i}}' "mapping/1.genmap.txt" > "mapping/kmers.conserv"
# perform k-means to pick distant k-mers
#Rscript -e 'd=read.table("mapping/kmers.conserv",h=1);d$C=kmeans(d,'$K')$cluster;d;' | awk '/^[0-9]/{print}' | sort -k3,3nr  > "mapping/top.kmers.cluster"
#sort -k2,2nr -k1,1n "mapping/kmers.conserv" | awk -vLEN=40 'BEGIN{ pos=-(LEN+1); }  (pos + LEN) < $1 || (pos - LEN) > $1 {print; pos=$1}' | head -n $K > "mapping/top.kmers"
#awk -vK=$K 'NF==4 && !($4 in f){print $2"\t"$3"\t"$4; f[$4]=1;}' "mapping/top.kmers.cluster" > "mapping/top.kmers.rep"
IFS=$'\n'
CNT=0

# create bed file 
awk -v ref="$REFSEQID2" -v problen="$PROBLEN" '{print ref"\t"($1-1)"\t"($1+problen-1)"\t"$2" "$3}' mapping/kmers.conserv >  mapping/kmers.conserv.bed
seqkit subseq --quiet --bed mapping/kmers.conserv.bed "input/1.fasta" > "prob.fa"

# get all tax. ids
#esearch -db taxonomy -query "txid${TAXID}[Organism:exp]" | efetch -format docsum | xtract -pattern DocumentSummary -element TaxId > taxids
mmseqs easy-search "prob.fa" $MMDB mmseqs.search tmp --spaced-kmer-mode 0 --search-type 3 --format-output query,target,pident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,taxid,taxname,taxlineage
awk -v taxid=$TAXID -F'\t' '$12 != taxid {print}' mmseqs.search > mmseqs.search.notax
awk 'FNR==NR{f[$1]=$12";"f[$1]; next } {print $0"\t"f[$1]}' <(cat mmseqs.search.notax) <(seqkit fx2tab  prob.fa ) > prob.tsv

#time blastn -query "prob.fa" -db nt -remote -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore staxids" > blast.nt.search
#awk 'NR==FNR{str=$2; for(i=3; i <=NF; i++) { str=str"\t"$i;} f[">"$1]=str; next} /^>/ {printf("%s\t%s\n", $0, f[$1]); getline; print}' <(awk -F'\t' '!($1":"$13 in a) {f[$1]=$13"\t"f[$1]; a[$1":"$13]=1} END{for(i in f){print i"\t"f[i];}}' blast.nt.search) "prob.fa"

