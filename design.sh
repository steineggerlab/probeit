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

# download ref. seq
for TID in $(echo $TAXID|awk -F',' '{for(i = 1; i <= NF; i++){print $i}}'); do 
  QUERY="txid${TID}[Organism:exp] and \"complete genome\"[title]"
  echo "Query NCBI ${QUERY}"
  # download genomes
  esearch -db nuccore -query "${QUERY}" | efetch -format fasta >> genomes.fa
done
cd ..

# cluster 
mmseqs easy-linclust input/genomes.fa input/genomes.clu  input/tmp --kmer-per-seq 1000 --min-seq-id 0.97 --cov-mode 1 -c 0.95

# id to genomes lookup
awk 'BEGIN{cnt=0}/^>/{gsub(">","",$1); print $1"\t"cnt; cnt++}' input/genomes.clu_rep_seq.fasta > id.lookup

# compute mappability
genmap index -F input/genomes.clu_rep_seq.fasta -I index
mkdir mapping
genmap map -E $ERRORINPROB --csv -K $PROBLEN -t -b --frequency-large -I index -O mapping

# setcover the k-mers  
mkdir setcover
setcover mapping/genomes.clu_rep_seq.genmap.csv $COVERAGE > setcover/result
awk -vproblen="$PROBLEN" 'FNR==NR{f[$2]=$1; next}  {split($0,a,";"); split(a[1],b,","); print f[b[1]]"\t"b[2]"\t"b[2]+problen }' id.lookup setcover/result > setcover/result.bed 
seqkit subseq --quiet --bed "setcover/result.bed" "input/genomes.clu_rep_seq.fasta" > "prob.fa"

# reformat genmap output
#awk -F';' 'BEGIN{getline; total= NF-1} {split($1,a,","); pos=a[2]; cons="1"; cnt=0; for(i = 2; i <= NF; i++) { if(length($i) != 0) { cnt++; cons=cons";"i} } print pos"\t"cnt/total"\t"cons}' "mapping/genomes.genmap.csv" > "mapping/kmers.conserv"

# get all tax. ids
#esearch -db taxonomy -query "txid${TAXID}[Organism:exp]" | efetch -format docsum | xtract -pattern DocumentSummary -element TaxId > taxids
#mmseqs easy-search "prob.fa" $MMDB mmseqs.search tmp --spaced-kmer-mode 0 --search-type 3 --format-output query,target,pident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,taxid,taxname,taxlineage
#awk -v taxid=$TAXID -F'\t' '$12 != taxid {print}' mmseqs.search > mmseqs.search.notax
#awk 'FNR==NR{f[$1]=$12";"f[$1]; next } {print $0"\t"f[$1]}' <(cat mmseqs.search.notax) <(seqkit fx2tab  prob.fa ) > prob.tsv

#time blastn -query "prob.fa" -db nt -remote -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore staxids" > blast.nt.search
#awk 'NR==FNR{str=$2; for(i=3; i <=NF; i++) { str=str"\t"$i;} f[">"$1]=str; next} /^>/ {printf("%s\t%s\n", $0, f[$1]); getline; print}' <(awk -F'\t' '!($1":"$13 in a) {f[$1]=$13"\t"f[$1]; a[$1":"$13]=1} END{for(i in f){print i"\t"f[i];}}' blast.nt.search) "prob.fa"

