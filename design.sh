#!/bin/bash -e
export PATH=/Users/mad/Documents/workspace/probeit/setcover:/Users/mad/Documents/workspace/mmseqs2/build/src:/Users/mad/Documents/workspace/genmap-build/bin:$PATH
PARAMCNT=0
PROBLEN1=40
PROBLEN2=20
ERRORINPROB=0
COVERAGE=1 
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
    PROBLEN1="$2"
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
    TAXATOIGNORE="$2"
    shift # past argument
    shift # past value
    ;;
    -m)
    ERRORINPROB="$2"
    shift # past argument
    shift # past value
    ;;
    -c)
    COVERAGE="$2"
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
    echo " -e exclude taxa from cross reaction check (default '')"
    echo " -m compute k-mer conservation with N mismatches (default 0)"
    echo " -c genome coverage by probes (default 1)"
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

# cluster the sequences to reduce highly abundant strains  
mkdir cluster
mmseqs easy-linclust input/genomes.fa cluster/genomes.clu cluster/tmp --kmer-per-seq 1000 --min-seq-id 0.97 --cov-mode 1 -c 0.95

# create a genome mapping file (genomeID numericalID)
awk 'BEGIN{cnt=0}/^>/{gsub(">","",$1); print $1"\t"cnt; cnt++}' cluster/genomes.clu_rep_seq.fasta > id.lookup

# search probs against MMDB database
mkdir mmseqs
awk -v problen="${PROBLEN1}" '/^>/{header=$1;} !/^>/{ genomeLen=length($1); for(i = 0; i < genomeLen-problen; i++){ print header"_"(i+1); print substr($0, (i+1), problen); } }' "cluster/genomes.clu_rep_seq.fasta" > "mmseqs/probes.fa"
mmseqs filtertaxseqdb "${MMDB}" "mmseqs/filterdb"  --taxon-list $(echo "!${TAXATOIGNORE}")
mmseqs easy-search "mmseqs/probes.fa" "mmseqs/filterdb" "mmseqs/mmseqs.search" "mmseqs/tmp" --spaced-kmer-mode 0 --mask 0 -c 0.9 --min-seq-id 0.9 --cov-mode 2 --alignment-mode 4 --search-type 3 --format-output query,target,pident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,taxid,taxname,taxlineage

# filter ranges with hits to hits in MMDB
# create bed file for clusters and substract the detected alignments  
mkdir filter
seqkit fx2tab "cluster/genomes.clu_rep_seq.fasta" -i -l -n |  awk '{print $1"\t"1"\t"$2}' > "filter/genomes.clu_rep_seq.bed"
bedtools subtract -a "filter/genomes.clu_rep_seq.bed" -b <(awk -v problen="${PROBLEN1}" '{split($1, a, "."); split(a[2], b, "_"); gsub(/_[0-9]+$/,"", $1); print $1"\t"b[2]"\t"b[2]+problen}' "mmseqs/mmseqs.search") > "filter/crosstaxa.bed"

# compute mappability
mkdir mapping
genmap index -F <(awk '/^>/{print $1} !/^>/{print}' cluster/genomes.clu_rep_seq.fasta) -I index
genmap map -E $ERRORINPROB -S filter/crosstaxa.bed --csv -K ${PROBLEN1} -t -b --frequency-large -I index -O mapping

# remove duplicate k-mers and skips first line
awk -F";" 'NR==1{next}{n=split($2, b, "|"); pg[$1]=1; prevK=0; for(i = 1; i<=n && prevK==0; i++){ if(b[i]!=$1 && b[i] in pg){ prevK=1 }} if(prevK == 0){ print} }' mapping/*.genmap.csv > mapping/uniq.genmap.csv

# setcover the k-mers  
mkdir setcover
setcover mapping/uniq.genmap.csv $COVERAGE "$PROBLEN1" > "setcover/result"

# compute mappability
mkdir mapping_20
genmap map -E $ERRORINPROB -S filter/crosstaxa.bed --csv -K 20 -t -b --frequency-large -I index -O mapping_20

buildOverlappingSets=$(cat << 'EOF'
NR==FNR{
   split($1, b, ",");
   genome=b[1];
   probCnt++;
   probePos=b[2];
   n=split($2 ,a, "|");  
   delete ids;
   for(p = probePos - 100; p < probePos+100; p++){
           if(p >= (probePos-problen2) && p <= probePos + problen1 || ignore[genome","p][idx][b[1]] == 1){
	      ignore[genome","p][idx][b[1]]=1;
              continue;
           }  
	   probesAtPos[genome][p]++;
	   idx = probesAtPos[genome][p];
	   probeName[genome][p][idx]=probCnt;
   	   for(i = 1; i<=n; i++) 
	   { 
	     split(a[i], b, ",");
	     lookup[genome","p][idx][b[1]]=1;
	   }
   }
   next;
}
$1 in lookup{
   split($1, b, ",");
   genome=b[1];
   probePos=b[2];
   probesAtPosCnt = 0
   for(i in lookup[$1]) probesAtPosCnt++
   n=split($2 , a, "|"); 
   for(probeIdx = 1; probeIdx <= probesAtPosCnt; probeIdx++){
      split($1, c, ",");
      str=probeName[genome][probePos][probeIdx]","c[2]";"
      str=str""a[i]
      for(i = 2; i<=n; i++) {
         split(a[i], b, ",");
         if(lookup[genome","probePos][probeIdx][b[1]] == 1){
            str=str"|"a[i]
         }
      }
      print str
   }
}
EOF
)

awk -vproblen1="$PROBLEN1" -vproblen2="$PROBLEN2" -F';' "$buildOverlappingSets" setcover/result mapping_20/*.genmap.csv 


awk -vproblen="$PROBLEN1" 'FNR==NR{f[$2]=$1; next}  {split($0,a,";"); split(a[1],b,","); print f[b[1]]"\t"b[2]"\t"b[2]+problen }' "id.lookup" "setcover/result" > "setcover/result.bed"
seqkit subseq --quiet --bed "setcover/result.bed" "cluster/genomes.clu_rep_seq.fasta" > "prob.fa"
