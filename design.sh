#!/bin/bash -ex
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

mkdir -p "$OUTPUT" && cd "$OUTPUT"
mkdir -p input && cd input

# download ref. seq
if [ ! -f "download.done" ]; then
  for TID in $(echo $TAXID|awk -F',' '{for(i = 1; i <= NF; i++){print $i}}'); do 
    QUERY="txid${TID}[Organism:exp] and \"complete genome\"[title]"
    echo "Query NCBI ${QUERY}"
    # download genomes
    esearch -db nuccore -query "${QUERY}" | efetch -format fasta >> genomes.fa
  done
  touch download.done
fi
cd ..

# cluster the sequences to reduce highly abundant strains  
mkdir -p cluster
if [ ! -f "cluster/genomes.clu_rep_seq.fasta" ]; then
  mmseqs easy-linclust input/genomes.fa cluster/genomes.clu cluster/tmp --kmer-per-seq 1000 --min-seq-id 0.97 --cov-mode 1 -c 0.95
fi
# create a genome mapping file (genomeID numericalID)
awk 'BEGIN{cnt=0}/^>/{gsub(">","",$1); print $1"\t"cnt; cnt++}' cluster/genomes.clu_rep_seq.fasta > id.lookup

# search probs against MMDB database
mkdir -p mmseqs
if [ ! -f "mmseqs/mmseqs.search" ]; then
  awk -v problen="${PROBLEN1}" '/^>/{header=$1;} !/^>/{ genomeLen=length($1); for(i = 0; i < genomeLen-problen; i++){ print header"_"(i+1); print substr($0, (i+1), problen); } }' "cluster/genomes.clu_rep_seq.fasta" > "mmseqs/probes.fa"
  mmseqs filtertaxseqdb "${MMDB}" "mmseqs/filterdb"  --taxon-list $(echo "${TAXATOIGNORE}")
  mmseqs easy-search "mmseqs/probes.fa" "mmseqs/filterdb" "mmseqs/mmseqs.search" "mmseqs/tmp" --spaced-kmer-mode 0 --mask 0 -c 0.9 --min-seq-id 0.9 --cov-mode 2 --alignment-mode 4 --search-type 3 --format-output query,target,pident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,taxid,taxname,taxlineage
fi

# sed 's/\(.*\)_/\1      /'
# filter ranges with hits to hits in MMDB
# create bed file for clusters and substract the detected alignments  
mkdir -p filter
seqkit fx2tab "cluster/genomes.clu_rep_seq.fasta" -i -l -n |  awk '{print $1"\t"1"\t"$2}' > "filter/genomes.clu_rep_seq.bed"
bedtools subtract -a "filter/genomes.clu_rep_seq.bed" -b <(awk -v problen="${PROBLEN1}" '{split($1, a, "."); split(a[2], b, "_"); gsub(/_[0-9]+$/,"", $1); print $1"\t"b[2]-1"\t"b[2]-1+problen}' "mmseqs/mmseqs.search") > "filter/crosstaxa.bed"

# check primer3 
awk -F'_' '/^>/{str=$1; for(i=2;i<NF;i++){str=str"_"$i} print str"\t"$NF} !/^>/{print}' "mmseqs/probes.fa" > filter/probes.fa

primer3call=$(cat << 'EOF'
from Bio.SeqIO.FastaIO import SimpleFastaParser
from Bio.Seq import Seq
from Bio.Seq import reverse_complement
import pandas as pd
import primer3 #These are 1000 delta g as calculated by primer3
import numpy as np
import fire

#Create list of problematic sequences
def gc_content(this_oligo):
    gcs = this_oligo.count('G') + this_oligo.count('C')
    return gcs/len(this_oligo)

#Create list of problematic sequences
five_reps = []
for each_seq in ["AAAAA","TTTTT","CCCCC","GGGGG"]:
    five_reps.append(Seq(each_seq))

problem_seqs = five_reps
problem_seqs = list(map(str,problem_seqs))

def has_intrinsic_probs(candidate_oligo):
    for s in problem_seqs:
        if s in candidate_oligo:
            return True
    return False

def run_cli( lig_probes_in = "", cap_probes_in = "", outfile = "designed_probes", save_cut_probes = False,
            min_GC=0.30, max_GC = 0.70,
            homo_dimer_tm_cutoff = 60, hairpin_tm_max = 60,
            min_probe_tm = 40):
    print("Minimum Tm: " + str(min_probe_tm))
    print("Minimum GC percentage: " + str(min_GC))
    print("Maximum GC percentage: " + str(max_GC))
    print("Homodimer maximum Tm: " + str(homo_dimer_tm_cutoff))
    print("Hairpin maximum Tm: " + str(hairpin_tm_max))
    #Capture probes
    cp = []
    if(cap_probes_in != ""):
        cap_probes_out = outfile + "_cap.tsv"
        with open(cap_probes_in) as fasta_file:  # Will close handle cleanly
            identifiers = []
            seqs = []
            for title, sequence in SimpleFastaParser(fasta_file):
                identifiers.append(title.split(None, 1)[0])  # First word is ID
                seqs.append(sequence)
                this_seq = str(reverse_complement(Seq(sequence)))
                cp.append(this_seq)
        cap_probes = pd.DataFrame(list(zip(identifiers, seqs,cp)), columns =['id', 'genome_segment','cp'])
        print(str(len(cap_probes)) + " capture probes inputted")

    #ligation probes
    with open(lig_probes_in) as fasta_file:  # Will close handle cleanly
        identifiers = []
        posStart = []
        posEnd = []
        seqs = []
        rc = []
        p1 = []
        p2 = []
        for title, sequence in SimpleFastaParser(fasta_file):
            split_name = title.split('\t',2)
            identifiers.append(split_name[0])
            if(len(split_name)>1):
                this_posStart = int(split_name[1])
            else:
                this_posStart = 0
            this_posEnd = this_posStart + 40
            posStart.append(this_posStart)
            posEnd.append(this_posEnd)
            seqs.append(sequence)
            this_seq = str(reverse_complement(Seq(sequence)))
            rc.append(this_seq)
            p1.append(this_seq[0:20])
            p2.append(this_seq[20:40])
    lig_probes = pd.DataFrame(list(zip(identifiers, posStart,posEnd,seqs,rc,p1,p2)), columns =['id', 'chromStart','chromEnd','genome_segment','rc','p1','p2'])
    print(str(len(lig_probes)) + " ligation probe sets inputted")

    unique_probe_set = set(p1 + p2 + cp)
    unique_probes = pd.DataFrame(list(unique_probe_set), columns = ['p'])
    print("There were " + str(len(unique_probes)) + " unique probes")

    unique_probes = (unique_probes.assign(ultimate_base = unique_probes['p'].str[19])
                     .assign(penultimate_base = unique_probes['p'].str[18]))
    print("Ultimate and penultimate bases assigned")
    unique_probes = (unique_probes.assign(tm = np.vectorize(primer3.calcTm)(unique_probes['p'])))
    unique_probes['hairpin_tm'] = list(map(primer3.calcHairpinTm, unique_probes['p']))
    unique_probes = (unique_probes.assign(homodimer_tm = np.vectorize(primer3.calcHomodimerTm)(unique_probes['p'])))
    unique_probes = (unique_probes.assign(intrinsic_probs = np.vectorize(has_intrinsic_probs)(unique_probes['p'])))
    print("Thermodynamic calculations complete")
    unique_probes = unique_probes.assign(GC_perc = np.vectorize(gc_content)(unique_probes['p']))
    print("GC perc assigned")
    lig_probes_calc = pd.merge(lig_probes, unique_probes.add_prefix('p1_'), how='left', left_on='p1', right_on='p1_p')
    lig_probes_calc = pd.merge(lig_probes_calc, unique_probes.add_prefix('p2_'), how='left', left_on='p2', right_on='p2_p')

    sorder = ['p1_hairpin_tm','p2_hairpin_tm','p1_homodimer_tm','p2_homodimer_tm']
    ascending = [True,True,True,True]
    out_statement = "Ligation probes output: " + outfile + '.tsv'

    if(cap_probes_in != ""):
        cap_probes_calc = pd.merge(cap_probes, unique_probes.add_prefix('cp_'), how='left', left_on='cp', right_on='cp_p')
        cap_probes_filtered = (cap_probes_calc
                       .query('cp_intrinsic_probs == False').query('cp_tm >= ' + str(min_probe_tm))
                       .query('cp_GC_perc >= ' + str(min_GC) + ' & cp_GC_perc <= ' + str(max_GC))
                       .query('cp_homodimer_tm <= ' + str(homo_dimer_tm_cutoff)).query('cp_hairpin_tm <= ' + str(hairpin_tm_max))
                        )
        cap_probes_sorted = cap_probes_filtered.sort_values(by = ['cp_hairpin_tm','cp_homodimer_tm'], ascending = [True,True])
        cap_probes_sorted.to_csv(cap_probes_out, sep = '\t')
        print('Capture probes that passed all filters: ' + str(len(cap_probes_filtered)))
        out_statement = out_statement + "\nCapture probes output: " + cap_probes_out

    lig_probes_filtered = (lig_probes_calc
                   .query('p1_intrinsic_probs == False & p2_intrinsic_probs == False')
                   .query('p1_tm >= ' + str(min_probe_tm) + ' & p2_tm >= ' + str(min_probe_tm))
                   .query('p1_GC_perc >= ' + str(min_GC) + ' & p2_GC_perc >= ' +  str(min_GC) + ' & p1_GC_perc <= ' + str(max_GC) + ' & p2_GC_perc <= ' +  str(max_GC))
                   .query('p1_homodimer_tm <= ' + str(homo_dimer_tm_cutoff) + ' & p2_homodimer_tm <= ' +  str(homo_dimer_tm_cutoff))
                   .query('p1_hairpin_tm <= ' + str(hairpin_tm_max) + ' & p2_hairpin_tm <= ' +  str(hairpin_tm_max))
                    )
    print('Ligation probe sets that passed all filters: ' + str(len(lig_probes_filtered)))
    lig_probes_sorted = lig_probes_filtered.sort_values(by = sorder, ascending = ascending)
    lig_probes_filtered.rename(columns={'id':'chrom'})[['chrom','chromStart','chromEnd']].to_csv(outfile +'_bed.bed',index = False,sep = '\t')
    lig_probes_sorted.to_csv((outfile + '_lig.tsv'), sep = '\t')
    if(save_cut_probes):
        lig_probes_calc.to_csv((outfile + '_uncut.tsv'), sep = '\t')

    return(out_statement)

if __name__ == '__main__':
  fire.Fire(run_cli)
EOF
)

if [ ! -f "filter/primer3_bed.bed" ]; then
  python -c "$primer3call" --lig_probes_in=filter/probes.fa --outfile=filter/primer3 --save_cut_probes=True
fi
# combine two filters: cross tax. hits and primer3 
bedtools subtract -a filter/genomes.clu_rep_seq.bed -b filter/primer3_bed.bed > filter/primer3.neg.bed
bedtools subtract -a filter/crosstaxa.bed -b filter/primer3.neg.bed > filter/crosstaxa.primer3.bed

# compute mappability
mkdir -p mapping
if [ ! -f "mapping/done" ]; then
  genmap index -F <(awk '/^>/{print $1} !/^>/{print}' cluster/genomes.clu_rep_seq.fasta) -I index
  genmap map --no-reverse-complement -E $ERRORINPROB -S filter/crosstaxa.primer3.bed --csv -K ${PROBLEN1} -t -b --frequency-large -I index -O mapping
  touch mapping/done
fi
# remove duplicate k-mers and skips first line
awk -F";" 'NR==1{next}{n=split($2, b, "|"); s=$1";"; delete found; for(i=1; i<=n; i++){ split(b[i], e, ","); if(!(e[1] in found)){ s=s""b[i]"|"; } found[e[1]]=1; } print substr(s, 1, length(s)-1); }' mapping/*.genmap.csv > mapping/no.double.entry.csv
awk -F";" '{n=split($2, b, "|"); pg[$1]=1; prevK=0; for(i = 1; i<=n && prevK==0; i++){ if(b[i]!=$1 && b[i] in pg){ prevK=1 }} if(prevK == 0){ print} }' mapping/no.double.entry.csv > mapping/uniq.genmap.csv

# setcover the k-mers  
mkdir -p setcover
if [ ! -f "setcover/result" ]; then
  setcover mapping/uniq.genmap.csv $COVERAGE "$PROBLEN1" 0.9 > "setcover/result"
fi

buildFastaForShortProb=$(cat << 'EOF'
BEGIN{
    cnt = 1;
}
FNR == NR {
    split($1, header, " ");
    seqname[NR-1]=header[1];
    seq[NR-1]=$2;
    next
}
{
    split($0,entry,";");
    n=split(entry[2], a, "|");
    for(i = 1; i <= n; i++){
        split(a[i], b,",");
        probeStart=b[2] + 1;
        probeEnd=probeStart + problen;
        probPos[cnt] = b[1]":"probeStart":"probeEnd
        cnt++;
    }
}
END{
    for(key in probPos){
        n=split(probPos[key], values, ":");
        seqId1=values[1];
        probeStart=values[2];
        probeEnd=values[3];
        seqLen=length(seq[seqId1]);
        start=(probeStart - window > 1) ? probeStart - window : 1;
        end=(probeStart+problen+window < seqLen) ? probeStart+problen+window : seqLen;
        delete maskMap;
        for(key2 in probPos){
            n=split(probPos[key2], values, ":");
            seqId2=values[1];
            if(seqId1 == seqId2){
                probeStartMask=values[2];
                probeEndMask=values[3];
                for(pos = probeStartMask; pos < probeEndMask; pos++){
                    maskMap[pos]=1;
                }
            }
        }
        str=""
        split(seq[seqId1], seqChar, "");
        for(pos=start; pos <= end ; pos++) {
            if(pos in maskMap){
                str=str"N"
            } else {
                str = str""seqChar[pos]
            }
        }
        print ">"seqname[seqId1]":"probeStart-1":"probeEnd-1"\t"start"\t"end"\t"seqLen;
        print str;
    }
}
EOF
)

mkdir -p input_20
awk -F'\t' -vproblen=$PROBLEN1 -vwindow=200 "$buildFastaForShortProb" <(seqkit fx2tab cluster/genomes.clu_rep_seq.fasta) setcover/result > "input_20/seq.fa"
genmap index -F input_20/seq.fa -I index_20
awk 'BEGIN{cnt=0}/^>/{gsub(">","",$1); print $1"\t"cnt; cnt++}' "input_20/seq.fa" > "id_20.lookup"

# compute mappability
mkdir -p mapping_20
genmap map --no-reverse-complement -E 1 --csv -K $PROBLEN2 -t -b --frequency-large -I index_20 -O mapping_20
# remove duplicate k-mers and skips first line
awk -F";" 'NR==1{next}{n=split($2, b, "|"); s=$1";"; delete found; for(i=1; i<=n; i++){ split(b[i], e, ","); if(!(e[1] in found)){ s=s""b[i]"|"; } found[e[1]]=1; } print substr(s, 1, length(s)-1); }' mapping_20/*.genmap.csv > mapping_20/no.double.entry.csv
awk -F";" '{n=split($2, b, "|"); pg[$1]=1; prevK=0; for(i = 1; i<=n && prevK==0; i++){ if(b[i]!=$1 && b[i] in pg){ prevK=1 }} if(prevK == 0){ print} }'  mapping_20/no.double.entry.csv > mapping_20/uniq.genmap.csv
mkdir -p setcover_20

# minimize 20-mers
setcover mapping_20/uniq.genmap.csv 1 1 0.99 > "setcover_20/result"

# extract probes
awk -vproblen="$PROBLEN1" 'FNR==NR{f[$2]=$1; next}  {split($0,a,";"); split(a[1],b,","); print f[b[1]]"\t"b[2]"\t"b[2]+problen"\t"a[2] }' "id.lookup" "setcover/result" > "setcover/result.bed"
seqkit subseq --quiet --bed "setcover/result.bed" "cluster/genomes.clu_rep_seq.fasta" > "prob40.fa"
# check if probe capute expected stuff
while read i; do awk '{n=split($2, s, ""); cnt = 0; for(i=0; i<=n; i++){if(s[i]=="|") {cnt++}} print cnt+1; } ' <(echo $i); echo $(grep -c $(echo $i|awk '{print $3}' ) cluster/genomes.clu_rep_seq.fasta) ; done < <(seqkit fx2tab prob40.fa) | awk '{first=$1; getnextline; if($1!=first){print "ERROR"}}'

awk -vproblen="$PROBLEN2" 'FNR==NR{f[$2]=$1; next} {split($0,a,";"); split(a[1],b,",");  print f[b[1]]"\t"b[2]"\t"b[2]+problen"\t"a[2] }' "id_20.lookup" "setcover_20/result" > "setcover_20/result.bed"
seqkit subseq --quiet --bed "setcover_20/result.bed" "input_20/seq.fa" > "prob20.fa"
while read i; do awk '{n=split($2, s, ""); cnt = 0; for(i=0; i<=n; i++){if(s[i]=="|") {cnt++}} print cnt+1; } ' <(echo $i); echo $(grep -c $(echo $i|awk '{print $3}' ) input_20/seq.fa) ; done < <(seqkit fx2tab prob20.fa) | awk '{first=$1; getnextline; if($1!=first){print "ERROR"}}'
