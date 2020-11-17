from Bio.SeqIO.FastaIO import SimpleFastaParser
from Bio.Seq import Seq
from Bio.Seq import reverse_complement
import pandas as pd
import primer3
import numpy as np
import sys
import os
import shutil
import re
import subprocess
import getopt

# Create list of problematic sequences
problem_seqs =  ["AAAAA","TTTTT","CCCCC","GGGGG"]
five_reps = [Seq(each_seq) for each_seq in problem_seqs]
# method for help message printing
def help():
    print("REQUIRED")
    print(" -p|--positive sequence set that should be covered")
    print(" -n|--negative sequence set that should be not contained")
    print(" -o|--output result output folder")
    print("OPTIONAL")
    print(" --seq-id-cluster clustering identity treshold (default 0.97)")
    print(" --seq-id-probe identity treshold to filter probes aligned to neg. set (default 0.90)")
    print(" --cluster cluster sequences (default 1)")
    print(" --probe-error1 error allowed in probe 1 (default 0)")
    print(" --probe-error2 error allowed in probe 2 (default 1)")
    print(" --probe-len1 length of first probe (default 40)")
    print(" -m compute k-mer conservation with N mismatches (default 0)")
    print(" -c genome coverage by probes (default 1)")
    quit()

# method for making directory
def make_dir(dir_name):
    if not os.path.exists(dir_name):
        os.makedirs(dir_name)

# method for file copying
def copy_file(original, copy):
    if not os.path.exists(copy):
        shutil.copy(original, copy)
        print(original+' is copied to '+os.getcwd()+os.path.sep+copy)

# method for making lookup file
def lookup_maker(fasta_file,out_file):
    if os.path.exists(out_file):
        return
    
    with open(fasta_file) as f:
        with open(out_file,'w') as w:
            cnt = 0
            for title,seq in SimpleFastaParser(f):
                header = title.split()[0].strip()
                w.write(header+'\t'+str(cnt)+'\n')
                cnt+=1
                    
# method for making no_double_entry file
def no_double_entry_maker(source,output):
    if os.path.exists(output):
           return
    for file in os.listdir(source):
        if file[-11:]=='.genmap.csv':
            with open(source+file) as f:
                with open(output,'w') as w:
                    is_header = True
                    for line in f:
                        if is_header:
                            is_header = False
                            continue
                        c1 = line.split(';')[0]
                        s = c1+';'
                        c2 = line.split(';')[1]
                        b = c2.split('|')
                        n=len(b)
                        found=[]
                        for i in range(n):
                            e = b[i].split(',')
                            if e[0] not in found:
                                s = s+b[i].strip()+"|"
                                found.append(e[0])
                        w.write(s[:-1]+'\n')
                        
# method for making uniq_genmap file
def uniq_genmap_maker(no_double_entry,uniq_genmap):
    if os.path.exists(uniq_genmap):
           return
    with open(no_double_entry) as f:
        with open(uniq_genmap,'w') as w:
            pg = []
            for line in f:
                c1 = line.split(';')[0]
                c2 = line.split(';')[1]
                b = c2.split('|')
                n = len(b)
                pg.append(c1)
                prevK = False
                for i in range(n):
                    if b[i].strip()!=c1 and b[i].strip() in pg:
                        prevK = True
                        break
                if not prevK:
                    w.write(line)

#method for result_bed file
def result_bed_maker(lookup,setcover_result,outfile, problen):
    if os.path.exists(outfile):
        return
    f=dict()
    with open(lookup) as f1:
        for line in f1:
            c1 = line.split()[0].strip()
            c2 = line.split()[1].strip()
            f[c2] = c1
    with open(setcover_result) as f2:
        with open(outfile,'w') as w:
            for line in f2:
                a=line.split(';')
                b=a[0].split(',')
                w.write(f[b[0]]+'\t'+b[1]+'\t'+str(int(b[1])+problen)+'\t'+a[1].strip()+'\n')
# PRIMER3CALL
#method for calculating gc ratio
def gc_content(this_oligo):
    gcs = this_oligo.count('G') + this_oligo.count('C')
    return gcs/len(this_oligo)

# method that return true if candidate aligo contains "AAAAA","TTTTT","CCCCC","GGGGG"
def has_intrinsic_probs(candidate_oligo):
    for s in problem_seqs:
        if s in candidate_oligo:
            return True
    return False

def primer3call(
        # input data for ligation probes: filter/probes.fa
        lig_probes_in = "",
        # input data for capture probes
        cap_probes_in = "",
        # name of out file
        outfile = "designed_probes",
        #if save_cut_probes is True then a dataframe lig_probes_calc will be saved as a csv file.
        save_cut_probes = False,
        #cutoff values
        min_GC=0.30,
        max_GC = 0.70,
        homo_dimer_tm_cutoff = 60,
        hairpin_tm_max = 60,
        min_probe_tm = 40,
        #probe size
        probe_size = 40
    ):
    
# PRINTOUT CUT-OFF VALUES
    print("Minimum Tm: " + str(min_probe_tm))
    print("Minimum GC percentage: " + str(min_GC))
    print("Maximum GC percentage: " + str(max_GC))
    print("Homodimer maximum Tm: " + str(homo_dimer_tm_cutoff))
    print("Hairpin maximum Tm: " + str(hairpin_tm_max))

# IF STATEMENT FOR cp and cap_probes
# cp        : list to contains capture probes sequence
# cap_probes: dataframe contains identifiers and capture probes sequence
    cp = []
    if(cap_probes_in != ""): # if input for cap probes exsists then
        cap_probes_out = outfile + "_cap.tsv" #<-default: designed_probes_cap.tsv
        with open(cap_probes_in) as fasta_file:  # Will close handle cleanly
            # list to contains identifiers of cap_probes_in fasta file
            identifiers = []
            # list to contains sequences of cap_probes_in fasta file
            seqs = []
            for title, sequence in SimpleFastaParser(fasta_file):
                identifiers.append(title.split(None, 1)[0])  # First word is ID
                seqs.append(sequence)
                this_seq = str(reverse_complement(Seq(sequence)))
                cp.append(this_seq)
        #TO MAKE cap_probes
        cap_probes = pd.DataFrame(list(zip(identifiers, seqs,cp)), columns =['id', 'genome_segment','cp'])
        print(str(len(cap_probes)) + " capture probes inputted")

# p1, p2 and lig_probes
# p1 and p2 : list for ligation probes sequences
# rc        : list for complementary sequence
# lig_probes: dataframe that contains identifier, startpos, endpos, sequence, complement sequence, p1 and p2 of lig_probes_in fasta file
    with open(lig_probes_in) as fasta_file:  # Will close handle cleanly
        # list for identifiers of lig_probes_in fasta file
        identifiers = []
        # list for start position  and end position
        posStart = []
        posEnd = []
        # list to contains sequence of lig_probes_in fasta file
        seqs = []
        rc = []
        p1 = []
        p2 = []
        for title, sequence in SimpleFastaParser(fasta_file):
            split_name = title.split('\t',2)
            identifiers.append(split_name[0])
            
            # IF STATEMENT FOR this_posStart
            if(len(split_name)>1):
                this_posStart = int(split_name[1])
            else:
                this_posStart = 0
            
            this_posEnd = this_posStart + probe_size
            posStart.append(this_posStart)
            posEnd.append(this_posEnd)
            seqs.append(sequence)
            this_seq = str(reverse_complement(Seq(sequence)))
            rc.append(this_seq)
            # p1 and p2
            mid_pos = round(probe_size/2)
            p1.append(this_seq[0:mid_pos])
            p2.append(this_seq[mid_pos:probe_size])
    # TO MAKE lig_probes
    lig_probes = pd.DataFrame(list(zip(identifiers, posStart,posEnd,seqs,rc,p1,p2)), columns =['id', 'chromStart','chromEnd','genome_segment','rc','p1','p2'])
    print(str(len(lig_probes)) + " ligation probe sets inputted")
    
# unique_probe_set and unique_probes
# unique_probe_set: common set of p1, p2 and cp
# unique_probes   : probes from unique_probe_set + thermodynamic values
    unique_probe_set = set(p1 + p2 + cp) #unique_probe_set
    unique_probes = pd.DataFrame(list(unique_probe_set), columns = ['p'])
    print("There were " + str(len(unique_probes)) + " unique probes")
    #unique_probes.ultimate_base   : the last base
    #unique_probes.penultimate_base: the second to last base
    unique_probes = (unique_probes.assign(ultimate_base = unique_probes['p'].str[-1]).assign(penultimate_base = unique_probes['p'].str[-2]))
    print("Ultimate and penultimate bases assigned")
    #unique_probes.tm tm value of probes
    unique_probes = (unique_probes.assign(tm = np.vectorize(primer3.calcTm)(unique_probes['p'])))
    #unique_probes.hairpin_tm hairpin tm value od probes
    unique_probes['hairpin_tm'] = list(map(primer3.calcHairpinTm, unique_probes['p']))
    #unique_probes.homodimer_tm homo dimer tm value od probes
    unique_probes = (unique_probes.assign(homodimer_tm = np.vectorize(primer3.calcHomodimerTm)(unique_probes['p'])))
    #unique_probes.intrinsic_probs if the probe contains "AAAAA","TTTTT","CCCCC","GGGGG" then true
    unique_probes = (unique_probes.assign(intrinsic_probs = np.vectorize(has_intrinsic_probs)(unique_probes['p'])))
    print("Thermodynamic calculations complete")
     #unique_probes.GC_perc gc content of the probes
    unique_probes = unique_probes.assign(GC_perc = np.vectorize(gc_content)(unique_probes['p']))
    print("GC perc assigned")
    
# lig_probes_calc: MERGED DF with lig_probes and unique_probes, which contains thermodynamic values
    #TO MAKE lig_probes_calc
    lig_probes_calc = pd.merge(lig_probes, unique_probes.add_prefix('p1_'), how='left', left_on='p1', right_on='p1_p')
    lig_probes_calc = pd.merge(lig_probes_calc, unique_probes.add_prefix('p2_'), how='left', left_on='p2', right_on='p2_p')
    
    # VALUES FOR SORTING
    sorder = ['p1_hairpin_tm','p2_hairpin_tm','p1_homodimer_tm','p2_homodimer_tm']
    ascending = [True,True,True,True]
    # RETURN STRING
    out_statement = "Ligation probes output: " + outfile + '.tsv'
    
# IF CONITIONAL STATEMENT FOR cap_probes_calc, cap_probes_filtered and cap_probes_sorted
# ONLY EXCUTED WHEN cap_probes_in IS NOT ""
# cap_probes_calc    : MERGED DF with cap_probes and unique_probes, which contains thermodynamic values
# cap_probes_filtered: FILTERED cap_probes_calc with thermodynamic thresholds
# cap_probes_sorted  : SORTED cap_probes_filtered
    if(cap_probes_in != ""):  # if input for cap probes exsists then
        # TO MAKE cap_probes_calc
        cap_probes_calc = pd.merge(cap_probes, unique_probes.add_prefix('cp_'), how='left', left_on='cp', right_on='cp_p')
        # TO MAKE cap_probes_filtered
        cap_probes_filtered = (
            cap_probes_calc
            .query('cp_intrinsic_probs == False')
            .query('cp_tm >= ' + str(min_probe_tm))
            .query('cp_GC_perc >= ' + str(min_GC) + ' & cp_GC_perc <= ' + str(max_GC))
            .query('cp_homodimer_tm <= ' + str(homo_dimer_tm_cutoff))
            .query('cp_hairpin_tm <= ' + str(hairpin_tm_max))
        )
        # TO MAKE cap_probes_sorted
        cap_probes_sorted = cap_probes_filtered.sort_values(by = ['cp_hairpin_tm','cp_homodimer_tm'], ascending = [True,True])
        # cap_probes_sorted to CSV
        cap_probes_sorted.to_csv(cap_probes_out, sep = '\t')
        print('Capture probes that passed all filters: ' + str(len(cap_probes_filtered)))
        out_statement = out_statement + "\nCapture probes output: " + cap_probes_out

# lig_probes_filtered: FILTERED lig_probes_calc with thermodynamic thresholds
# lig_probes_sorted  : SORTED lig_probes_filtered
    # TO MAKE lig_probes_filtered
    lig_probes_filtered = (
        lig_probes_calc
        .query('p1_intrinsic_probs == False & p2_intrinsic_probs == False')
        .query('p1_tm >= ' + str(min_probe_tm) + ' & p2_tm >= ' + str(min_probe_tm))
        .query('p1_GC_perc >= ' + str(min_GC) + ' & p2_GC_perc >= ' +  str(min_GC) + ' & p1_GC_perc <= ' + str(max_GC) + ' & p2_GC_perc <= ' +  str(max_GC))
        .query('p1_homodimer_tm <= ' + str(homo_dimer_tm_cutoff) + ' & p2_homodimer_tm <= ' +  str(homo_dimer_tm_cutoff))
        .query('p1_hairpin_tm <= ' + str(hairpin_tm_max) + ' & p2_hairpin_tm <= ' +  str(hairpin_tm_max))
    )
    print('Ligation probe sets that passed all filters: ' + str(len(lig_probes_filtered)))
    # TO MAKE lig_probes_sorted
    lig_probes_sorted = lig_probes_filtered.sort_values(by = sorder, ascending = ascending)
    # lig_probes_filtered to CSV
    lig_probes_filtered.rename(columns={'id':'chrom'})[['chrom','chromStart','chromEnd']].to_csv(outfile +'_bed.bed',index = False,sep = '\t')
    # lig_probes_sorted to CSV
    lig_probes_sorted.to_csv((outfile + '_lig.tsv'), sep = '\t')
    
    # OPTIONAL: sace dataframe lig_probes_calc as csv file
    if(save_cut_probes):
        # if true then save lig_probes_calc as tsv file
        lig_probes_calc.to_csv((outfile + '_uncut.tsv'), sep = '\t')
    print(out_statement)

# ARGUMENTS PROCESSING
args = sys.argv[1:]
optlist,args = getopt.getopt(args,'hp:n:o:c:',['positive=','negative=','probe-len1=','probe-len2=','output=','seq-id-cluster=','seq-id-probe=','probe-error1','probe-error2=','help'])
not_get_p = True
not_get_n = True
not_get_o = True
PROBLEN1=40 # PROBLEN1=40
PROBLEN2=20 # PROBLEN2=20
ERRORINPROB1=0 # ERRORINPROB1=0
ERRORINPROB2=1 # ERRORINPROB2=1
COVERAGE=1 # COVERAGE=1
CLUSTER=1 # CLUSTER=1
SEQID=0.97 # SEQID="0.97"
SEQIDPROBE=0.90 # SEQIDPROBE="0.90"
INPUT = ""
NEGATIVE = ""
OUTPUT = ""
POSITIONAL = []

for opt, value in optlist:
    try:
        if opt in ('-h','--help'):
            help()
        elif opt in ('-p','--positive'):
            INPUT = str(value)
            not_get_p = False
        elif opt in ('-n','--negative'):
            NEGATIVE = str(value)
            not_get_n = False
        elif opt == '--probe-len1':
            PROBLEN1 = int(value)
        elif opt == '--probe-len2':
            PROBLEN2 = int(value)
        elif opt in ('-o','--output'):
            OUTPUT = str(value)
            not_get_o = False
        elif opt == '--seq-id-cluster':
            SEQID = float(value)
        elif opt == '--seq-id-probe':
            SEQIDPROBE = float(value)
        elif opt == '--probe-error1':
            ERRORINPROB1 = int(value)
        elif opt == '--probe-error2':
            ERRORINPROB2 = int(value)
        elif opt == '-c':
            COVERAGE = int(value)
        elif opt == '--cluster':
            CLUSTER = int(value)
    except:
        help()
    
# CODE FOR ESSENTIAL ARGUMENTS -p -n -o
# PRINT OUT ALERTS AND INFO
if not_get_p or not_get_n or not_get_o :
    help()
# ABSOLUTE PATH of INPUT(-p)
ABSINPUT = os.getcwd()+os.path.sep+INPUT
# ABSOLUTE PATH of NEGATIVE(-n)
ABSNEGATIVE = os.getcwd()+os.path.sep+NEGATIVE

# make a working directory named like OUTPUT(-o)
make_dir(OUTPUT)
# go into working directory
os.chdir(OUTPUT)
# make a new directory named 'input'
make_dir('input')
# make a copy file into 'input' directory
copy_file(ABSINPUT,'input/genomes.fa')

# cluster the sequences to reduce highly abundant strains
# make a new directory 'cluster'
make_dir('cluster')

# if CLUSTER equals 1 then do mmesqs easy-lincluster with input/genomes.fa(copy of inputfile -p)
# TO MAKE cluster/genomes.clu_rep_seq.fasta
if CLUSTER == 1 :
    print("cluster input sequences")
    if not os.path.exists("cluster/genomes.clu_rep_seq.fasta" ):
        os.system('mmseqs easy-linclust input/genomes.fa cluster/genomes.clu cluster/tmp -v 0 --kmer-per-seq-scale 0.5 --kmer-per-seq 1000 --min-seq-id {} --cov-mode 1 -c 0.95'.format(SEQID))
# if CLUSTER doesn't equal 1 then just make a copy of input file(-p)
else:
    copy_file(ABSINPUT, 'cluster/genomes.clu_rep_seq.fasta')
#mkae a lookup file of cluster/genomes.clu_rep_seq.fasta
lookup_maker('cluster/genomes.clu_rep_seq.fasta','id.lookup')

                
# search probs against MMDB database

# make a new directory 'mmseqs'
make_dir('mmseqs')

# if mmseqs/mmseqs.search dosen't exist then make a file 'mmseqs/probes.fa' and then make a file "mmseqs/mmseqs.search" using mmseqseasy-search between 'mmseqs/probes.fa' and negative(-n)
# mmseqs/probes.fa will be used to make filter/probes.fa
if not os.path.exists("mmseqs/mmseqs.search"):
    print("remove probes that aligns to negative set")
    with open('cluster/genomes.clu_rep_seq.fasta') as f:
        with open('mmseqs/probes.fa','w') as w:
            for title,seq in SimpleFastaParser(f):
                header = '>'+title.split()[0].strip()
                genomeLen=len(seq.strip())
                for i in range(genomeLen-PROBLEN1):
                    w.write(header+'_'+str(i+1)+'\n')
                    w.write(seq[i:i+PROBLEN1]+'\n')
    
    os.system('mmseqs easy-search "mmseqs/probes.fa" {} "mmseqs/mmseqs.search" "mmseqs/tmp" -v 0 --spaced-kmer-mode 0 -k 13 --mask 0 -c 0.9 --min-seq-id {} --cov-mode 2 --alignment-mode 4 --search-type 3 --format-output query,target,pident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue'.format(ABSNEGATIVE,SEQIDPROBE))

# filter ranges with hits to hits in MMDB
# create bed file for clusters and substract the detected alignments  

make_dir('filter')
with open('cluster/genomes.clu_rep_seq.fasta') as f:
    with open('filter/genomes.clu_rep_seq.bed','w') as w:
        for title,seq in SimpleFastaParser(f):
            header = title.split()[0].strip()
            w.write(header + '\t1\t'+str(len(seq.strip()))+'\n')

with open('mmseqs/mmseqs.search') as f:
    with open('mmseqs/mmseqs.txt','w')as w:
        for line in f:
            c1 = line.split()[0].strip()
            c1 = re.sub('_[0-9]+$','',c1)
            b=c1.split('_')
            n=len(b)
            w.write(c1+'\t'+str(int(b[n-1])-1+PROBLEN1)+'\n')

os.system('bedtools subtract -a "filter/genomes.clu_rep_seq.bed" -b mmseqs/mmseqs.txt > filter/crosstaxa.bed')


# TO MAKE A FILE 'filter/probes.fa'
# filter/probes.fa: file used as lig_probes_in
with open('mmseqs/probes.fa') as f:
    with open('filter/probes.fa','w') as w:
        for title,seq in SimpleFastaParser(f):
            header = ('>'+title).split('_')
            w.write('_'.join(header[:-1])+'\t'+header[-1]+'\n')
            w.write(seq+'\n')
                
# To call primer3call
if not os.path.exists("filter/primer3_bed.bed"):
    print("filter probes based on primer3")
    primer3call(lig_probes_in='filter/probes.fa', outfile='filter/primer3', save_cut_probes=True, probe_size = PROBLEN1)
    
# combine two filters: cross tax. hits and primer3
os.system('bedtools subtract -a filter/genomes.clu_rep_seq.bed -b filter/primer3_bed.bed > filter/primer3.neg.bed')
os.system('bedtools subtract -a filter/crosstaxa.bed -b filter/primer3.neg.bed > filter/crosstaxa.primer3.bed')

if not os.path.exists("mapping_probe1/done"):
    print('compute mappability')
    make_dir('mapping_probe1')
    with open('cluster/genomes.clu_rep_seq.fasta') as f:
        with open('cluster/genmap.clu_rep_seq.fasta','w') as w:
            for title,seq in SimpleFastaParser(f):
                w.write(('>'+title).split()[0].strip()+'\n')
                w.write(seq+'\n')
            
    os.system("genmap index -F cluster/genmap.clu_rep_seq.fasta -I index_probe1 > /dev/null")
    os.system('genmap map --no-reverse-complement -E {} -S filter/crosstaxa.primer3.bed --csv -K {} -t -b --frequency-large -I index_probe1 -O mapping_probe1 > /dev/null'.format(ERRORINPROB1,PROBLEN1))
    #touch mapping_probe1/done
    with open('mapping_probe1/done','w') as w:
        w.write('')

no_double_entry_maker('mapping_probe1/','mapping_probe1/no.double.entry.csv')
uniq_genmap_maker('mapping_probe1/no.double.entry.csv','mapping_probe1/uniq.genmap.csv')
# setcover the k-mers
#mkdir -p setcover_probe1
make_dir('setcover_probe1')
if not os.path.exists('setcover_probe1/result'):
    print("minimize probe set")
    os.system('../setcover/setcover -c {} -l {} -p 0.9 -d 11 -i 1 mapping_probe1/uniq.genmap.csv cluster/genomes.clu_rep_seq.fasta > "setcover_probe1/result"'.format(COVERAGE,PROBLEN1))

# extract probes
# TO MAKE setcover_probe1/result.bed

result_bed_maker('id.lookup', 'setcover_probe1/result','setcover_probe1/result.bed',PROBLEN1)

#FINAL OUTPUT FILE1 setcover_probe1/result.bed cluster/genomes.clu_rep_seq.fasta > prob40.fa
if not os.path.exists('prob40.fa'):
    os.system('seqkit subseq --quiet --bed "setcover_probe1/result.bed" "cluster/genomes.clu_rep_seq.fasta" > "prob40.fa"')
# 'prob'+str(PROBLEN1)+'.fa'
# 'prob1.fa'
# check if probe capute expected stuff

#If the -r option is given, this signifies 'raw' input, and backslash escaping is disabled.

# DEBUGGING
# while read -r i; 
#     do awk '{
#         n=split($2, s, ""); 
#         cnt = 0; 
#         for(i=0; i<=n; i++){
#             if(s[i]=="|") {
#                 cnt++
#             }
#         } 
#         print cnt+1; 
#     } ' <(echo "$i"); 
#     echo $(grep -c $(echo "$i"|awk '{print $3}' ) cluster/genomes.clu_rep_seq.fasta) ; 
#     done < <(seqkit fx2tab prob40.fa) | awk '{first=$1; getnextline; if($1!=first){print "ERROR"}}'

# PROBE2
make_dir('input_probe2')
if PROBLEN2 != -1:
    cnt=1
    with open('cluster/genomes.clu_rep_seq.fasta') as f1:
        seqname=dict()
        seq=dict()
        nr=0
        for title,sequence in SimpleFastaParser(f1):
            header=('>'+title).split()
            seqname[nr]=header[0]
            seq[nr]=sequence
            nr+=1
    with open('setcover_probe1/result') as f2:
        lines = f2.readlines()
        probPos = dict()
        for line in lines:
            entry=line.split(';')
            a=entry[1].split('|')
            n=len(a)
            for i in range(n):
                b=a[i].split(',')
                probeStart=int(b[1])+1
                probeEnd=probeStart+PROBLEN1
                probPos[cnt]=b[0]+':'+str(probeStart)+':'+str(probeEnd)
                cnt+=1
    with open('input_probe2/seq.fa','w') as w:
        for key in probPos:
            values = probPos[key].split(':')
            n = len(values)
            seqId1 = int(values[0])
            probeStart = int(values[1])
            probeEnd = int(values[2])
            seqLen = len(''.join(seq[seqId1].split()))
            start = probeStart-200 if probeStart > 200 else 1
            end =  probeStart+PROBLEN1+200 if (probeStart+PROBLEN1+200<seqLen) else seqLen
            maskMap = []
            for key2 in probPos:
                values=probPos[key2].split(':')
                n = len(values)
                seqId2 = int(values[0])
                if seqId1==seqId2:
                    probeStartMask = int(values[1])-1
                    probeEndMask = int(values[2])-1
                    for pos in range(probeStartMask,probeEndMask):
                        maskMap.append(pos)

            STR=''
            seqChar=''.join(seq[seqId1].split())
            for pos in range(start-1,end):
                if pos in maskMap:
                    STR=STR+'N'
                else:
                    STR=STR+seqChar[pos]

            w.write(seqname[seqId1]+':'+str(probeStart-1)+':'+str(probeEnd-1)+'\t'+str(start)+'\t'+str(end)+'\t'+str(seqLen)+'\n')
            w.write(STR+'\n')
            
    os.system('genmap index -F input_probe2/seq.fa -I index_probe2 > /dev/null')
    lookup_maker('input_probe2/seq.fa','id_probe2.lookup')

    # compute mappability
    # mkdir -p mapping_probe2
    make_dir('mapping_probe2')
    os.system('genmap map --no-reverse-complement -E {} --csv -K {} -t -b --frequency-large -I index_probe2 -O mapping_probe2 > /dev/null'.format(ERRORINPROB2,PROBLEN2))

    # remove duplicate k-mers and skips first line
    no_double_entry_maker('mapping_probe2/','mapping_probe2/no.double.entry.csv')
    uniq_genmap_maker('mapping_probe2/no.double.entry.csv','mapping_probe2/uniq.genmap.csv')
    # mkdir -p setcover_probe2
    make_dir('setcover_probe2')

    # minimize 20-mers
    os.system('../setcover/setcover -c 1 -l 1 -p 0.99 -d 20 -i 10 mapping_probe2/uniq.genmap.csv input_probe2/seq.fa > "setcover_probe2/result"')
    result_bed_maker('id_probe2.lookup', 'setcover_probe2/result','setcover_probe2/result.bed',PROBLEN2)
    #FINAL OUTPUT FILE 2
    os.system('seqkit subseq --quiet --bed "setcover_probe2/result.bed" "input_probe2/seq.fa" > "prob2.fa"')

# DEBUGGING
# while read -r i;
# do awk '{
#     n=split($2, s, "");
#     cnt = 0;
#     for(i=0; i<=n; i++){
#         if(s[i]=="|") {
#             cnt++
#         }
#     }
#     print cnt+1;
# } '<(echo "$i");
# echo $(grep -c $(echo "$i"|awk '{print $3}' ) input_probe2/seq.fa) ;
# done < <(seqkit fx2tab prob2.fa) | awk '{
#     first=$1;
#     getnextline;
#     if($1!=first){
#         print "ERROR"
#     }
# }'

