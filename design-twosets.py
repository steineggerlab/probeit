from Bio.SeqIO.FastaIO import SimpleFastaParser
from Bio.Seq import Seq
from Bio.Seq import reverse_complement
import pandas as pd
import primer3
import numpy as np
import sys
import os
import shutil
import getopt
import subprocess
# Create list of problematic sequences
problem_seqs = ['A'*5, 'T'*5, 'C'*5, 'C'*5]
five_reps = [Seq(each_seq) for each_seq in problem_seqs]


# method for help message printing
def printUsage():
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


def printAndLog(message):
    if message == '':
        return
    print(message.strip())
    with open('log.txt', 'a')as log:
        log.write(message.strip() + '\n')


def cmdAndLog(commandList, output=''):
    sp = subprocess.Popen(commandList, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = sp.communicate()
    printAndLog(stdout.decode('UTF-8'))
    printAndLog(stderr.decode('UTF-8'))
    if output != '':
        with open(output, 'a') as out:
            out.write(stdout.decode('UTF-8').strip() + '\n')
            out.write(stderr.decode('UTF-8').strip() + '\n')


# method for absoulte path
def makeAbsPath(fileName, directory=''):
    return os.path.abspath(directory+fileName)


# method for making directory
def makeDir(directory):
    if not os.path.exists(directory):
        os.makedirs(directory)
    return directory+os.path.sep


# method for file copying
def copyFile(inputFile, output):
    if not os.path.exists(output):
        shutil.copy(inputFile, output)
        printAndLog(inputFile + ' is copied to ' + output)


# method for making lookup file
def makeLookup(inputFasta, output):
    if os.path.exists(output):
        return output
    with open(inputFasta) as f:
        with open(output, 'w') as w:
            cnt = 0
            for title, seq in SimpleFastaParser(f):
                header = title.split()[0].strip()
                w.write(header + '\t' + str(cnt) + '\n')
                cnt += 1


# method for making no_double_entry file
def makeNoDoubleEntryCSV(directory):
    absOutput = makeAbsPath(directory=directory, fileName='no.double.entry.csv')
    if os.path.exists(absOutput):
        return absOutput
    for file in os.listdir(directory):
        if file[-11:] == '.genmap.csv':
            with open(makeAbsPath(directory=directory, fileName=file)) as f:
                with open(absOutput, 'w') as w:
                    is_header = True
                    for line in f:
                        if is_header:
                            is_header = False
                            continue
                        c1 = line.split(';')[0]
                        s = c1 + ';'
                        c2 = line.split(';')[1]
                        b = c2.split('|')
                        n = len(b)
                        found = []
                        for i in range(n):
                            e = b[i].split(',')
                            if e[0] not in found:
                                s = s + b[i].strip() + "|"
                                found.append(e[0])
                        w.write(s[:-1] + '\n')
    return absOutput


# method for making uniq_genmap file
def makeUniqGenmapCSV(directory, inputCSV):
    absOutput = makeAbsPath(directory=directory, fileName='uniq.genmap.csv')
    if os.path.exists(absOutput):
        return absOutput
    with open(inputCSV) as f:
        with open(absOutput, 'w') as w:
            pg = []
            for line in f:
                c1 = line.split(';')[0]
                c2 = line.split(';')[1]
                b = c2.split('|')
                n = len(b)
                pg.append(c1)
                prevK = False
                for i in range(n):
                    if b[i].strip() != c1 and b[i].strip() in pg:
                        prevK = True
                        break
                if not prevK:
                    w.write(line)
    return absOutput


# method for result_bed file
def makeResultBed(lookup, setcoverResult, output, problen):
    if os.path.exists(output):
        return output
    f = dict()
    with open(lookup) as f1:
        for line in f1:
            c1 = line.split()[0].strip()
            c2 = line.split()[1].strip()
            f[c2] = c1
    with open(setcoverResult) as f2:
        with open(output, 'w') as w:
            for line in f2:
                a = line.split(';')
                b = a[0].split(',')
                w.write(f[b[0]] + '\t' + b[1] + '\t' + str(int(b[1]) + problen) + '\t' + a[1].strip() + '\n')


#
def copyInputFasta(directory, inputFasta, output):
    absOutput = makeAbsPath(directory=directory, fileName=output)
    copyFile(inputFasta, absOutput)
    return absOutput


# method for reducing redundancy using mmseqs cluster
def reduceRedundancy(doCluster, seqIdentity, directory, inputFasta, output, lookup):
    outputFileName = makeAbsPath(directory=directory, fileName=output)
    outputFasta = outputFileName + '_rep_seq.fasta'
    absLookup = makeAbsPath(fileName=lookup)
    if doCluster:
        printAndLog("cluster input sequences")
        command = (
            'mmseqs easy-linclust {} {} cluster/tmp -v 3 --kmer-per-seq-scale 0.5 --kmer-per-seq 1000 --min-seq-id {}' +
            ' --cov-mode 1 -c 0.95'
        )
        if not os.path.exists(outputFasta):
            cmdAndLog(command.format(inputFasta, outputFileName, seqIdentity).split())
        else:
            copyFile(inputFasta, outputFasta)
    makeLookup(outputFasta, absLookup)
    return outputFasta, absLookup


# method for making probes
def makeInitialProbesMaskedProbes(directory, inputFasta, negative, output, maskedOutput, probeLen, seqIdProbe):
    absOutput = makeAbsPath(directory=directory, fileName=output)
    printAndLog("remove probes that aligns to negative set")
    with open(inputFasta) as f:
        with open(absOutput, 'w') as w:
            for title, seq in SimpleFastaParser(f):
                header = '>' + title.split()[0].strip()
                genomeLen = len(seq.strip())
                for i in range(genomeLen - probeLen):
                    w.write(header + '_' + str(i + 1) + '\n')
                    w.write(seq[i:i + probeLen] + '\n')
    absMaskedOutput = makeAbsPath(directory=directory, fileName=maskedOutput)
    if not os.path.exists(absMaskedOutput):
        command = (
            'mmseqs easy-search {} {} {} mmseqs/tmp -v 3 --spaced-kmer-mode 0 -k 13 --mask 0 -c 0.9 --min-seq-id {}' +
            ' --cov-mode 2 --alignment-mode 4 --search-type 3 ' +
            '--format-output query,target,pident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue'
        )
        cmdAndLog(command.format(absOutput, negative, absMaskedOutput, seqIdProbe).split())
    return absOutput, absMaskedOutput


# method for filter/genomes.clu_rep_seq.bed file
def makeRangeFileForFasta(directory, inputFasta, maskResult, maskOutput, problen):
    absOutput = makeAbsPath(directory=directory, fileName=inputFasta.split(os.path.sep)[-1].split('.fasta')[0] + '.bed')
    absTempFile = makeAbsPath(directory=directory, fileName='mmseqs.txt')
    with open(inputFasta) as f:
        with open(absOutput, 'w') as w:
            for title, seq in SimpleFastaParser(f):
                header = title.split()[0].strip()
                w.write(header + '\t1\t' + str(len(seq.strip())) + '\n')
    with open(maskResult) as f:
        with open(absTempFile, 'w')as w:
            for line in f:
                c1 = line.split()[0].strip()
                b = c1.split('_')[-1]
                c1 = c1.split('_')[0]
                w.write(c1 + '\t' + str(int(b) - 1) + '\t' + str(int(b) - 1 + problen) + '\n')
    absMaskOutput = makeSubtractBed(directory=directory, inputBed=absOutput, negativeBed=absTempFile, output=maskOutput)
    return absOutput, absMaskOutput


#
def makeFastaForThermoFilter(directory, inputFasta, output):
    absOutput = makeAbsPath(directory=directory, fileName=output)
    with open(inputFasta) as f:
        with open(absOutput, 'w') as w:
            for title, seq in SimpleFastaParser(f):
                header = ('>' + title).split('_')
                w.write('_'.join(header[:-1]) + '\t' + header[-1] + '\n')
                w.write(seq + '\n')
    return absOutput


#
def makeIndexFastaForComputingMappability(fastaDirectory, computingDirectory, inputFasta, output):
    absDone = makeAbsPath(directory=computingDirectory, fileName='done')
    absOutput = makeAbsPath(directory=fastaDirectory, fileName=output)
    if os.path.exists(absDone):
        return absOutput
    printAndLog('compute mappability')
    with open(inputFasta) as f:
        with open(absOutput, 'w') as w:
            for title, seq in SimpleFastaParser(f):
                w.write(('>' + title).split()[0].strip() + '\n')
                w.write(seq + '\n')
    with open(absDone, 'w') as w:
        w.write('')
    return absOutput


# method for cluster/genmap.clu_rep_seq.fasta file
def computeMappability(directory, inputFasta, indexFile, error, kmer, selector=''):
    command1 = 'genmap index -F {} -I {} > /dev/null'
    command2 = 'genmap map --no-reverse-complement -E {} {} --csv -K {} -t -b --frequency-large -I {} -O {} > /dev/null'
    selectorOpt = '' if selector == '' else '-S '+selector
    os.system(command1.format(inputFasta, indexFile))
    os.system(command2.format(error, selectorOpt, kmer, indexFile, directory))
    noDoubleEntryCSV = makeNoDoubleEntryCSV(directory=directory)
    uniqGenmapCSV = makeUniqGenmapCSV(inputCSV=noDoubleEntryCSV, directory=directory)
    return uniqGenmapCSV


#
def setcover(inputCSV, inputFasta, output, coverage, length, proportion, distance, iteration):
    command = '../setcover/setcover -c {} -l {} -p {} -d {} -i {} {} {}'
    cmdAndLog(
        command.format(coverage, length, proportion, distance, iteration, inputCSV, inputFasta).split(), output=output
    )


#
def makeFinalProbeFasta(inputBed, inputFasta, output):
    absOutput = makeAbsPath(fileName=output)
    os.system('seqkit subseq --quiet --bed "{}" "{}" > "{}"'.format(inputBed, inputFasta, output))
    return absOutput


# method for input_probe2/seq.fa file
def makeProbe2InputFasta(problen, window, directory, inputFasta, inputSetcoverResult, output, lookup):
    absOutput = makeAbsPath(directory=directory, fileName=output)
    absLookup = makeAbsPath(fileName=lookup)
    cnt = 1
    with open(inputFasta) as f1:
        seqname = dict()
        seq = dict()
        nr = 0
        for title, sequence in SimpleFastaParser(f1):
            header = ('>' + title).split()
            seqname[nr] = header[0]
            seq[nr] = sequence
            nr += 1
    with open(inputSetcoverResult) as f2:
        lines = f2.readlines()
        probPos = dict()
        for line in lines:
            entry = line.split(';')
            a = entry[1].split('|')
            n = len(a)
            for i in range(n):
                b = a[i].split(',')
                probeStart = int(b[1]) + 1
                probeEnd = probeStart + problen
                probPos[cnt] = b[0] + ':' + str(probeStart) + ':' + str(probeEnd)
                cnt += 1
    with open(absOutput, 'w') as w:
        for key in probPos:
            values = probPos[key].split(':')
            # n = len(values)
            seqId1 = int(values[0])
            probeStart = int(values[1])
            probeEnd = int(values[2])
            seqLen = len(''.join(seq[seqId1].split()))
            start = probeStart - window if probeStart > window else 1
            end = probeStart + problen + window if (probeStart + problen + window < seqLen) else seqLen
            maskMap = []
            for key2 in probPos:
                values = probPos[key2].split(':')
                # n = len(values)
                seqId2 = int(values[0])
                if seqId1 == seqId2:
                    probeStartMask = int(values[1]) - 1
                    probeEndMask = int(values[2]) - 1
                    for pos in range(probeStartMask, probeEndMask):
                        maskMap.append(pos)

            STR = ''
            seqChar = ''.join(seq[seqId1].split())
            for pos in range(start - 1, end):
                if pos in maskMap:
                    STR = STR + 'N'
                else:
                    STR = STR + seqChar[pos]

            w.write(
                seqname[seqId1] + ':' + str(probeStart - 1) + ':' + str(probeEnd - 1) + '\t' + str(start) + '\t' +
                str(end) + '\t' + str(seqLen) + '\n')
            w.write(STR + '\n')
    makeLookup(absOutput, absLookup)
    return absOutput, absLookup


#
def minimzeProbeSet(
    directory,
    inputCSV,
    inputFasta,
    inputLookup,
    outputResult,
    outputBed,
    coverage,
    length,
    proportion,
    distance,
    iteration,
    probLen
):
    absOutputBed = makeAbsPath(directory=directory, fileName=outputBed)
    absOutputResult = makeAbsPath(directory=directory, fileName=outputResult)
    if os.path.exists(absOutputResult):
        return absOutputResult, absOutputBed
    printAndLog("minimize probe set")
    setcover(
        inputCSV=inputCSV,
        inputFasta=inputFasta,
        output=absOutputResult,
        coverage=coverage,
        length=length,
        proportion=proportion,
        distance=distance,
        iteration=iteration,
    )
    with open(absOutputResult) as f:
        resultList = f.readlines()
    with open(absOutputResult, 'w') as w:
        for line in resultList:
            if ';' in line:
                w.write(line)
    makeResultBed(
        lookup=inputLookup,
        setcoverResult=absOutputResult,
        output=absOutputBed,
        problen=probLen,
    )
    return absOutputResult, absOutputBed


# method for calculating gc ratio
def gc_content(this_oligo):
    gcs = this_oligo.count('G') + this_oligo.count('C')
    return gcs / len(this_oligo)


# method that return true if candidate aligo contains "AAAAA","TTTTT","CCCCC","GGGGG"
def has_intrinsic_probs(candidate_oligo):
    for s in problem_seqs:
        if s in candidate_oligo:
            return True
    return False


#
def filterProbesByThermoProperty(
        directory, inputBed, output,
        # input
        lig_probes_in="", cap_probes_in="",
        # option
        save_cut_probes=False,
        # cutoff values
        minGC=0.30, maxGC=0.70, homo_dimer_tm_cutoff=60, hairpin_tm_max=60, min_probe_tm=40,
        # probe size
        probe_size=40
):
    absOutput = makeAbsPath(directory=directory, fileName=output)
    absPositiveOutput = absOutput + '_bed.bed'
    absNegativeOutput = makeAbsPath(directory=directory, fileName=output+'.neg.bed')

    if os.path.exists(absNegativeOutput):
        return absPositiveOutput, absNegativeOutput
    printAndLog("filter probes based on primer3")

    # PRINTOUT CUT-OFF VALUES
    printAndLog("Minimum Tm: " + str(min_probe_tm))
    printAndLog("Minimum GC percentage: " + str(minGC))
    printAndLog("Maximum GC percentage: " + str(maxGC))
    printAndLog("Homodimer maximum Tm: " + str(homo_dimer_tm_cutoff))
    printAndLog("Hairpin maximum Tm: " + str(hairpin_tm_max))
    # TO MAKE DF FOR CAPTURE PROBES
    cp = []
    if cap_probes_in != "":
        cap_probes_out = absOutput + "_cap.tsv"
        with open(cap_probes_in) as f:
            identifiers = []
            seqs = []
            for title, sequence in SimpleFastaParser(f):
                identifiers.append(title.split(None, 1)[0])  # First word is ID
                seqs.append(sequence)
                this_seq = str(reverse_complement(Seq(sequence)))
                cp.append(this_seq)
        cap_probes = pd.DataFrame(list(zip(identifiers, seqs, cp)), columns=['id', 'genome_segment', 'cp'])
        printAndLog(str(len(cap_probes)) + " capture probes inputted")
    #  TO MAKE DF FOR LIGATION PROBES
    with open(lig_probes_in) as f:
        identifiers = []
        posStart = []
        posEnd = []
        seqs = []
        rc = []
        p1 = []
        p2 = []
        for title, sequence in SimpleFastaParser(f):
            split_name = title.split('\t', 2)
            identifiers.append(split_name[0])
            if len(split_name) > 1:
                this_posStart = int(split_name[1])
            else:
                this_posStart = 0
            this_posEnd = this_posStart + probe_size
            posStart.append(this_posStart)
            posEnd.append(this_posEnd)
            seqs.append(sequence)
            this_seq = str(reverse_complement(Seq(sequence)))
            rc.append(this_seq)
            mid_pos = round(probe_size / 2)
            p1.append(this_seq[0:mid_pos])
            p2.append(this_seq[mid_pos:probe_size])
    lig_probes = pd.DataFrame(list(zip(identifiers, posStart, posEnd, seqs, rc, p1, p2)),
                              columns=['id', 'chromStart', 'chromEnd', 'genome_segment', 'rc', 'p1', 'p2'])
    printAndLog(str(len(lig_probes)) + " ligation probe sets inputted")
    # TO MAKE DF unique_probes CONTAINS SEQUENCES AND THERMODYANMIC DATA OF p1, p2 AND cp WITHOUT REDUNDAUCY
    unique_probe_set = set(p1 + p2 + cp)
    unique_probes = pd.DataFrame(list(unique_probe_set), columns=['p'])
    printAndLog("There were " + str(len(unique_probes)) + " unique probes")
    unique_probes = (unique_probes.assign(ultimate_base=unique_probes['p'].str[-1]).assign(
        penultimate_base=unique_probes['p'].str[-2]))
    printAndLog("Ultimate and penultimate bases assigned")
    unique_probes = (unique_probes.assign(tm=np.vectorize(primer3.calcTm)(unique_probes['p'])))
    unique_probes['hairpin_tm'] = list(map(primer3.calcHairpinTm, unique_probes['p']))
    unique_probes = (unique_probes.assign(homodimer_tm=np.vectorize(primer3.calcHomodimerTm)(unique_probes['p'])))
    unique_probes = (unique_probes.assign(intrinsic_probs=np.vectorize(has_intrinsic_probs)(unique_probes['p'])))
    printAndLog("Thermodynamic calculations complete")
    unique_probes = unique_probes.assign(GC_perc=np.vectorize(gc_content)(unique_probes['p']))
    printAndLog("GC perc assigned")
    # TO MAKE lig_probes_calc
    lig_probes_calc = pd.merge(lig_probes, unique_probes.add_prefix('p1_'), how='left', left_on='p1', right_on='p1_p')
    lig_probes_calc = pd.merge(lig_probes_calc, unique_probes.add_prefix('p2_'), how='left', left_on='p2',
                               right_on='p2_p')
    # VALUES FOR SORTING
    sorder = ['p1_hairpin_tm', 'p2_hairpin_tm', 'p1_homodimer_tm', 'p2_homodimer_tm']
    ascending = [True, True, True, True]
    # RETURN STRING
    out_statement = "Ligation probes output: " + absOutput + '.tsv'
    # TO FILTER CAPTURE PROBES USING THERMODYNAMIC FEATURES AND SAVE FILTERED CAPTURE PROBES AS TSV FILE
    if cap_probes_in != "":
        # TO MAKE cap_probes_calc
        cap_probes_calc = (
            pd.merge(cap_probes, unique_probes.add_prefix('cp_'), how='left', left_on='cp', right_on='cp_p')
        )
        # FILTER
        cap_probes_filtered = (
            cap_probes_calc
            .query('cp_intrinsic_probs == False')
            .query('cp_tm >= ' + str(min_probe_tm))
            .query('cp_GC_perc >= ' + str(minGC) + ' & cp_GC_perc <= ' + str(maxGC))
            .query('cp_homodimer_tm <= ' + str(homo_dimer_tm_cutoff))
            .query('cp_hairpin_tm <= ' + str(hairpin_tm_max))
        )
        # SORT AND SAVE
        cap_probes_sorted = cap_probes_filtered.sort_values(by=['cp_hairpin_tm', 'cp_homodimer_tm'],
                                                            ascending=[True, True])
        cap_probes_sorted.to_csv(cap_probes_out, sep='\t')
        printAndLog('Capture probes that passed all filters: ' + str(len(cap_probes_filtered)))
        out_statement = out_statement + "\nCapture probes output: " + cap_probes_out
    # TO FILTER LIGATION PROBES USING THERMODYNAMIC FEATURES AND SAVE FILTERED LIGATION PROBES AS TSV FILE
    # FILTER
    lig_probes_filtered = (
        lig_probes_calc
        .query('p1_intrinsic_probs == False & p2_intrinsic_probs == False')
        .query('p1_tm >= ' + str(min_probe_tm) + ' & p2_tm >= ' + str(min_probe_tm))
        .query(
            'p1_GC_perc >= ' +
            str(minGC) +
            ' & p2_GC_perc >= ' +
            str(minGC) +
            ' & p1_GC_perc <= ' +
            str(maxGC) +
            ' & p2_GC_perc <= ' +
            str(maxGC)
        )
        .query('p1_homodimer_tm <= ' + str(homo_dimer_tm_cutoff) + ' & p2_homodimer_tm <= ' + str(homo_dimer_tm_cutoff))
        .query('p1_hairpin_tm <= ' + str(hairpin_tm_max) + ' & p2_hairpin_tm <= ' + str(hairpin_tm_max))
    )
    # SORT AND SAVE
    printAndLog('Ligation probe sets that passed all filters: ' + str(len(lig_probes_filtered)))
    lig_probes_sorted = lig_probes_filtered.sort_values(by=sorder, ascending=ascending)
    lig_probes_filtered.rename(columns={'id': 'chrom'})[['chrom', 'chromStart', 'chromEnd']].to_csv(
        absPositiveOutput, index=False, sep='\t')
    lig_probes_sorted.to_csv((absOutput + '_lig.tsv'), sep='\t')
    # OPTIONAL: save dataframe lig_probes_calc as tsv file
    if save_cut_probes:
        lig_probes_calc.to_csv((absOutput + '_uncut.tsv'), sep='\t')
    printAndLog(out_statement)
    absNegativeOutput = (
        makeSubtractBed(
            directory=directory, inputBed=inputBed, negativeBed=absPositiveOutput, output=output + '.neg.bed'
        )
    )
    return absPositiveOutput, absNegativeOutput


#
def makeSubtractBed(directory, inputBed, negativeBed, output):
    absOutput = makeAbsPath(directory=directory, fileName=output)
    os.system('bedtools subtract -a {} -b {} > {}'.format(inputBed, negativeBed, absOutput))
    return absOutput


# MAIN
def main(pars):
    # INITIATION VARIABLES
    notGetPos = True
    notGetNeg = True
    notGetOut = True
    probLen1 = 40
    probLen2 = 20
    errorInProb1 = 0
    errorInProb2 = 1
    coverage = 1
    cluster = 1
    seqId = 0.97
    seqIdProbe = 0.90
    inputFasta = ""
    negativeFasta = ""
    output = ""
    # ARGUMENT PARSING
    for opt, value in pars:
        try:
            inputFasta = str(value) if opt in ('-p', '--positive') else inputFasta
            notGetPos = False if opt in ('-p', '--positive') else notGetPos
            negativeFasta = str(value) if opt in ('-n', '--negative') else negativeFasta
            notGetNeg = False if opt in ('-n', '--negative') else notGetNeg
            probLen1 = int(value) if opt == '--probe-len1' else probLen1
            probLen2 = int(value) if opt == '--probe-len2' else probLen2
            output = str(value) if opt in ('-o', '--output') else output
            notGetOut = False if opt in ('-o', '--output') else notGetOut
            seqId = float(value) if opt == '--seq-id-cluster' else seqId
            seqIdProbe = float(value) if opt == '--seq-id-probe' else seqIdProbe
            errorInProb1 = int(value) if opt == '--probe-error1' else errorInProb1
            errorInProb2 = int(value) if opt == '--probe-error2' else errorInProb2
            coverage = int(value) if opt == '-c' else coverage
            cluster = int(value) if opt == '--cluster' else cluster
            if opt in ('-h', '--help'):
                printUsage()
        except:
            printUsage()
    # PRINT OUT ALERTS AND INFO WHEN AGUMENTS NOT PROPER
    if notGetPos or notGetNeg or notGetOut:
        printUsage()
    # MAKE WORKING DIRECTORIES
    # absInputFasta
    absInputFasta = makeAbsPath(inputFasta)
    absNegativeFasta = makeAbsPath(negativeFasta)
    # MAKE WORKING DIRECTORY
    workDir = makeDir(output)
    # MOVE TO WORKING DIRECTORY
    os.chdir(workDir)
    # MAKE DIRECTORIES
    inputDir = makeDir('input')
    clusterDir = makeDir('cluster')
    mmseqsDir = makeDir('mmseqs')
    filterDir = makeDir('filter')
    mappingProbe1Dir = makeDir('mapping_probe1')
    setcoverProbe1Dir = makeDir('setcover_probe1')
    # COPY INPUT FASTA FILE
    clonedInputFasta = copyInputFasta(directory=inputDir, inputFasta=absInputFasta, output='genomes.fa')
    # REDUCE REDUNDANCY OF INPUT FASTA FILE
    redundancyReducedFasta, nameLookupFile1 = (
        reduceRedundancy(
            doCluster=(cluster == 1),
            seqIdentity=seqId,
            directory=clusterDir,
            inputFasta=clonedInputFasta,
            output='genomes.clu',
            lookup='id.lookup',        # name.lookup.1
        )
    )
    # MAKE INITIAL PROBES FASTA AND MASK PROBES FASTA
    #                   mmseqs.search  #sesrchAgainstNegativeFasta
    initialProbesFasta, maskedProbes = (
        makeInitialProbesMaskedProbes(
            directory=mmseqsDir,
            inputFasta=redundancyReducedFasta,
            negative=absNegativeFasta,
            output='probes.fa',
            maskedOutput='mmseqs.search',
            probeLen=probLen1,
            seqIdProbe=seqIdProbe,
        )
    )
    # WRITE MASK RESULT BED FILE
    redundancyReducedBed, crosstaxaBed = (
        makeRangeFileForFasta(
            directory=filterDir,
            inputFasta=redundancyReducedFasta,
            maskResult=maskedProbes,
            maskOutput='crosstaxa.bed',
            problen=probLen1,
        )
    )
    # MAKE A PROBES FASTA FILE FOR THERMO FILTER
    thermoFilterInputFasta = (
        makeFastaForThermoFilter(directory=filterDir, inputFasta=initialProbesFasta, output='probes.fa')
    )
    # TO FILTER PROBES by THERMODYNAMIC FEATURES
    positiveProbesBed, negativeProbesBed = (
        filterProbesByThermoProperty(
            directory=filterDir,
            output='primer3',
            inputBed=redundancyReducedBed,
            lig_probes_in=thermoFilterInputFasta,
            save_cut_probes=True,
            probe_size=probLen1,
        )
    )
    thermoCrosstaxaBed = (
        makeSubtractBed(
            directory=filterDir,
            inputBed=crosstaxaBed,
            negativeBed=negativeProbesBed,
            output='crosstaxa.{}.bed'.format('primer3'),      # thermo
        )
    )
    # function for combile negative and crosstaxa
    # TO COMPUTE MAPPABILITY
    indexFastaForComputingMappability = (
        makeIndexFastaForComputingMappability(
            fastaDirectory=clusterDir,
            computingDirectory=mappingProbe1Dir,
            inputFasta=redundancyReducedFasta,
            output='63.fasta',
        )
    )
    uniqeGenmapCSV1 = (
        computeMappability(
            directory=mappingProbe1Dir,
            inputFasta=indexFastaForComputingMappability,
            indexFile='index_probe1',
            error=errorInProb1,
            kmer=probLen1,
            selector=thermoCrosstaxaBed,
        )
    )
    # MINIMIZE 1ST PROBE SET
    minimizedProbeSetResult1, minimizedProbeSetBed1 = (
        minimzeProbeSet(
            directory=setcoverProbe1Dir,
            inputCSV=uniqeGenmapCSV1,
            inputFasta=redundancyReducedFasta,
            inputLookup=nameLookupFile1,
            outputResult='result',
            outputBed='result.bed',
            coverage=coverage,
            length=probLen1,
            proportion=0.9,
            distance=11,
            iteration=1,
            probLen=probLen1,
        )
    )
    # FINAL OUTPUT FASTA 1
    # if not os.path.exists(makeAbsPath(fileName='prob40.fa')):
    # probe1.fa
    makeFinalProbeFasta(
        inputBed=minimizedProbeSetBed1, inputFasta=redundancyReducedFasta, output=makeAbsPath(fileName='prob40.fa')
    )
    # TODO: DEBUGGING CODE DEBUG(finalProbeFasta1)
    # PROBE2
    inputProbe2Dir = makeDir('input_probe2')
    mappingProbe2Dir = makeDir('mapping_probe2')
    setcoverProbe2Dir = makeDir('setcover_probe2')
    if probLen2 != -1:  # probLen2
        # TO MAKE A NEW FASTA FILE input_probe2/seq.fa WHICH CONTAINS PROBE2 CANDIDATES
        probe2InputFasta, nameLookupFile2 = (
            makeProbe2InputFasta(
                problen=probLen1,
                window=200,
                directory=inputProbe2Dir,
                inputFasta=redundancyReducedFasta,
                inputSetcoverResult=minimizedProbeSetResult1, output='seq.fa',
                lookup='id_probe2.lookup',       # name.lookup.2
            )
        )
        # TO COMPUTE MAPPABILITY
        uniqeGenmapCSV2 = (
            computeMappability(
                directory=mappingProbe2Dir,
                inputFasta=probe2InputFasta,
                indexFile='index_probe2',
                error=errorInProb2,
                kmer=probLen2,
            )
        )
        # MINIMIZED 2ND PROBE SET
        minimizedProbeSetResult2, minimizedProbeSetBed2 = (
            minimzeProbeSet(
                directory=setcoverProbe2Dir,
                inputCSV=uniqeGenmapCSV2,
                inputFasta=probe2InputFasta,
                inputLookup=nameLookupFile2,
                outputResult='result',
                outputBed='result.bed',
                coverage=1,
                length=1,
                proportion=0.99,
                distance=20,
                iteration=10,
                probLen=probLen2,
            )
        )
        # FINAL OUTPUT FASTA 2
        makeFinalProbeFasta(
            inputBed=minimizedProbeSetBed2, inputFasta=probe2InputFasta, output=makeAbsPath(fileName='prob2.fa')
        )
    # TODO: DEBUGGING CODE DEBUG(finalProbeFasta2)


if __name__ == '__main__':
    # ARGUMENTS PROCESSING
    args = sys.argv[1:]
    optlist, args = getopt.getopt(
        args,
        'hp:n:o:c:',
        [
            'positive=',
            'negative=',
            'probe-len1=',
            'probe-len2=',
            'output=',
            'seq-id-cluster=',
            'seq-id-probe=',
            'probe-error1=',
            'probe-error2=',
            'help',
        ],
    )
    # TO CALL MAIN METHOD
    main(optlist)
