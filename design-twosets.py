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
# rpetitive sequences
problematicSeqs = ['A' * 5, 'T' * 5, 'C' * 5, 'G' * 5]
repetitiveSeqs = [Seq(each_seq) for each_seq in problematicSeqs]

# names of directories
dirForInputFasta1 = 'input'
dirForReducingRedundancy = 'cluster'
dirForMakingAndMaskingProbes = 'mmseqs'
dirForFilteringProbesWithThermoProperty = 'filter'
dirForMappingProbe1 = 'mapping_probe1'
dirForMappingIndex1 = 'index_probe1'
dirForSetcover1 = 'setcover_probe1'
dirForInputFasta2 = 'input_probe2'
dirForMappingProbe2 = 'mapping_probe2'
dirForMappingIndex2 = 'index_probe2'
dirForSetcover2 = 'setcover_probe2'

# names of files
nameOfLogFile = 'log.txt'
nameOfCopiedInputFasta = 'genomes.fa'
nameOfRedundancyReduced = 'genomes.clu'
nameOfLookup1 = 'id.lookup'  # name.lookup.1
nameOfProbes = 'probes.fa'
nameOfNegativeProbes = 'mmseqs.search'
nameOfNegativeRemovedBed = 'crosstaxa.bed'
nameOfThermoProper = 'primer3'
nameOfNegativeRemovedThermoProperBed = 'crosstaxa.primer3.bed'
nameOfDeduplicatedGenmapCSV = 'no.double.entry.csv'
nameOfUniqueGenmapCSV = 'uniq.genmap.csv'
nameOfSetcoverResult = 'result'
nameOf2ndInputFasta = 'seq.fa'
nameOfLookup2 = 'id_probe2.lookup'  # name.lookup.2
nameOfFinalProbe1 = 'prob40.fa'  # probe1.fa
nameOfFinalProbe2 = 'prob2.fa'


# method for usage message printing
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


# method for printing out and recording messages
def printAndLog(message):
    if message == '':
        return
    print(message.strip())
    with open(nameOfLogFile, 'a')as log:
        log.write(message.strip() + '\n')


# method for executing command and recording messages
def cmdAndLog(commandList, output=''):
    sp = subprocess.Popen(commandList, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = sp.communicate()
    printAndLog(stdout.decode('UTF-8'))
    printAndLog(stderr.decode('UTF-8'))
    if output != '':
        with open(output, 'a') as out:
            out.write(stdout.decode('UTF-8').strip() + '\n')
            out.write(stderr.decode('UTF-8').strip() + '\n')


# method for getting absolute paths
def makeAbsPath(fileName, directory=''):
    return os.path.abspath(directory+fileName)


# method for making directories
def makeDir(directory):
    if not os.path.exists(directory):
        os.makedirs(directory)
    return directory+os.path.sep


# method for copying files
def copyFile(inputFile, output):
    if not os.path.exists(output):
        shutil.copy(inputFile, output)
        printAndLog(inputFile + ' is copied to ' + output)


# method for making lookup files(index files)
def makeLookup(inputFasta, output):
    if os.path.exists(output):
        return output
    with open(inputFasta) as f:
        headers = [title.split()[0].strip() for title, seq in SimpleFastaParser(f)]
        lookupLines = [headers[i] + '\t' + str(i) + '\n' for i in range(len(headers))]
    with open(output, 'w') as w:
        w.writelines(lookupLines)


# method for setcover result bed files
def makeResultBed(lookup, setcoverResult, output, problen):
    if os.path.exists(output):
        return output
    genomeAndIdx = dict()
    with open(lookup) as f1:
        for line in f1:
            genome = line.split()[0].strip()
            idx = line.split()[1].strip()
            genomeAndIdx[idx] = genome
    with open(setcoverResult) as f2:
        with open(output, 'w') as w:
            for line in f2:
                mstchedKmers = line.split(';')[1].strip()
                genome = line.split(';')[0].split(',')[0]
                pos = line.split(';')[0].split(',')[1]
                w.write('\t'.join([genomeAndIdx[genome], pos, str(int(pos) + problen), mstchedKmers])+'\n')


# method for reducing redundancy using mmseqs easy-linclust
def reduceRedundancy(needCluster, seqIdentity, directory, inputFasta):
    outputFileName = makeAbsPath(directory=directory, fileName=nameOfRedundancyReduced)
    outputFasta = outputFileName + '_rep_seq.fasta'
    absLookup = makeAbsPath(fileName=nameOfLookup1)
    command = (
        'mmseqs easy-linclust {} {} cluster/tmp -v 3 --kmer-per-seq-scale 0.5 --kmer-per-seq 1000 --min-seq-id {}' +
        ' --cov-mode 1 -c 0.95 --remove-tmp-files 0'
    )
    if needCluster:
        printAndLog("cluster input sequences")
        if not os.path.exists(outputFasta):
            cmdAndLog(command.format(inputFasta, outputFileName, seqIdentity).split())
            os.system('seqkit sort {} > {} -w 0'.format(outputFasta, 'cluster/sorted.fasta'))
            os.system('seqkit rename {} > {} -w 0'.format('cluster/sorted.fasta', outputFasta))
        else:
            copyFile(inputFasta, outputFasta)
    makeLookup(outputFasta, absLookup)
    return outputFasta, absLookup


# method for making probes while split sequences into probeLen-mer(default 40mer)
def makeInitialProbesMaskedProbes(directory, inputFasta, negative, probeLen, seqIdProbe):
    absOutput = makeAbsPath(directory=directory, fileName=nameOfProbes)
    printAndLog("remove probes that aligns to negative set")
    f = open(inputFasta)
    w = open(absOutput, 'w')
    for title, seq in SimpleFastaParser(f):
        header = '>' + title.split()[0].strip()
        genomeLen = len(seq.strip())
        lines = [header + '_' + str(i + 1) + '\n' + seq[i:i + probeLen] + '\n' for i in range(genomeLen - probeLen)]
        w.writelines(lines)
    f.close()
    w.close()
    # for i in range(genomeLen - probeLen):
    #     w.write(header + '_' + str(i + 1) + '\n')
    #     w.write(seq[i:i + probeLen] + '\n')
    absMaskedOutput = makeAbsPath(directory=directory, fileName=nameOfNegativeProbes)
    if not os.path.exists(absMaskedOutput):
        command = (
            'mmseqs easy-search {} {} {} mmseqs/tmp -v 3 --spaced-kmer-mode 0 -k 13 --mask 0 -c 0.9 --min-seq-id {}' +
            ' --cov-mode 2 --alignment-mode 4 --search-type 3 ' +
            '--format-output query,target,pident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue'
        )
        cmdAndLog(command.format(absOutput, negative, absMaskedOutput, seqIdProbe).split())
    return absOutput, absMaskedOutput


# method for making the input fasta file into a bed file
def makeRangeFileFromFasta(directory, inputFasta, negativeProbes, problen):
    absOutput = makeAbsPath(directory=directory, fileName=nameOfRedundancyReduced+'_rep_seq.bed')
    absNegativeTempFile = makeAbsPath(directory=directory, fileName='mmseqs.txt')
    with open(inputFasta) as f:
        with open(absOutput, 'w') as w:
            for title, seq in SimpleFastaParser(f):
                header = title.split()[0].strip()
                w.write(header + '\t1\t' + str(len(seq.strip())) + '\n')
    with open(negativeProbes) as f:
        with open(absNegativeTempFile, 'w')as w:
            for line in f:
                c1 = line.split()[0].strip()
                b = c1.split('_')[-1]
                c1 = c1.split('_')[0]
                w.write(c1 + '\t' + str(int(b) - 1) + '\t' + str(int(b) - 1 + problen) + '\n')
    absNegativeRemovedOutput = (
        makeSubtractBed(
            directory=directory,
            inputBed=absOutput,
            negativeBed=absNegativeTempFile,
            output=nameOfNegativeRemovedBed
        )
    )
    return absOutput, absNegativeRemovedOutput


# method for making a fasta file for thermodynamic property filtering
def makeFastaForThermoFilter(directory, inputFasta):
    absOutput = makeAbsPath(directory=directory, fileName=nameOfProbes)
    with open(inputFasta) as f:
        with open(absOutput, 'w') as w:
            for title, seq in SimpleFastaParser(f):
                header = ('>' + title).split('_')
                w.write('_'.join(header[:-1]) + '\t' + header[-1] + '\n')
                w.write(seq + '\n')
    return absOutput


#
def makeIndexFastaForComputingMappability(fastaDirectory, computingDirectory, inputFasta):
    absDone = makeAbsPath(directory=computingDirectory, fileName='done')
    absOutput = makeAbsPath(directory=fastaDirectory, fileName='63.fasta')
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
def computeMappability(directory, inputFasta, indexDir, error, kmer, selector=''):
    command1 = 'genmap index -F {} -I {} > /dev/null'
    command2 = 'genmap map --no-reverse-complement -E {} {} --csv -K {} -t -b --frequency-large -I {} -O {} > /dev/null'
    selectorOpt = '' if selector == '' else '-S '+selector
    os.system(command1.format(inputFasta, indexDir))
    os.system(command2.format(error, selectorOpt, kmer, indexDir, directory))
    # Deduplicated Genmap CSV
    absDeduplicatedGenmapCSV = makeAbsPath(directory=directory, fileName=nameOfDeduplicatedGenmapCSV)
    if not os.path.exists(absDeduplicatedGenmapCSV):
        w = open(absDeduplicatedGenmapCSV, 'w')
        originalCSVs = []
        for file in os.listdir(directory):
            if file[-11:] != '.genmap.csv':
                continue
            originalCSVs.append(pd.read_csv(makeAbsPath(fileName=file, directory=directory), sep=';'))
        originalCSV = np.array(pd.concat(originalCSVs))
        for row in originalCSV:
            indexKmer = row[0]
            matchedKmers = row[1].split('|')
            foundGenomes = []
            line = ''
            for kmer in matchedKmers:
                genome = kmer.split(',')[0]
                if genome not in foundGenomes:
                    line = line + '|' + kmer
                    foundGenomes.append(genome)
            w.write(indexKmer + ';' + line[1:] + '\n')
        w.close()
    # Unique Genmap CSV
    absUniqueGenmapCSV = makeAbsPath(directory=directory, fileName=nameOfUniqueGenmapCSV)
    if os.path.exists(absUniqueGenmapCSV):
        return absUniqueGenmapCSV
    with open(absDeduplicatedGenmapCSV) as f:
        with open(absUniqueGenmapCSV, 'w') as w:
            foundKmers = []
            for line in f:
                indexKmer = line.split(';')[0]
                matchedKmers = line.split(';')[1].split('|')
                foundKmers.append(indexKmer)
                isPrevKmer = False
                for kmer in matchedKmers:
                    if kmer.strip() != indexKmer and kmer.strip() in foundKmers:
                        isPrevKmer = True
                        break
                if not isPrevKmer:
                    w.write(line)
    return absUniqueGenmapCSV


# method for executing setcover
def setcover(inputCSV, inputFasta, output, coverage, length, proportion, distance, iteration):
    command = '../setcover/setcover -c {} -l {} -p {} -d {} -i {} {} {}'
    cmdAndLog(
        command.format(coverage, length, proportion, distance, iteration, inputCSV, inputFasta).split(), output=output
    )


# make output probes fasta using seqkt subseq
def makeFinalProbeFasta(inputBed, inputFasta, output):
    absOutput = makeAbsPath(fileName=output)
    os.system('seqkit subseq --quiet --bed "{}" "{}" > "{}"'.format(inputBed, inputFasta, output))
    return absOutput


# method for input_probe2/seq.fa file
def makeProbe2InputFasta(problen, directory, inputFasta, inputSetcoverResult, lookup, output, window=200):
    absOutput = makeAbsPath(directory=directory, fileName=output)
    absLookup = makeAbsPath(fileName=lookup)
    #
    with open(inputFasta) as f1:
        seqName = dict()
        seqs = dict()
        seqIdx = 0
        for title, sequence in SimpleFastaParser(f1):
            header = ('>' + title).split()
            seqName[seqIdx] = header[0]
            seqs[seqIdx] = sequence
            seqIdx += 1
    #
    with open(inputSetcoverResult) as f2:
        probPos = dict()
        probIdx = 1
        matchedKmers = sum([line.split(';')[1].split('|') for line in f2.readlines()], [])
        for kmer in matchedKmers:
            genome = int(kmer.split(',')[0])
            probeStart = int(kmer.split(',')[1]) + 1
            probeEnd = probeStart + problen
            probPos[probIdx] = [genome, probeStart, probeEnd]
            probIdx += 1
    #
    with open(absOutput, 'w') as w:
        for key in probPos:
            seqId = probPos[key][0]
            probeStart = probPos[key][1]
            probeEnd = probPos[key][2]
            inputSequence = seqs[seqId].strip()
            seqLen = len(inputSequence)
            start = probeStart - window if probeStart > window else 1
            end = probeStart + problen + window if (probeStart + problen + window < seqLen) else seqLen
            maskMap = (
                sum([list(range(probPos[k][1]-1, probPos[k][2]-1)) for k in probPos if probPos[k][0] == seqId], [])
            )
            outputSequence = ''.join(['N' if pos in maskMap else inputSequence[pos] for pos in range(start-1, end)])
            outputHeader = ':'.join([seqName[seqId], str(probeStart - 1), str(probeEnd - 1)])
            outputHeader = '\t'.join([outputHeader, str(start), str(end), str(seqLen)])
            w.write(outputHeader + '\n')
            w.write(outputSequence + '\n')
    makeLookup(absOutput, absLookup)
    return absOutput, absLookup


#
def minimzeProbeSet(
    directory,
    inputCSV,
    inputFasta,
    inputLookup,
    coverage,
    length,
    proportion,
    distance,
    iteration,
    probLen,
):
    absOutputBed = makeAbsPath(directory=directory, fileName=nameOfSetcoverResult+'.bed')
    absOutputResult = makeAbsPath(directory=directory, fileName=nameOfSetcoverResult)
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
        w.writelines([line for line in resultList if ';' in line])
    makeResultBed(
        lookup=inputLookup,
        setcoverResult=absOutputResult,
        output=absOutputBed,
        problen=probLen,
    )
    return absOutputResult, absOutputBed


# method for calculating gc ratio
def getContentGC(oligo):
    gcs = oligo.count('G') + oligo.count('C')
    return gcs / len(oligo)


# method for
def hasLowComplexity(candidate_oligo):
    for s in problematicSeqs:
        if s in candidate_oligo:
            return True
    return False


#
def filterProbesByThermoProperty(
        directory, inputBed, ligInput="", capInput="",
        needSaveUnfiltered=True,
        # cutoff values
        minGC=0.30, maxGC=0.70, homoDimerTmCutoff=60, hairpinTmMax=60, minProbeTm=40,
        # probe size
        probeLen=40
):
    absOutput = makeAbsPath(directory=directory, fileName=nameOfThermoProper)
    absThermoProperOutput = absOutput + '_bed.bed'  # primer3_bed.bed
    absThermoImproperOutput = absOutput + '.neg.bed'  # primer3.neg.bed

    if os.path.exists(absThermoImproperOutput):
        return absThermoProperOutput, absThermoImproperOutput
    printAndLog("filter probes based on primer3")

    # PRINTOUT CUT-OFF VALUES
    printAndLog("Minimum Tm: " + str(minProbeTm))
    printAndLog("Minimum GC percentage: " + str(minGC))
    printAndLog("Maximum GC percentage: " + str(maxGC))
    printAndLog("Homodimer maximum Tm: " + str(homoDimerTmCutoff))
    printAndLog("Hairpin maximum Tm: " + str(hairpinTmMax))
    # TO MAKE DF FOR CAPTURE PROBES
    cp = []
    if capInput != "":
        cap_probes_out = absOutput + "_cap.tsv"
        with open(capInput) as f:
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
    with open(ligInput) as f:
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
            this_posEnd = this_posStart + probeLen
            posStart.append(this_posStart)
            posEnd.append(this_posEnd)
            seqs.append(sequence)
            this_seq = str(reverse_complement(Seq(sequence)))
            rc.append(this_seq)
            mid_pos = round(probeLen / 2)
            p1.append(this_seq[0:mid_pos])
            p2.append(this_seq[mid_pos:probeLen])
    ligProbes = pd.DataFrame(
        list(zip(identifiers, posStart, posEnd, seqs, rc, p1, p2)),
        columns=['id', 'chromStart', 'chromEnd', 'genome_segment', 'rc', 'p1', 'p2'],
    )
    printAndLog(str(len(ligProbes)) + " ligation probe sets inputted")
    # TO MAKE DF uniqueProbs CONTAINS SEQUENCES AND THERMODYANMIC DATA OF p1, p2 AND cp WITHOUT REDUNDAUCY
    unique_probe_set = set(p1 + p2 + cp)
    uniqueProbs = pd.DataFrame(list(unique_probe_set), columns=['p'])
    printAndLog("There were " + str(len(uniqueProbs)) + " unique probes")
    uniqueProbs = (
        uniqueProbs
        .assign(ultimate_base=uniqueProbs['p'].str[-1])
        .assign(penultimate_base=uniqueProbs['p'].str[-2])
    )
    printAndLog("Ultimate and penultimate bases assigned")
    uniqueProbs = (uniqueProbs.assign(tm=np.vectorize(primer3.calcTm)(uniqueProbs['p'])))
    uniqueProbs['hairpin_tm'] = list(map(primer3.calcHairpinTm, uniqueProbs['p']))
    uniqueProbs = (uniqueProbs.assign(homodimer_tm=np.vectorize(primer3.calcHomodimerTm)(uniqueProbs['p'])))
    uniqueProbs = (uniqueProbs.assign(intrinsic_probs=np.vectorize(hasLowComplexity)(uniqueProbs['p'])))
    printAndLog("Thermodynamic calculations complete")
    uniqueProbs = uniqueProbs.assign(GC_perc=np.vectorize(getContentGC)(uniqueProbs['p']))
    printAndLog("GC perc assigned")
    # TO MAKE ligProbsAndThermos
    ligProbsAndThermos = pd.merge(ligProbes, uniqueProbs.add_prefix('p1_'), how='left', left_on='p1', right_on='p1_p')
    ligProbsAndThermos = (
        pd.merge(ligProbsAndThermos, uniqueProbs.add_prefix('p2_'), how='left', left_on='p2', right_on='p2_p')
    )
    # VALUES FOR SORTING
    sorder = ['p1_hairpin_tm', 'p2_hairpin_tm', 'p1_homodimer_tm', 'p2_homodimer_tm']
    ascending = [True, True, True, True]
    # RETURN STRING
    outputStatement = "Ligation probes output: " + absOutput + '.tsv'  # primer3.tsv
    # TO FILTER CAPTURE PROBES USING THERMODYNAMIC FEATURES AND SAVE FILTERED CAPTURE PROBES AS TSV FILE
    if capInput != "":
        # TO MAKE cap_probes_calc
        cap_probes_calc = (
            pd.merge(cap_probes, uniqueProbs.add_prefix('cp_'), how='left', left_on='cp', right_on='cp_p')
        )
        # FILTER
        cap_probes_filtered = (
            cap_probes_calc
            .query('cp_intrinsic_probs == False')
            .query('cp_tm >= ' + str(minProbeTm))
            .query('cp_GC_perc >= ' + str(minGC) + ' & cp_GC_perc <= ' + str(maxGC))
            .query('cp_homodimer_tm <= ' + str(homoDimerTmCutoff))
            .query('cp_hairpin_tm <= ' + str(hairpinTmMax))
        )
        # SORT AND SAVE
        cap_probes_sorted = (
            cap_probes_filtered.sort_values(by=['cp_hairpin_tm', 'cp_homodimer_tm'], ascending=[True, True])
        )
        cap_probes_sorted.to_csv(cap_probes_out, sep='\t')
        printAndLog('Capture probes that passed all filters: ' + str(len(cap_probes_filtered)))
        outputStatement = outputStatement + "\nCapture probes output: " + cap_probes_out
    # TO FILTER LIGATION PROBES USING THERMODYNAMIC FEATURES AND SAVE FILTERED LIGATION PROBES AS TSV FILE
    # FILTER
    ligProbsFiltered = (
        ligProbsAndThermos
        .query('p1_intrinsic_probs == False & p2_intrinsic_probs == False')
        .query('p1_tm >= ' + str(minProbeTm) + ' & p2_tm >= ' + str(minProbeTm))
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
        .query('p1_homodimer_tm <= ' + str(homoDimerTmCutoff) + ' & p2_homodimer_tm <= ' + str(homoDimerTmCutoff))
        .query('p1_hairpin_tm <= ' + str(hairpinTmMax) + ' & p2_hairpin_tm <= ' + str(hairpinTmMax))
    )
    # SORT AND SAVE
    printAndLog('Ligation probe sets that passed all filters: ' + str(len(ligProbsFiltered)))
    ligProbsFilteredAndSorted = ligProbsFiltered.sort_values(by=sorder, ascending=ascending)
    ligProbsFiltered.rename(
        columns={'id': 'chrom'}
    )[['chrom', 'chromStart', 'chromEnd']].to_csv(absThermoProperOutput, index=False, sep='\t')
    ligProbsFilteredAndSorted.to_csv((absOutput + '_lig.tsv'), sep='\t')  # primer3_lig.tsv
    # OPTIONAL: save dataframe ligProbsAndThermos as tsv file
    if needSaveUnfiltered:
        ligProbsAndThermos.to_csv((absOutput + '_uncut.tsv'), sep='\t')  # primer3_uncut.tsv
    printAndLog(outputStatement)
    absThermoImproperOutput = (
        makeSubtractBed(
            directory=directory,
            inputBed=inputBed,
            negativeBed=absThermoProperOutput,
            output=nameOfThermoProper + '.neg.bed',  # primer3.neg.bed
        )
    )
    return absThermoProperOutput, absThermoImproperOutput


#
def makeSubtractBed(directory, inputBed, negativeBed, output):
    absOutput = makeAbsPath(directory=directory, fileName=output)
    os.system('bedtools subtract -a {} -b {} > {}'.format(inputBed, negativeBed, absOutput))
    return absOutput


#
def countErrors(seq1, seq2):
    if len(seq1) != len(seq2):
        return -1
    cnt = 0
    for i in range(len(seq1)):
        if seq1[i] != seq2[i]:
            cnt += 1
    return cnt


def confirmFinalProbFasta(inputFasta, finalProbFasta, probLen, error):
    f1 = open(finalProbFasta)
    for t1, s1 in SimpleFastaParser(f1):
        cnt1 = t1.count('|') + 1
        cnt2 = 0
        with open(inputFasta) as f2:
            foundTitles = []
            for t2, s2 in SimpleFastaParser(f2):
                for i in range(len(s2) - probLen + 1):
                    if countErrors(s1, s2[i:i + probLen]) in range(error + 1):
                        cnt2 += 1
                        foundTitles.append(t2.split()[0])
                        break
            if cnt1 != cnt2:
                printAndLog('An Error was found in {}.'.format(finalProbFasta))
                printAndLog(t1)
                printAndLog('Seq {} was expected to be found {} times.'.format(s1, cnt1))
                printAndLog('However it was found {} times'.format(cnt2))
                printAndLog(', '.join(foundTitles))
                quit()
    f1.close()
    printAndLog('Probe file {} has no error.'.format(finalProbFasta))


# MAIN METHOD
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
    cluster = True
    seqId = 0.97
    seqIdProbe = 0.90
    inputFasta = ""
    negativeFasta = ""
    output = ""
    # ARGUMENT PARSING
    for opt, value in pars:
        if opt in ('-h', '--help'):
            printUsage()
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
            cluster = int(value) == 1 if opt == '--cluster' else cluster
        except Exception:
            printUsage()
    # PRINT OUT ALERTS AND INFO WHEN AGUMENTS NOT PROPER
    if notGetPos or notGetNeg or notGetOut:
        printUsage()
    # get global variable
    # MAKE WORKING DIRECTORIES
    # absInputFasta
    absInputFasta = makeAbsPath(inputFasta)
    absNegativeFasta = makeAbsPath(negativeFasta)
    # MAKE WORKING DIRECTORY
    workDir = makeDir(output)
    # MOVE TO WORKING DIRECTORY
    os.chdir(workDir)
    printAndLog(' '.join(['{} {}'.format(i[0], i[1]).strip() for i in optlist]))
    # MAKE DIRECTORIES
    inputDir = makeDir(dirForInputFasta1)
    deduplicatingDir = makeDir(dirForReducingRedundancy)
    maskingDir = makeDir(dirForMakingAndMaskingProbes)
    thermoFilteringDir = makeDir(dirForFilteringProbesWithThermoProperty)
    mappingDir1 = makeDir(dirForMappingProbe1)
    minimizingDir1 = makeDir(dirForSetcover1)
    # COPY INPUT FASTA FILE
    clonedInputFasta = makeAbsPath(directory=inputDir, fileName=nameOfCopiedInputFasta)
    copyFile(inputFile=absInputFasta, output=clonedInputFasta)
    # REDUCE REDUNDANCY OF INPUT FASTA FILE
    redundancyReducedFasta, nameLookupFile1 = (
        reduceRedundancy(
            needCluster=cluster,
            seqIdentity=seqId,
            directory=deduplicatingDir,
            inputFasta=clonedInputFasta,
        )
    )
    # MAKE INITIAL PROBES FASTA AND MASK PROBES FASTA
    #                   mmseqs.search  #sesrchAgainstNegativeFasta
    initProbesFasta, negativeProbes = (
        makeInitialProbesMaskedProbes(
            directory=maskingDir,
            inputFasta=redundancyReducedFasta,
            negative=absNegativeFasta,
            probeLen=probLen1,
            seqIdProbe=seqIdProbe,
        )
    )
    # WRITE MASK RESULT BED FILE
    redundancyReducedBed, negativeRemovedBed = (
        makeRangeFileFromFasta(
            directory=thermoFilteringDir,
            inputFasta=redundancyReducedFasta,
            negativeProbes=negativeProbes,
            problen=probLen1,
        )
    )
    # MAKE A PROBES FASTA FILE FOR THERMO FILTER
    thermoFilterInputFasta = makeFastaForThermoFilter(directory=thermoFilteringDir, inputFasta=initProbesFasta)
    # TO FILTER PROBES by THERMODYNAMIC FEATURES
    thermoProperProbesBed, thermoImproperProbesBed = (
        filterProbesByThermoProperty(
            directory=thermoFilteringDir,
            inputBed=redundancyReducedBed,
            ligInput=thermoFilterInputFasta,
            probeLen=probLen1,
        )
    )
    negativeRemovedThermoProperBed = (
        makeSubtractBed(
            directory=thermoFilteringDir,
            inputBed=negativeRemovedBed,
            negativeBed=thermoImproperProbesBed,
            output=nameOfNegativeRemovedThermoProperBed
        )
    )
    # function for combile negative and crosstaxa
    # TO COMPUTE MAPPABILITY
    indexFastaForComputingMappability = (
        makeIndexFastaForComputingMappability(
            fastaDirectory=deduplicatingDir, computingDirectory=mappingDir1, inputFasta=redundancyReducedFasta
        )
    )
    uniqeGenmapCSV1 = (
        computeMappability(
            directory=mappingDir1,
            inputFasta=indexFastaForComputingMappability,
            indexDir=dirForMappingIndex1,
            error=errorInProb1,
            kmer=probLen1,
            selector=negativeRemovedThermoProperBed,
        )
    )
    # MINIMIZE 1ST PROBE SET
    minimizedProbeSetResult1, minimizedProbeSetBed1 = (
        minimzeProbeSet(
            directory=minimizingDir1,
            inputCSV=uniqeGenmapCSV1,
            inputFasta=redundancyReducedFasta,
            inputLookup=nameLookupFile1,
            coverage=coverage,
            length=probLen1,
            proportion=0.9,  #
            distance=11,  #
            iteration=1,  #
            probLen=probLen1,
        )
    )
    # FINAL OUTPUT FASTA 1
    # if not os.path.exists(makeAbsPath(fileName='prob40.fa')):
    # probe1.fa
    finalProbeFasta1 = makeFinalProbeFasta(
        inputBed=minimizedProbeSetBed1,
        inputFasta=redundancyReducedFasta,
        output=makeAbsPath(fileName=nameOfFinalProbe1),
    )
    confirmFinalProbFasta(redundancyReducedFasta, finalProbeFasta1, probLen1, errorInProb1)
    # PROBE2
    inputDir2 = makeDir(dirForInputFasta2)
    mappingDir2 = makeDir(dirForMappingProbe2)
    minimizingDir2 = makeDir(dirForSetcover2)
    if probLen2 == -1:  # probLen2
        return
    # TO MAKE A NEW FASTA FILE input_probe2/seq.fa WHICH CONTAINS PROBE2 CANDIDATES
    probe2InputFasta, nameLookupFile2 = (
        makeProbe2InputFasta(
            problen=probLen1,
            directory=inputDir2,
            inputFasta=redundancyReducedFasta,
            inputSetcoverResult=minimizedProbeSetResult1,
            lookup=nameOfLookup2,
            output=nameOf2ndInputFasta,
        )
    )
    # TO COMPUTE MAPPABILITY
    uniqeGenmapCSV2 = (
        computeMappability(
            directory=mappingDir2,
            inputFasta=probe2InputFasta,
            indexDir=dirForMappingIndex2,
            error=errorInProb2,
            kmer=probLen2,
        )
    )
    # MINIMIZED 2ND PROBE SET
    minimizedProbeSetResult2, minimizedProbeSetBed2 = (
        minimzeProbeSet(
            directory=minimizingDir2,
            inputCSV=uniqeGenmapCSV2,
            inputFasta=probe2InputFasta,
            inputLookup=nameLookupFile2,
            coverage=1,  #
            length=1,  #
            proportion=0.99,  #
            distance=20,  #
            iteration=10,  #
            probLen=probLen2,
        )
    )
    # FINAL OUTPUT FASTA 2
    finalProbeFasta2 = makeFinalProbeFasta(
        inputBed=minimizedProbeSetBed2, inputFasta=probe2InputFasta, output=makeAbsPath(fileName=nameOfFinalProbe2)
    )
    confirmFinalProbFasta(probe2InputFasta, finalProbeFasta2, probLen2, errorInProb2)


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
