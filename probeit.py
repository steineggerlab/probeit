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
import re


class ParaSeqs:
    probLen = -1
    ntSnp = ''
    aaSnp = ''
    wtSeq = ''
    stSeq = ''
    mutLoc = -1
    found = -1

    def __init__(self, ntSnp, aaSnp, wtSeq, stSeq, probLen=40, mutLoc=-1, found=-1):
        self.ntSnp = ntSnp
        self.aaSnp = aaSnp
        self.wtSeq = wtSeq
        self.stSeq = stSeq
        self.probLen = probLen
        self.mutLoc = mutLoc
        self.found = found

    def __eq__(self, other):
        isSameClass = self.__class__ == other.__class__
        isSameNtSNP = self.ntSnp == other.ntSnp
        isSameWtSeq = self.wtSeq == other.wtSeq
        isSameStSeq = self.stSeq == other.stSeq
        return isSameClass and isSameNtSNP and isSameWtSeq and isSameStSeq

    def getProbsWithPos(self, pos):
        return (
            self.wtSeq[self.mutLoc - pos + 1:self.mutLoc - pos + self.probLen + 1],
            self.stSeq[self.mutLoc - pos + 1:self.mutLoc - pos + self.probLen + 1]
        )


class ProbeitUtils:
    def runCommand(self, command, verbose=False):
        if verbose:
            commandList = command.split()
            sp = subprocess.Popen(commandList, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            stdout, stderr = sp.communicate()
            return stdout.decode('UTF-8'), stderr.decode('UTF-8')
        else:
            os.system(command)

    def delDir(self, waste):
        command = "rm -rf " + waste
        self.runCommand(command)

    # TO CALL SEQKIT MODULES
    def sortFasta(self, fasta, sortedFasta):
        command = "seqkit sort {} > {} -w 0".format(fasta, sortedFasta)
        self.runCommand(command)

    def renameFasta(self, oldFasta, newFasta):
        command = "seqkit rename {} > {} -w 0".format(oldFasta, newFasta)
        self.runCommand(command)

    def getSubseqFasta(self, coordinateBed, oldFasta, newFasta):
        command = "seqkit subseq --quiet --bed {} {} > {}".format(coordinateBed, oldFasta, newFasta)
        self.runCommand(command)

    def getPatternPosition(self, patternFasta, genomeFasta, positonsTSV):
        command = "seqkit locate -f {} {} > {}".format(patternFasta, genomeFasta, positonsTSV)
        self.runCommand(command)
        return positonsTSV

    # TO CALL GENMAP MODULES
    # COMPUTE MAPPABILITY
    def comMappability(self, fasta, indexDir, error, kmer, outputDir, selectingOption=''):
        self.delDir(indexDir)
        command0 = " > /dev/null"
        command1 = "genmap index -F {} -I {}".format(fasta, indexDir)
        self.runCommand(command1 + command0)
        command2 = (
            "genmap map --no-reverse-complement -E {} {} --csv -K {} -t -b --frequency-large -I {} -O {}"
            .format(error, selectingOption, kmer, indexDir, outputDir)
        )
        self.runCommand(command2 + command0)
        return outputDir + '.'.join(fasta.split(os.path.sep)[-1].split('.')[:-1] + ['genmap', 'csv'])

    # TO CALL BEDTOOLS MODULES
    def getSubtractedBed(self, positiveBed, negativeBed, outputBed):
        print(positiveBed, negativeBed, outputBed)
        command = "bedtools subtract -a {} -b {} > {}".format(positiveBed, negativeBed, outputBed)
        self.runCommand(command)
        return outputBed

    def getWindowFasta(self, genomeFasta, maskingBed, maskedGenomeFasta, windowBed, windowFasta):
        command1 = "bedtools maskfasta -fi {} -bed {} -fo {}".format(genomeFasta, maskingBed, maskedGenomeFasta)
        self.runCommand(command1)
        inputDF = pd.read_csv(maskingBed, sep='\t', header=None)
        inputDF[1] = inputDF[1].apply(lambda x: x - 200)
        inputDF[2] = inputDF[2].apply(lambda x: x + 200)
        inputDF.to_csv(windowBed, sep='\t', header=False, index=False)
        command2 = "bedtools getfasta -fi {} -bed {} > {}".format(maskedGenomeFasta, windowBed, windowFasta)
        self.runCommand(command2)
        return windowFasta

    def makeLookup(self, windowFasta, lookup):
        cnt = 0
        with open(windowFasta) as f:
            with open(lookup, 'w') as w:
                for h, s in SimpleFastaParser(f):
                    w.write('{}\t{}\n'.format(h, cnt))
                    cnt += 1
        return lookup

    def doBlastSearch(self, workDir, inputFasta, strGenomeFasta, resultTSV):
        searchDir = workDir + 'search' + os.path.sep
        tempDir = searchDir + 'temp' + os.path.sep
        resultTSV = workDir + resultTSV
        if not os.path.isdir(searchDir):
            os.makedirs(searchDir)
        searchdb = searchDir + 'searchDB'
        strdb = searchDir + 'strainDB'
        aln = searchDir + 'mmseqs.aln'
        command1 = 'mmseqs createdb {} {}'
        command2 = 'mmseqs createdb {} {}'
        command3 = 'mmseqs search {} {} {} {} --search-type 3'
        command4 = 'mmseqs convertalis {} {} {} {} --format-output target,query,tseq,tstart,tend --search-type 3'
        self.runCommand(command1.format(inputFasta, searchdb), verbose=True)
        self.runCommand(command2.format(strGenomeFasta, strdb), verbose=True)
        self.runCommand(command3.format(searchdb, strdb, aln, tempDir), verbose=True)
        self.runCommand(command4.format(searchdb, strdb, aln, resultTSV), verbose=True)
        df = pd.read_csv(resultTSV, sep='\t', header=None)
        df.columns=['substr', 'snp', 'strseq', 'start', 'end']
        df['aln'] = df.apply(lambda x: x[2][int(x[3]-1):int(x[4])], axis=1)
        df['len'] = df.aln.apply(lambda x: len(x)-1)
        df = df[['substr', 'snp', 'len', 'aln']]
        df.to_csv(resultTSV, header=False, index=False, sep='\t')
        return resultTSV

    def clusterGenome(self, inputFasta, outputFasta, outputDir, seqIdentity):
        command = ' '.join(
            [
                "mmseqs easy-linclust", inputFasta, outputFasta, outputDir,
                "-v 3 --kmer-per-seq-scale 0.5 --kmer-per-seq 1000 --min-seq-id", str(seqIdentity),
                "--cov-mode 1 -c 0.95 --remove-tmp-files 0"
            ]
        )
        return self.runCommand(command, verbose=True)

    def searchNegative(self, output, negative, maskOutput, outputDir, seqInProbe):
        command = ' '.join(
            [
                "mmseqs easy-search", output, negative, maskOutput, outputDir+"tmp",
                "-v 3 --spaced-kmer-mode 0 -k 13 --mask 0 -c 0.9 --min-seq-id",
                str(seqInProbe), "--cov-mode 2 --alignment-mode 4 --search-type 3",
                "--format-output query,target,pident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue"
            ]
        )
        stdOut, stdErr = self.runCommand(command, verbose=True)
        return stdOut, stdErr

    def deDuplicateProbesCSV(self, inputCSV, outputCSV):
        def deduplicate(kmers):
            prevIndex = []
            newKmers = []
            for k in kmers:
                index = k.split(',')[0]
                if index not in prevIndex:
                    prevIndex.append(index)
                    newKmers.append(k)
            return '|'.join(newKmers)

        df = pd.read_csv(inputCSV, sep=';')
        df.columns = ['repKmer', 'kmers']
        # DEDUPLICATE EACH ROW
        df['kmers'] = df.kmers.apply(lambda x: x.split('|'))
        df['kmers'] = df.kmers.apply(lambda x: deduplicate(x))
        # DEDUPLICATE WHOLE ROWS
        repKmers = list(df['repKmer'])
        kmerLists = list(df['kmers'])
        repKmers = ['|{}|'.format(i) for i in repKmers]
        kmerLists = ['|{}|'.format(i) for i in kmerLists]
        newKmerLists = []
        maskingKmers = []
        for i in range(len(kmerLists)):
            if i == 0:
                newKmerLists.append(kmerLists[i])
                continue
            else:
                kmerList = kmerLists[i]
                maskingKmers.append(repKmers[i - 1])
                for j in maskingKmers:
                    if j in kmerList:
                        newKmerLists.append(None)
                        break
                else:
                    newKmerLists.append(kmerList)

        repKmers = [i[1:-1] if i else None for i in repKmers]
        newKmerLists = [i[1:-1] if i else None for i in newKmerLists]
        uniqGenmapDF = pd.DataFrame({'k-mer': repKmers, 'k-mers': newKmerLists})
        uniqGenmapDF.dropna().to_csv(outputCSV, header=False, index=False, sep=';')
        return outputCSV

    def setCover(self, coverage, length, proportion, distance, iteration, deDuplicatedCSV, windowFasta):
        command = (
            "./setcover/setcover -c {} -l {} -p {} -d {} -i {} {} {}"
            .format(coverage, length, proportion, distance, iteration, deDuplicatedCSV, windowFasta)
        )
        stdOut, stdErr = self.runCommand(command, verbose=True)
        print(stdOut, stdErr)
        return stdOut, stdErr

    def makeSetcoverResultBed(self, lookup, setcoverResult, setcoverResultBed):
        genomeAndIdx = dict()
        with open(lookup) as f1:
            for line in f1:
                genome = line.split()[0].strip()
                idx = line.split()[1].strip()
                genomeAndIdx[idx] = genome
        with open(setcoverResult) as f2:
            with open(setcoverResultBed, 'w') as w:
                for line in f2:
                    matchedKmers = line.split(';')[1].strip()
                    genome = line.split(';')[0].split(',')[0]
                    pos = line.split(';')[0].split(',')[1]
                    w.write('\t'.join([genomeAndIdx[genome], pos, str(int(pos) + 20), matchedKmers]) + '\n')
        return setcoverResultBed

    def getUserArgs(self, args):
        return ' '.join(['{} {}'.format(i[0], i[1]) if len(i) == 2 else i for i in args])


class Probeit:
    args = []
    subWork = None

    def __init__(self, args):
        self.args = args

    def checkArgs(self):
        if self.args == [] or self.args[0] == '-h' or self.args[0] == '--help':
            self.printUsage()
            return
        elif self.args[0] == 'posnegset':
            print('CURRENT: ', os.getcwd())
            self.subWork = PosNegSet(self.args[1:])
            self.subWork.excute()
            return
        elif self.args[0] == 'snp':
            self.subWork = SNP(self.args[1:])
            self.subWork.excute()
            return
        else:
            self.printUsage()
            return

    def printUsage(self):
        print("PROBEIT")
        print('usage ./probeit.py <module> [<args>]')
        print("WORKFLOWS")
        print("posnegset: make two-sets probes with positive and negative sets")
        print("snp: make two-sets probes with wildtype genome, strain genome and SNPs")
        quit()


class PosNegSet:
    args = []
    shortParams = 'hp:n:o:c:'
    longParams = [
        'positive=', 'negative=', 'output=', 'probe-len1=', 'probe-len2=', 'seq-id-cluster=',
        'seq-id-probe=', 'probe-error1=', 'probe-error2=', 'help',
    ]
    inputGenome = ''
    negGenome = ''
    workDir = ''
    utils = ProbeitUtils()
    probLen1 = 40
    probLen2 = 20
    seqIdentity = 0.97 # for mmseqs linclust
    seqIdProbe = 0.90 # for mmseqs easy-search
    error1 = 0
    error2 = 1
    coverage = 1
    needCluster = True
    window = 200

    def __init__(self, args):
        self.args = getopt.getopt(args, self.shortParams, self.longParams)[0]

    def getPars(self):
        for opt, value in self.args:
            if opt in ('-h', '--help'):
                self.printUsage()
            try:
                self.inputGenome = str(value) if opt in ('-p', '--positive') else self.inputGenome
                self.negGenome = str(value) if opt in ('-n', '--negative') else self.negGenome
                self.probLen1 = int(value) if opt == '--probe-len1' else self.probLen1
                self.probLen2 = int(value) if opt == '--probe-len2' else self.probLen2
                self.workDir = str(value) if opt in ('-o', '--output') else self.workDir
                self.seqIdentity = float(value) if opt == '--seq-id-cluster' else self.seqIdentity
                self.seqIdProbe = float(value) if opt == '--seq-id-probe' else self.seqIdProbe
                self.error1 = int(value) if opt == '--probe-error1' else self.error1
                self.error2 = int(value) if opt == '--probe-error2' else self.error2
                self.coverage = int(value) if opt == '-c' else self.coverage
                self.needCluster = int(value) == 1 if opt == '--cluster' else self.needCluster
            except Exception:
                self.printUsage()

    def checkArgs(self):
        print("Your arguments: {}".format('snp ' + self.utils.getUserArgs(self.args)))
        message = "{}You didn't input a proper argument for '{}' parameter or missed it. "
        isBadArguments = False
        if not self.inputGenome:
            print(message.format('[ERR]', '--positive'))
            isBadArguments = True
        if not self.negGenome:
            print(message.format('[ERR]', '--negative'))
            isBadArguments = True
        if not self.workDir:
            print(message.format('[ERR]', '--output'))
            isBadArguments = True
        if isBadArguments:
            self.printUsage()
        else:
            print('You input proper arguments.')

    def makeWorkDir(self):
        self.workDir = self.workDir + os.path.sep
        self.inputDir = self.workDir + 'input' + os.path.sep
        self.dedupDir = self.workDir + 'cluster' + os.path.sep
        self.maskingDir = self.workDir + 'mmseqs' + os.path.sep
        self.thermoFilteringDir = self.workDir + 'filter' + os.path.sep
        self.mapDir1 = self.workDir + 'mapping_probe1' + os.path.sep
        self.indexDir1 = self.workDir + 'index_probe1' + os.path.sep
        self.minimizingDir1 = self.workDir + 'setcover_probe1' + os.path.sep
        self.inputDir2 = self.workDir + 'input_probe2' + os.path.sep
        self.mapDir2 = self.workDir + 'mapping_probe2' + os.path.sep
        self.indexDir2 = self.workDir + 'index_probe2' + os.path.sep
        self.minimizingDir2 = self.workDir + 'setcover_probe2' + os.path.sep
        self.log = self.workDir + 'log.txt'
        dirList = [
            self.workDir, self.inputDir, self.dedupDir, self.maskingDir,
            self.thermoFilteringDir, self.mapDir1, self.minimizingDir1,
            self.inputDir2, self.mapDir2, self.minimizingDir2
        ]
        for directory in dirList:
            if not os.path.exists(directory):
                os.makedirs(directory)
        return

    def logUpdate(self, msg):
        with open(self.log, 'a') as w:
            w.write(msg + '\n')
        print(msg)

    def copyFile(self, original, copy):
        if not os.path.exists(copy):
            shutil.copy(original, copy)
            self.logUpdate('{} is copied to {}'.format(original, copy))
        return copy

    def makeLookup(self, inputFasta, output):
        if os.path.exists(output):
            return output
        with open(inputFasta) as f:
            headers = [title.split()[0].strip() for title, seq in SimpleFastaParser(f)]
            lookupLines = [headers[i] + '\t' + str(i) + '\n' for i in range(len(headers))]
        with open(output, 'w') as w:
            w.writelines(lookupLines)
        return output

    def thrmoFilter(
            self, directory, inputBed, ligInput, capInput="", doSaveUnfiltered=True,
    ):
        probeLen = self.probLen1
        minGC = 0.30
        maxGC = 0.70
        homoDimerTmCutoff = 60
        hairpinTmMax = 60
        minProbeTm = 40
        absOutput = directory + 'primer3'
        absThermoProperOutput = absOutput + '_bed.bed'  # primer3_bed.bed
        thermoImproperOutput = absOutput + '.neg.bed'  # primer3.neg.bed
        problematicSeqs = ['A' * 5, 'T' * 5, 'C' * 5, 'G' * 5]

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

        if os.path.exists(thermoImproperOutput):
            return thermoImproperOutput
        self.logUpdate("filter probes based on primer3")
        # PRINTOUT CUT-OFF VALUES
        self.logUpdate("Minimum Tm: " + str(minProbeTm))
        self.logUpdate("Minimum GC percentage: " + str(minGC))
        self.logUpdate("Maximum GC percentage: " + str(maxGC))
        self.logUpdate("Homodimer maximum Tm: " + str(homoDimerTmCutoff))
        self.logUpdate("Hairpin maximum Tm: " + str(hairpinTmMax))
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
            self.logUpdate(str(len(cap_probes)) + " capture probes inputted")
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
        self.logUpdate(str(len(ligProbes)) + " ligation probe sets inputted")
        # TO MAKE DF uniqueProbs CONTAINS SEQUENCES AND THERMODYANMIC DATA OF p1, p2 AND cp WITHOUT REDUNDAUCY
        unique_probe_set = set(p1 + p2 + cp)
        uniqueProbs = pd.DataFrame(list(unique_probe_set), columns=['p'])
        self.logUpdate("There were " + str(len(uniqueProbs)) + " unique probes")
        uniqueProbs = (
            uniqueProbs.assign(ultimate_base=uniqueProbs['p'].str[-1]).assign(penultimate_base=uniqueProbs['p'].str[-2])
        )
        self.logUpdate("Ultimate and penultimate bases assigned")
        uniqueProbs = (uniqueProbs.assign(tm=np.vectorize(primer3.calcTm)(uniqueProbs['p'])))
        uniqueProbs['hairpin_tm'] = list(map(primer3.calcHairpinTm, uniqueProbs['p']))
        uniqueProbs = (uniqueProbs.assign(homodimer_tm=np.vectorize(primer3.calcHomodimerTm)(uniqueProbs['p'])))
        uniqueProbs = (uniqueProbs.assign(intrinsic_probs=np.vectorize(hasLowComplexity)(uniqueProbs['p'])))
        self.logUpdate("Thermodynamic calculations complete")
        uniqueProbs = uniqueProbs.assign(GC_perc=np.vectorize(getContentGC)(uniqueProbs['p']))
        self.logUpdate("GC perc assigned")
        # TO MAKE ligProbsAndThermos
        ligProbsAndThermos = pd.merge(
            ligProbes, uniqueProbs.add_prefix('p1_'), how='left', left_on='p1', right_on='p1_p'
        )
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
            self.logUpdate('Capture probes that passed all filters: ' + str(len(cap_probes_filtered)))
            outputStatement = outputStatement + "\nCapture probes output: " + cap_probes_out
        # TO FILTER LIGATION PROBES USING THERMODYNAMIC FEATURES AND SAVE FILTERED LIGATION PROBES AS TSV FILE
        ligProbsFiltered = (
            ligProbsAndThermos
            .query('p1_intrinsic_probs == False & p2_intrinsic_probs == False')
            .query('p1_tm >= {' + str(minProbeTm) + '} & p2_tm >= {}'.format(minProbeTm, minProbeTm))
            .query('p1_GC_perc >= {} & p2_GC_perc >= {}& p1_GC_perc <= {} & p2_GC_perc <= {}'.format(
                minGC, minGC, maxGC, maxGC)
            )
            .query('p1_homodimer_tm <= {} & p2_homodimer_tm <= {}'.format(homoDimerTmCutoff, homoDimerTmCutoff))
            .query('p1_hairpin_tm <= {} & p2_hairpin_tm <= {}'.format(hairpinTmMax, hairpinTmMax))
        )
        # SORT AND SAVE
        self.logUpdate('Ligation probe sets that passed all filters: ' + str(len(ligProbsFiltered)))
        ligProbsFilteredAndSorted = ligProbsFiltered.sort_values(by=sorder, ascending=ascending)
        ligProbsFiltered.rename(
            columns={'id': 'chrom'}
        )[['chrom', 'chromStart', 'chromEnd']].to_csv(absThermoProperOutput, index=False, sep='\t')
        ligProbsFilteredAndSorted.to_csv((absOutput + '_lig.tsv'), sep='\t')  # primer3_lig.tsv
        # OPTIONAL: save dataframe ligProbsAndThermos as tsv file
        if doSaveUnfiltered:
            ligProbsAndThermos.to_csv((absOutput + '_uncut.tsv'), sep='\t')  # primer3_uncut.tsv
        self.logUpdate(outputStatement)
        return self.utils.getSubtractedBed(inputBed, absThermoProperOutput, directory + 'primer3.neg.bed')

    def filterInputData(self):
        clusterResult = self.dedupDir + 'genomes.clu'
        # COPY POSITIVE INPUT FASTA
        posGenome = self.copyFile(original=self.inputGenome, copy=self.inputDir + 'genome.fa')
        print(posGenome)
        # CLUTER INPUT FASTA
        self.deDupGenome = clusterResult + '_rep_seq.fasta'
        if self.needCluster:
            self.logUpdate('[DEDUPLICATE POSITIVE FASTA]')
            if not os.path.exists(self.deDupGenome):
                self.utils.clusterGenome(posGenome, clusterResult, self.dedupDir +'temp' + os.path.sep, self.seqIdentity)
                self.utils.sortFasta(self.deDupGenome, self.dedupDir + 'sorted.fasta')
                self.utils.renameFasta(self.dedupDir + 'sorted.fasta', self.deDupGenome)
            else:
                self.copyFile(posGenome, self.deDupGenome)
        # MAKE A LOOKUP
        self.lookup1 = self.makeLookup(self.deDupGenome, self.workDir + 'id.lookup')
        # MAKE 40MER FROM POSITIVE FASTA
        posProbes = self.maskingDir + 'probes.fa'
        self.logUpdate("remove probes that aligns to negative set")
        self.logUpdate('[GET {}mer PROBES FROM POSITIVE GENOME]'.format(self.probLen1))
        f = open(self.deDupGenome)
        w = open(posProbes, 'w')
        for title, seq in SimpleFastaParser(f):
            header = '>' + title.split()[0].strip()
            genomeLen = len(seq.strip())
            lines = [header + '_' + str(i + 1) + '\n' + seq[i:i + self.probLen1] + '\n' for i in
                     range(genomeLen - self.probLen1)]
            w.writelines(lines)
        f.close()
        w.close()
        # MAKE 40MER BED FROM NEGATIVE FASTA
        self.logUpdate('[REMOVE PROBES FOUND IN NEGATIVE GENOME]')
        negProbesCoords = self.maskingDir + 'mmseqs.search'
        if not os.path.exists(negProbesCoords):
            self.utils.searchNegative(posProbes, self.negGenome, negProbesCoords, self.maskingDir, self.seqIdProbe)
        # MAKE DEDUPLICATE GENOME POSITIONS
        deDupGenomeCoords = clusterResult + '_rep_seq.bed'
        with open(self.deDupGenome) as f:
            with open(deDupGenomeCoords, 'w') as w:
                for title, seq in SimpleFastaParser(f):
                    header = title.split()[0].strip()
                    w.write(header + '\t1\t' + str(len(seq.strip())) + '\n')
        #
        simpleNegProbesCoords = self.thermoFilteringDir + 'mmseqs.txt'
        with open(negProbesCoords) as f:
            with open(simpleNegProbesCoords, 'w')as w:
                for line in f:
                    c1 = line.split()[0].strip()
                    b = c1.split('_')[-1]
                    c1 = c1.split('_')[0]
                    w.write(c1 + '\t' + str(int(b) - 1) + '\t' + str(int(b) - 1 + self.probLen1) + '\n')
        negRemCoords = self.utils.getSubtractedBed(deDupGenomeCoords, simpleNegProbesCoords, self.thermoFilteringDir + 'crosstaxa.bed')
        self.logUpdate('[FILTER PROBES WITH THERMODYNAMIC FEATURES]')
        probesForTheromFilter = self.thermoFilteringDir + 'probes.fa'
        with open(posProbes) as f:
            with open(probesForTheromFilter, 'w') as w:
                for title, seq in SimpleFastaParser(f):
                    header = ('>' + title).split('_')
                    w.write('_'.join(header[:-1]) + '\t' + header[-1] + '\n')
                    w.write(seq + '\n')
        # TO FILTER PROBES by THERMODYNAMIC FEATURES
        thermoImpCoords = self.thrmoFilter(self.thermoFilteringDir, deDupGenomeCoords, probesForTheromFilter)
        self.thermoProCoords = self.utils.getSubtractedBed(negRemCoords, thermoImpCoords, self.thermoFilteringDir + 'crosstaxa.primer3.bed')
        self.mapGenome1 = self.dedupDir + 'seq.fasta'
        if not os.path.exists(self.mapGenome1):
            self.logUpdate('compute mappability')
            with open(self.deDupGenome) as f:
                with open(self.mapGenome1, 'w') as w:
                    for title, seq in SimpleFastaParser(f):
                        w.write(('>' + title).split()[0].strip() + '\n')
                        w.write(seq + '\n')

    def make1stProbe(self):
        # COMPUTEMAPPABILITY
        select = '-S ' + self.thermoProCoords
        self.logUpdate('[COMPUTE MAPPABILITY]')
        comMap1 = self.utils.comMappability(self.mapGenome1, self.indexDir1, self.error1, self.probLen1, self.mapDir1, select)
        self.logUpdate('[DEDUPLICATE GENMAP RESULT]')
        uniqComMap1 = self.utils.deDuplicateProbesCSV(comMap1, self.mapDir1 + 'uniq.genmap.csv')
        # MINIMIZE 1ST PROBE SET
        self.minimizedProbeSetResult1 = self.minimizingDir1 + 'result'
        minimizedProbeSetBed1 = self.minimizingDir1 + 'result.bed'
        if not os.path.exists(self.minimizedProbeSetResult1):
            self.logUpdate("[minimize probe set]")
            msg, err = self.utils.setCover(self.coverage, self.probLen1, 0.9, 11, 1, uniqComMap1, self.deDupGenome)
            print(msg)
            with open(self.minimizedProbeSetResult1, 'w') as w:
                w.write(msg)
            if not os.path.exists(minimizedProbeSetBed1):
                genomeAndIdx = dict()
                with open(self.lookup1) as f1:
                    for line in f1:
                        genome = line.split()[0].strip()
                        idx = line.split()[1].strip()
                        genomeAndIdx[idx] = genome
                with open(self.minimizedProbeSetResult1) as f2:
                    with open(minimizedProbeSetBed1, 'w') as w:
                        for line in f2:
                            kmers = line.split(';')[1].strip()
                            idx = line.split(';')[0].split(',')[0]
                            pos = line.split(';')[0].split(',')[1]
                            w.write('\t'.join([genomeAndIdx[idx], pos, str(int(pos)+self.probLen1), kmers])+'\n')
        self.utils.getSubseqFasta(minimizedProbeSetBed1, self.deDupGenome, self.workDir + 'probe1.fa')

    def make2ndWindow(self):
        probe2InputFasta = self.inputDir2 + 'seq.fa'
        with open(self.deDupGenome) as f1:
            seqName = dict()
            seqs = dict()
            seqIdx = 0
            for title, sequence in SimpleFastaParser(f1):
                header = ('>' + title).split()
                seqName[seqIdx] = header[0]
                seqs[seqIdx] = sequence
                seqIdx += 1
        #
        with open(self.minimizedProbeSetResult1) as f2:  # 1st setcover result
            probPos = dict()
            probIdx = 1
            matchedKmers = sum([line.split(';')[1].split('|') for line in f2.readlines()], [])
            for kmer in matchedKmers:
                genome = int(kmer.split(',')[0])
                probeStart = int(kmer.split(',')[1]) + 1
                probeEnd = probeStart + self.probLen1
                probPos[probIdx] = [genome, probeStart, probeEnd]
                probIdx += 1
        #
        with open(probe2InputFasta, 'w') as w:  # seq.fa
            for key in probPos:
                seqIdx = probPos[key][0]
                probeStart = probPos[key][1]
                probeEnd = probPos[key][2]
                inputSequence = seqs[seqIdx].strip()
                seqLen = len(inputSequence)
                start = probeStart - self.window if probeStart > self.window else 1
                end = probeStart + self.probLen1 + self.window if (probeStart + self.probLen1 + self.window < seqLen) else seqLen
                maskMap = sum(
                    [list(range(probPos[k][1]-1, probPos[k][2]-1)) for k in probPos if probPos[k][0] == seqIdx], []
                )
                outputSequence = ''.join(
                    ['N' if pos in maskMap else inputSequence[pos] for pos in range(start - 1, end)]
                )
                outputHeader = ':'.join([seqName[seqIdx], str(probeStart - 1), str(probeEnd - 1)])
                outputHeader = '\t'.join([outputHeader, str(start), str(end), str(seqLen)])
                w.write(outputHeader + '\n')
                w.write(outputSequence + '\n')
        self.lookup2 = self.makeLookup(probe2InputFasta, self.workDir + 'name.lookup.2')

    def make2ndProbe(self):
        # COMPUTEMAPPABILITY
        comMap2 = self.utils.comMappability(
            self.deDupGenome, self.indexDir2, self.error2, self.probLen2, self.mapDir2,  '-S ' + self.thermoProCoords
        )
        uniqComMap2 = self.utils.deDuplicateProbesCSV(comMap2, self.mapDir2 + 'uniq.genmap.csv')
        # MINIMIZE 1ST PROBE SET
        self.minimizedProbeSetResult2 = self.minimizingDir2 + 'result'
        minimizedProbeSetBed2 = self.minimizingDir2 + 'result.bed'
        self.minimizedProbeSetResult2 = self.minimizingDir2 + 'result'
        if not os.path.exists(self.minimizedProbeSetResult2):
            self.logUpdate("minimize probe set")
            # msg, err = self.utils.setCover(self.coverage, self.probLen2, 0.99, 20, 10, uniqComMap2, self.deDupGenome)
            msg, err = self.utils.setCover(self.coverage, 1, 0.99, 20, 10, uniqComMap2, self.deDupGenome)
            with open(self.minimizedProbeSetResult2, 'w') as w:
                w.write(msg)
            if not os.path.exists(minimizedProbeSetBed2):
                genomeAndIdx = dict()
                with open(self.lookup1) as f1:
                    for line in f1:
                        genome = line.split()[0].strip()
                        idx = line.split()[1].strip()
                        genomeAndIdx[idx] = genome
                with open(self.minimizedProbeSetResult2) as f2:
                    with open(minimizedProbeSetBed2, 'w') as w:
                        for line in f2:
                            kmers = line.split(';')[1].strip()
                            idx = line.split(';')[0].split(',')[0]
                            pos = line.split(';')[0].split(',')[1]
                            w.write('\t'.join([genomeAndIdx[idx], pos, str(int(pos)+self.probLen2), kmers])+'\n')
        self.utils.getSubseqFasta(minimizedProbeSetBed2, self.deDupGenome, self.workDir + 'probe2.fa')
        return

    def excute(self):
        self.getPars()
        self.checkArgs()
        self.makeWorkDir()
        self.logUpdate(' '.join(['{} {}'.format(i[0], i[1]).strip() for i in self.args]))
        self.filterInputData()
        self.make1stProbe()
        self.make2ndWindow()
        self.make2ndProbe()
        return

    def printUsage(self):
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


class SNP:
    args = []
    shortParmas = 'hr:a:s:p:m:o:'
    longParams = ['reference=', 'annotation=', 'strain-fasta=', 'positions=', 'mutation-list=', 'output=', 'help']
    utils = ProbeitUtils()
    refGenome = ''
    refGenomeAnnot = ''
    strGenome = ''
    posList = []
    snpList = []
    workDir = ''
    window = 200

    probLen1 = 40
    probLen2 = 20

    # 1, 1, 0.99, 20, 10

    # seqIdentity = 0.97
    # seqIdProbe = 0.90
    # error1 = 0
    # error2 = 1
    # coverage = 1
    # needCluster = True
    # window = 200

    def __init__(self, args):
        self.args = getopt.getopt(args, self.shortParmas, self.longParams)[0]

    def getArgList(self, value, isInt=False):
        if isInt:
            return [int(i.strip()) for i in value.split(',')]
        else:
            return value.split(',')

    def getPars(self):
        for opt, val in self.args:
            if opt in ('-h', '--help'):
                self.printUsage()
            try:
                self.refGenome = val if opt in ('-r', '--reference') else self.refGenome
                self.refGenomeAnnot = val if opt in ('-a', '--annotation') else self.refGenomeAnnot
                self.strGenome = val if opt in ('-s', '--strain-fasta') else self.strGenome
                self.workDir = val if opt in ('-o', '--output') else self.workDir
                self.posList = self.getArgList(val, isInt=True) if opt in ('-p', '--positions') else self.posList
                self.snpList = self.getArgList(val) if opt in ('-m', '--mutation-list') else self.snpList
            except Exception:
                print("Your arguments: snp {}".format(self.utils.getUserArgs(self.args)))
                self.printUsage()
        return

    def checkArgs(self):
        print("Your arguments: {}".format('snp ' + self.utils.getUserArgs(self.args)))
        message = "{}You didn't input a proper argument for '{}' parameter or missed it. "
        isBadArguments = False
        if not self.refGenome:
            print(message.format('[ERR]', '--reference'))
            isBadArguments = True
        if not self.strGenome:
            print(message.format('[ERR]', '--strain-fasta'))
            isBadArguments = True
        if not self.posList:
            print(message.format('[ERR]', '--positions'))
            isBadArguments = True
        if not self.snpList:
            print(message.format('[ERR]', '--mutaion-list'))
            isBadArguments = True
        if not self.workDir:
            print(message.format('[ERR]', '--output'))
            isBadArguments = True
        if isBadArguments:
            self.printUsage()
        else:
            print('You input proper arguments.' if self.refGenomeAnnot else message.format('[warn]', '--annotation'))

    def logUpdate(self, msg):
        with open(self.log, 'a') as w:
            w.write(msg+'\n')
        print(msg)

    def makeWorkDir(self):
        self.workDir = self.workDir + os.path.sep
        self.tempDir = self.workDir + 'temp' + os.path.sep
        self.log = self.workDir + 'log.txt'
        self.indexDir = self.workDir + 'index'
        print(os.path.exists(self.workDir))
        self.utils.delDir(self.workDir)
        print(os.path.exists(self.workDir))
        os.makedirs(self.workDir)
        os.makedirs(self.tempDir)

    def getStrKmerNearSNP(self, mutation, seqWithSNP):
        searchProbe = '{}blast.fa'.format(self.workDir)
        blastOutput = '{}blast.tsv'.format(self.workDir)
        with open(searchProbe, 'w') as w:
            w.write('>{}\n{}\n'.format(mutation, seqWithSNP))
        blastOutput = self.utils.doBlastSearch(self.workDir, searchProbe, self.strGenome, 'blast.tsv')
        return blastOutput

    def parseBlastResult(self, blastResult):
        found = -1
        groupedDf = blastResult.groupby(['WTsequence', 'STsequence', 'locSNP', 'SNPbyNT'])
        wtSequence, stSequence, locSNP, ntSNP = '', '', '', ''
        for i in groupedDf:
            if len(i[1]) > found and i[0][1][i[0][2]] == i[0][3][-1]:
                found = len(i[1])
                wtSequence = i[0][0]
                stSequence = i[0][1]
                locSNP = i[0][2]
                ntSNP = i[0][3]
            return wtSequence, stSequence, ntSNP, locSNP, found

    def parseMutation(self, mutation):
        parsedMutation = mutation.split(':')
        if len(parsedMutation) == 2 and 'nt' == parsedMutation[0]:
            mutType, snp = tuple(parsedMutation)
            orf = None
        elif len(parsedMutation) == 3 and 'aa' == parsedMutation[0]:
            mutType, orf, snp = tuple(parsedMutation)
            orf = orf.replace('orf', 'ORF').replace('Orf', 'ORF')
        else:
            mutType, orf, snp = (None, None, None)
        return mutType, orf, snp

    def getOrfStartPos(self, annotation, orf):
        with open(annotation) as annotationFile:
            for line in annotationFile:
                if '#' not in line and 'gene=' + orf in line and line.split()[2] == 'gene':
                    return int(line.split()[3])
            else:
                return -1

    def getReferenceSeq(self, refGenomeFasta):
        for h, s in SimpleFastaParser(open(refGenomeFasta)):
            return s
        return ''

    def checkCodon(self, threeMer):
        try:
            Seq(threeMer).translate()
            return True
        except Exception:
            return False

    def setProbesByPos(self):
        self.probesByPos = {pos: [] for pos in self.posList + [-1]}
        minPos = min(self.posList)
        maxPos = max(self.posList)
        refSeq = self.getReferenceSeq(self.refGenome)
        if not refSeq:
            self.logUpdate('[ERR]Failure to get reference sequence from reference genome.')
            self.printUsage()
        for snp in self.snpList:
            self.logUpdate('')
            self.logUpdate('SNP {}'.format(snp))
            mutType, orf, mutation = self.parseMutation(snp)
            print(snp, mutType, orf, mutation)
            if mutType not in ['aa', 'nt']:
                self.logUpdate('[ERR]SNP {} has a wrong format.'.format(snp))
                continue
            if mutType == 'aa':
                if self.refGenomeAnnot == '':
                    self.logUpdate('[ERR]For Amino Acid based SNPs reference annotation needed.')
                    continue
                orfStartPos = self.getOrfStartPos(self.refGenomeAnnot, orf)
                if orfStartPos == -1:
                    continue
                    self.logUpdate('[ERR]Failure to find snp {} in reference annotaion.'.format(snp))
                aa1, aa2, mutPos = mutation[0], mutation[-1], int(mutation[1:-1])
                codonStartPos = orfStartPos + (mutPos - 1) * 3 - 1
                codonEndPos = orfStartPos + mutPos * 3 - 1
                refCodon = refSeq[codonStartPos: codonEndPos]
                if aa1 != Seq(refCodon).translate():
                    self.logUpdate('Failure to find SNP {} in reference genome'.format(snp))
                    continue
                seqWithSNP = refSeq[codonStartPos - (maxPos-1): codonEndPos + (self.probLen1 - minPos)]
                strainKmerNearSNP = self.getStrKmerNearSNP(mutation, seqWithSNP)
                df = pd.read_csv(strainKmerNearSNP, sep='\t', header=None)
                df.columns = ['subGenome', 'SNPbyAA', 'match', 'STsequence']
                try:
                    df = df[df.STsequence.apply(lambda x: len(x) == len(seqWithSNP))]
                    df['STcodon'] = df.STsequence.apply(lambda x: x[maxPos - 1:maxPos + 2])
                    df = df[df.STcodon.apply(lambda x: self.checkCodon(x))]
                    df = df[df.STcodon.apply(lambda x: Seq(x).translate() == aa2)]
                    df['WTcodon'] = refCodon
                    df['diffNT'] = df.apply(lambda x: [i for i in range(len(x[4])) if x[4][i] != x[5][i]], axis=1)
                    df['diffNT'] = df.diffNT.apply(lambda x: x[0] if len(x) == 1 else -1)
                    df = df[df['diffNT'].apply(lambda x: x in (0, 1, 2))]
                    df['locSNP'] = df['diffNT'].apply(lambda x: x + maxPos - 1)
                    df['SNPbyNT'] = df.apply(
                        lambda x: '{}{}{}'.format(x[5][x[6]], codonStartPos + x[6], x[4][x[6]]), axis=1
                    )
                    df['WTsequence'] = seqWithSNP
                except:
                    self.logUpdate('[ERR]Failure to find snp {} in the strain genome or reference genome'.format(snp))
                    continue
                wtSequence, stSequence, ntSNP, locSnp, found = self.parseBlastResult(blastResult=df)
                self.logUpdate('aa:{}:{} converted nt:{}'.format(orf, mutation, ntSNP))
                mutSeqs = ParaSeqs(ntSNP, '{}:{}'.format(orf, mutation), wtSequence, stSequence, mutLoc=locSnp)
            else:
                nt1, nt2, snpPos = mutation[0], mutation[-1], int(mutation[1:-1])
                refNT = refSeq[snpPos]
                if nt1 != refNT:
                    self.logUpdate('[ERR]Failure to find SNP {} in reference genome'.format(snp))
                    continue
                seqWithSNP = refSeq[snpPos - (maxPos - 1):snpPos + 1 + (self.probLen1 - minPos)]
                blastOutput = self.getStrKmerNearSNP(mutation, seqWithSNP)
                df = pd.read_csv(blastOutput, sep='\t', header=None)
                df.columns = ['subGenome', 'SNPbyNT', 'match', 'STsequence']
                df['WTsequence'] = seqWithSNP
                df['locSNP'] = maxPos - 1
                df = df[df.STsequence.apply(lambda x: len(x) == len(seqWithSNP))]
                df = df[df.STsequence.apply(lambda x: x[maxPos - 1]) == nt2]
                wtSequence, stSequence, ntSNP, locSnp, found = self.parseBlastResult(blastResult=df)
                mutSeqs = ParaSeqs(ntSNP, '', wtSequence, stSequence, mutLoc=locSnp)
            if found == 0:
                self.logUpdate('[ERR]Failure to find SNP {} in strain genome'.format(snp))
                continue
            self.probesByPos[-1].append(mutSeqs)
            for pos in self.posList:
                wtProbe, stProbe = mutSeqs.getProbsWithPos(pos)
                paraSeq = ParaSeqs(mutSeqs.ntSnp, mutSeqs.aaSnp, wtProbe, stProbe, found=found, probLen=self.probLen1)
                self.probesByPos[pos].append(paraSeq)

    def make1stProbe(self):
        probeLines = []
        # print(self.probesByPos)
        for pos in self.posList:
            probCSV = '{}pos{}.csv'.format(self.workDir, pos)
            csvWriter = open(probCSV, 'w')
            csvWriter.write('WT sequence,ST sequence,found,ntSNP,aaSNP\n')

            for probs in self.probesByPos[pos]:
                csvWriter.write(
                    '{},{},{},{},{}\n'.format(probs.wtSeq, probs.stSeq, probs.found, probs.ntSnp, probs.aaSnp))
                probeLines.append(
                    '>{}{};{}\n'.format(probs.ntSnp, '=' + probs.aaSnp if probs.aaSnp else '', pos))
                probeLines.append('{}\n'.format(probs.stSeq))
            csvWriter.close()
        with open('{}{}mer.fa'.format(self.workDir, self.probLen1), 'w') as fastaWriter:
            fastaWriter.writelines(probeLines)

    def makeProbMaskCoords(self, inputDF, outputBed):
        inputDF['seqID'] = inputDF['seqID'].apply(lambda x: x.split(';')[0])
        inputDF.to_csv(outputBed, sep='\t', header=False, index=False)
        return outputBed

    def makeWindowCoordBed(self, inputDF, outputBed):
        inputDF[1] = inputDF[1].apply(lambda x: x - self.window)
        inputDF[2] = inputDF[2].apply(lambda x: x + self.window)
        inputDF.to_csv(outputBed, sep='\t', header=False, index=False)

    def make2ndWindow(self):
        # FILE NAMES
        snpNearprobes = '{}SNP.pattern.fasta'.format(self.workDir)
        snpMaskedBed = '{}SNP.masked.bed'.format(self.workDir)
        # USING PROBEs AND STRAIN GENOME MAKE PROBEs MASK TSV
        with open(snpNearprobes, 'w') as w:
            w.writelines(['>{}{}\n{}\n'.format(p.ntSnp, '=' + p.aaSnp if p.aaSnp else '', p.stSeq) for p in self.probesByPos[-1]])
        self.lookupTSV = self.utils.getPatternPosition(snpNearprobes, self.strGenome, '{}lookup.tsv'.format(self.workDir))
        snpMaskedBed = self.makeProbMaskCoords(
            pd.read_csv(self.lookupTSV, sep='\t')[['seqID', 'start', 'end']], snpMaskedBed
        )
        # MAKE PROBEs MAKSED FASTA
        self.probe2InputFasta = self.utils.getWindowFasta(
            self.strGenome, snpMaskedBed,
            '{}SNP.masked.fasta'.format(self.workDir),
            '{}window.bed'.format(self.workDir),
            '{}2nd.input.fasta'.format(self.workDir)
        )
        lookup = self.utils.makeLookup(self.probe2InputFasta, self.workDir + 'name.lookup')
        self.utils.delDir(self.indexDir)
        mappedCSV = self.utils.comMappability(self.probe2InputFasta, self.indexDir, 1, self.probLen2, self.workDir)
        dedupMappedCSV = self.utils.deDuplicateProbesCSV(mappedCSV, '{}uniq.genmap.csv'.format(self.workDir))
        stdOut, stdErr = self.utils.setCover(1, 1, 0.99, 20, 10, dedupMappedCSV, self.probe2InputFasta)
        setcoverResult = self.workDir + 'result'
        self.logUpdate(stdOut)
        with open(setcoverResult, 'w') as w:
            w.write(stdOut)
        self.setcoverResultBed = self.utils.makeSetcoverResultBed(lookup, setcoverResult, self.workDir + 'result.bed')

    def set2ndProbe(self):
        tempSecondProbeFasta = '{}temp.{}mer.fasta'.format(self.workDir, self.probLen2)
        secondProbeFasta = '{}{}mer.fa'.format(self.workDir, self.probLen2)
        self.utils.getSubseqFasta(self.setcoverResultBed, self.probe2InputFasta, tempSecondProbeFasta)
        maskDF = pd.read_csv(self.lookupTSV, sep='\t')
        maskDF['lookup'] = maskDF.apply(lambda x: '{}:{}-{}'.format(x[0], x[4] - self.window, x[5] + self.window), axis=1)
        kmers = list(maskDF['patternName'])
        # PARSING TEMP 2ND FASTA
        w = open(secondProbeFasta, 'w')
        for h, s in SimpleFastaParser(open(tempSecondProbeFasta)):
            p = re.compile('[0-9]+,')
            kmerIndex = p.findall(h)
            soveredSNPs = ':'.join(list((set([kmers[int(i[:-1])] for i in kmerIndex]))))
            w.write('>{}\n'.format(soveredSNPs))
            w.write(s + '\n')

    def excute(self):
        self.getPars()
        self.checkArgs()
        self.makeWorkDir()
        self.logUpdate('Your arguments: snp ' + self.utils.getUserArgs(self.args) + '\n')
        self.logUpdate("MAKE 1ST PROBES")
        self.setProbesByPos()
        self.make1stProbe()
        self.logUpdate("MAKE 2ND PROBES")
        self.make2ndWindow()
        self.set2ndProbe()
        self.logUpdate("DONE")
        return

    def printUsage(self):
        print(" -r|--reference Reference Genome [FASTA FILE]")
        print(" -s|--strain-fasta Strain Genome [FASTA FILE]")
        print(" -o|--output Directory for Output Data [DIR NAME]")
        print(" -p|--positions Positions for First Probes [Comma Separated INT ARRAY]")
        print(" -m|--mutation-list SNPs Interest Strain has [Comma Separated SNP ARRAY]")
        print("OPTIONAL")
        print(" -a|--annotation Annotation File of Reference Genome [GFF FILE]")
        quit()


def main(args):
    probeit = Probeit(args)
    probeit.checkArgs()
    quit()


if __name__ == '__main__':
    main(sys.argv[1:])
