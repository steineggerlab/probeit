#!/usr/bin/env python
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

    def getProbesWithPos(self, pos):
        return (
            self.wtSeq[self.mutLoc - pos + 1:self.mutLoc - pos + self.probLen + 1],
            self.stSeq[self.mutLoc - pos + 1:self.mutLoc - pos + self.probLen + 1]
        )


class ProbeitUtils:
    @staticmethod
    def runCommand(command, verbose=False):
        print('[CLI] '+command)
        if verbose:
            commandList = command.split()
            sp = subprocess.Popen(commandList, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            stdout, stderr = sp.communicate()
            return stdout.decode('UTF-8'), stderr.decode('UTF-8')
        else:
            os.system(command)
            return

    @classmethod
    def delDir(cls, waste):
        command = "rm -rf " + waste
        cls.runCommand(command)

    # TO CALL SEQKIT MODULES
    @classmethod
    def sortFasta(cls, inputFasta, outputFasta):
        command = "seqkit sort {} > {} -w 0".format(inputFasta, outputFasta)
        cls.runCommand(command)

    @classmethod
    def renameFasta(cls, oldFasta, newFasta):
        command = "seqkit rename {} > {} -w 0".format(oldFasta, newFasta)
        cls.runCommand(command)

    @classmethod
    def getSubseqFasta(cls, coordinateBed, inputFasta, outputFasta):
        command1 = "seqkit subseq --quiet --bed {} {} > {}".format(coordinateBed, inputFasta, outputFasta)
        cls.runCommand(command1)

    @classmethod
    def getPatternPosition(cls, patternFasta, genomeFasta, positonsTSV):
        command = "seqkit locate -f {} {} > {}".format(patternFasta, genomeFasta, positonsTSV)
        cls.runCommand(command)
        df = pd.read_csv(positonsTSV, sep='\t')
        df.sort_values(by=list(df.columns), ignore_index=True).to_csv(positonsTSV, sep='\t', index=False)
        return positonsTSV

    # TO CALL GENMAP MODULES
    # COMPUTE MAPPABILITY
    @classmethod
    def computeMappability(cls, fasta, indexDir, error, kmer, outputDir, selectingOption=''):
        cls.delDir(indexDir)
        command0 = " > /dev/null"
        command1 = "genmap index -F {} -I {}".format(fasta, indexDir)
        cls.runCommand(command1 + command0)
        command2 = (
            "genmap map --no-reverse-complement -E {} {} --csv -K {} -t -b --frequency-large -I {} -O {}"
            .format(error, selectingOption, kmer, indexDir, outputDir)
        )
        cls.runCommand(command2 + command0)
        return outputDir + '.'.join(fasta.split(os.path.sep)[-1].split('.')[:-1] + ['genmap', 'csv'])

    # TO CALL BEDTOOLS MODULES
    @classmethod
    def getSubtractedBed(cls, positiveBed, negativeBed, outputBed):
        command = "bedtools subtract -a {} -b {} > {}".format(positiveBed, negativeBed, outputBed)
        cls.runCommand(command)
        return outputBed

    @classmethod
    def getWindowFasta(cls, genomeFasta, maskingBed, maskedGenomeFasta, windowBed, windowFasta, windowSize):
        command1 = "bedtools maskfasta -fi {} -bed {} -fo {}".format(genomeFasta, maskingBed, maskedGenomeFasta)
        cls.runCommand(command1)
        inputDF = pd.read_csv(maskingBed, sep='\t', header=None)
        inputDF[1] = inputDF[1].apply(lambda x: x - windowSize)
        inputDF[2] = inputDF[2].apply(lambda x: x + windowSize)
        inputDF.to_csv(windowBed, sep='\t', header=False, index=False)
        command2 = "bedtools getfasta -fi {} -bed {} > {}".format(maskedGenomeFasta, windowBed, windowFasta)
        cls.runCommand(command2)
        return windowFasta

    @classmethod
    def makeLookup(cls, windowFasta, lookup, simpleGenome=False):
        if os.path.exists(lookup):
            return lookup
        with open(windowFasta) as f:
            if simpleGenome:
                headers = [title.split()[0].strip() for title, seq in SimpleFastaParser(f)]
            else:
                headers = [title.strip() for title, seq in SimpleFastaParser(f)]
        lookupLines = [headers[i] + '\t' + str(i) + '\n' for i in range(len(headers))]
        with open(lookup, 'w') as w:
            w.writelines(lookupLines)
        return lookup

    @classmethod
    def searchSNPs(cls, workDir, inputFasta, strGenomeFasta, resultTSV, kmer=12):
        searchDir = workDir + 'search' + os.path.sep
        tempDir = searchDir + 'temp' + os.path.sep
        resultTSV = workDir + resultTSV
        if os.path.isdir(searchDir):
            ProbeitUtils.delDir(searchDir)
        os.makedirs(searchDir)
        searchdb = searchDir + 'searchDB'
        strdb = searchDir + 'strainDB'
        aln = searchDir + 'mmseqs.aln'
        command1 = 'mmseqs createdb {} {}'
        command2 = 'mmseqs createdb {} {}'
        command3 = 'mmseqs search {} {} {} {} --search-type 3 -k {}'
        command4 = 'mmseqs convertalis {} {} {} {} --format-output target,query,tseq,tstart,tend --search-type 3'
        out1, err1 = cls.runCommand(command1.format(inputFasta, searchdb), verbose=True)
        out2, err2 = cls.runCommand(command2.format(strGenomeFasta, strdb), verbose=True)
        out3, err3 = cls.runCommand(command3.format(searchdb, strdb, aln, tempDir, kmer), verbose=True)
        out4, err4 = cls.runCommand(command4.format(searchdb, strdb, aln, resultTSV), verbose=True)
        df = pd.read_csv(resultTSV, sep='\t', header=None)
        df.columns = ['substr', 'snp', 'strseq', 'start', 'end']
        df['aln'] = df.apply(lambda x: x[2][int(x[3]-1):int(x[4])], axis=1)
        df['len'] = df.aln.apply(lambda x: len(x)-1)
        df = df[['substr', 'snp', 'len', 'aln']]
        df.to_csv(resultTSV, header=False, index=False, sep='\t')
        print(err1 + err2 + err3 + err4)
        return resultTSV, out1 + out2 + out3 + out4

    @classmethod
    def clusterGenome(cls, inputFasta, outputFasta, outputDir, seqIdentity):
        command = ' '.join(
            [
                "mmseqs easy-linclust", inputFasta, outputFasta, outputDir,
                "-v 3 --kmer-per-seq-scale 0.5 --kmer-per-seq 1000 --min-seq-id", str(seqIdentity),
                "--cov-mode 1 -c 0.95 --remove-tmp-files 0"
            ]
        )
        stdout, stderr = cls.runCommand(command, verbose=True)
        return stdout, stderr

    @classmethod
    def searchNegative(cls, output, negative, maskOutput, outputDir, seqInProbe):
        command = ' '.join(
            [
                "mmseqs easy-search", output, negative, maskOutput, outputDir+"tmp",
                "-v 3 --spaced-kmer-mode 0 -k 13 --mask 0 -c 0.9 --min-seq-id",
                str(seqInProbe), "--cov-mode 2 --alignment-mode 4 --search-type 3",
                "--format-output query,target,pident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue"
            ]
        )
        stdOut, stdErr = cls.runCommand(command, verbose=True)
        print(stdOut)
        print(stdErr)
        return stdOut, stdErr

    @classmethod
    def deDuplicateMapCSV(cls, inputCSV, outputCSV):
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

    @classmethod
    def setCover(cls, coverage, length, eStop, dist, reps, mapCSV, genome):
        filePath = os.path.sep.join(os.path.realpath(__file__).split(os.path.sep)[:-1])
        setcoverPath = '{}{}{}{}{}'.format(filePath, os.path.sep, 'setcover', os.path.sep, 'setcover')
        command = " -c {} -l {} -p {} -d {} -i {} {} {}".format(coverage, length, eStop, dist, reps, mapCSV, genome)
        stdOut, stdErr = cls.runCommand(setcoverPath + command, verbose=True)
        return stdOut, stdErr

    @classmethod
    def makeMinzProbeBed(cls, lookup, setcoverResult, setcoverResultBed, probeLen):
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
                    w.write('\t'.join([genomeAndIdx[genome], pos, str(int(pos) + probeLen), matchedKmers]) + '\n')
        return setcoverResultBed

    @classmethod
    def getUserArgs(cls, args):
        return ' '.join(['{} {}'.format(i[0], i[1]) if len(i) == 2 else i for i in args])

    @classmethod
    def sortFile(cls, inputFile, outputFile):
        command = 'sort {} > {}'.format(inputFile, outputFile)
        cls.runCommand(command)
        cls.runCommand(command)

    @staticmethod
    def simplifyFastaHeaders(inputFasta, outputFasta):
        with open(inputFasta) as f:
            with open(outputFasta, 'w') as w:
                for title, seq in SimpleFastaParser(f):
                    w.write(('>' + title).split()[0].strip() + '\n')
                    w.write(seq + '\n')

    @classmethod
    def makeProbe(
            cls, indexDir, mapDir, mapError, probLen, window, genome, lookup, probe, minzResult, minzBed,
            minzCovered, minzLen, minzEarlyStop, minzSimScore, minzRepeats, selector='',
    ):
        # COMPUTEMAPPABILITY
        message = ''
        message += '[INFO]compute mappability'
        comMap = ProbeitUtils.computeMappability(window, indexDir, mapError, probLen, mapDir, selector)
        message += '[INFO]deduplicate genmap result'
        uniqComMap = ProbeitUtils.deDuplicateMapCSV(comMap, mapDir + 'uniq.genmap.csv')

        # SETCOVER
        message += "[INFO]minimize probe set"
        msg, err = ProbeitUtils.setCover(
            minzCovered,
            minzLen,
            minzEarlyStop,
            minzSimScore,
            minzRepeats,
            uniqComMap,
            genome
        )
        message += msg
        message += err
        with open(minzResult, 'w') as w:
            w.write(msg)

        # MAKE PROBEs
        ProbeitUtils.makeMinzProbeBed(lookup, minzResult, minzBed, probLen)
        ProbeitUtils.getSubseqFasta(minzBed, genome, probe)
        return message

    @staticmethod
    def parseGenmapPattern(header):
        p1 = re.compile('[0-9]+,')
        p2 = re.compile(',[0-9]+')
        return [int(i[:-1]) for i in p1.findall(header)], [int(i[1:]) for i in p2.findall(header)]


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
            print('CURRENT: ', os.getcwd())
            self.subWork = SNP(self.args[1:])
            self.subWork.excute()
            return
        else:
            self.printUsage()
            return

    @staticmethod
    def printUsage():
        print("PROBEIT")
        print('usage ./probeit.py <module> [<args>]')
        print("WORKFLOWS")
        print("posnegset: make two-sets probes with positive and negative sets")
        print("snp: make two-sets probes with wildtype genome, strain genome and SNPs")
        quit()


class PosNegSet:
    args = []
    shortParams = 'hp:n:o:'
    longParams = [
        # usage
        'help',
        # required
        'positive=', 'negative=', 'output=',
        # optional
        'not-cluster-positive', 'probe-len1=', 'probe-len2=', 'probe-error1=', 'probe-error2=',
        'minimizing-covered1=', 'minimizing-covered2=', 'minimizing-repeats1=', 'minimizing-repeats2',
        'minimizing-earlystop-criteria1=', 'minimizing-earlystop-criteria2=', 'window-size=',
        # hidden
        'dedup-genome-identity=', 'rid-negative-identity=', 'minimizing-simscore1=', 'minimizing-simscore2='
    ]
    inputGenome = ''
    negGenome = ''
    workDir = ''
    probeLen1 = 40
    probeLen2 = 20
    deDupGenomeIdentity = 0.97  # for mmseqs linclust
    ridNegIdentity = 0.90   # for mmseqs easy-search
    mapError1 = 0
    mapError2 = 1
    minzCovered1 = 1
    minzCovered2 = 1
    needCluster = True
    windowSize = 200
    minzEarlyStop1 = 0.9
    minzSimScore1 = 11
    minzRepeats1 = 1
    minzEarlyStop2 = 0.99
    minzSimScore2 = 20
    minzRepeats2 = 10

    def __init__(self, args):
        self.args = getopt.getopt(args, self.shortParams, self.longParams)[0]

    def getPars(self):
        for opt, val in self.args:
            if opt in ('-h', '--help'):
                self.printUsage()
            try:
                self.inputGenome = str(val) if opt in ('-p', '--positive') else self.inputGenome
                self.negGenome = str(val) if opt in ('-n', '--negative') else self.negGenome
                self.workDir = str(val) if opt in ('-o', '--output') else self.workDir
                # optional args
                self.windowSize - int(val) if opt == '--window-size' else self.windowSize
                self.needCluster = False if opt == '--not-cluster-positive' else self.needCluster
                self.probeLen1 = int(val) if opt == '--probe-len1' else self.probeLen1
                self.probeLen2 = int(val) if opt == '--probe-len2' else self.probeLen2
                self.mapError1 = int(val) if opt == '--probe-error1' else self.mapError1
                self.mapError2 = int(val) if opt == '--probe-error2' else self.mapError2
                self.minzCovered1 = int(val) if opt == '--minimizing-covered1' else self.minzCovered1
                self.minzCovered2 = int(val) if opt == '--minimizing-covered2' else self.minzCovered2
                self.minzRepeats1 = int(val) if opt == '--minimizing-repeats1' else self.minzRepeats1
                self.minzRepeats2 = int(val) if opt == '--minimizing-repeats2' else self.minzRepeats2
                self.minzEarlyStop1 = float(val) if opt == '--minimizing-earlystop-criteria1' else self.minzEarlyStop1
                self.minzEarlyStop2 = float(val) if opt == '--minimizing-earlystop-criteria2' else self.minzEarlyStop2
                # HIDDEN args
                # identity in mmseqs cluster
                self.deDupGenomeIdentity = float(val) if opt == '--dedup-genome-identity' else self.deDupGenomeIdentity
                # identity in mmeqs seqrch
                self.ridNegIdentity = float(val) if opt == '--rid-negative-identity' else self.ridNegIdentity
                # setcover similarity
                self.minzSimScore1 = int(val) if opt == '--minimizing-simscore1' else self.minzSimScore1
                self.minzSimScore2 = int(val) if opt == '--minimizing-simscore2' else self.minzSimScore2

                
            except Exception as e:
                print(e)
                print("Your arguments: posnegset {}".format(ProbeitUtils.getUserArgs(self.args)))
                self.printUsage()

    def checkArgs(self):
        print("Your arguments: {}".format('snp ' + ProbeitUtils.getUserArgs(self.args)))
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
        self.minzDir1 = self.workDir + 'setcover_probe1' + os.path.sep
        self.inputDir2 = self.workDir + 'input_probe2' + os.path.sep
        self.mapDir2 = self.workDir + 'mapping_probe2' + os.path.sep
        self.indexDir2 = self.workDir + 'index_probe2' + os.path.sep
        self.minzDir2 = self.workDir + 'setcover_probe2' + os.path.sep
        self.log = self.workDir + 'log.txt'
        dirList = [
            self.workDir, self.inputDir, self.dedupDir, self.maskingDir,
            self.thermoFilteringDir, self.mapDir1, self.minzDir1,
            self.inputDir2, self.mapDir2, self.minzDir2
        ]
        for directory in dirList:
            if not os.path.exists(directory):
                os.makedirs(directory)
        return

    def logUpdate(self, msg):
        print(msg)
        with open(self.log, 'a') as w:
            w.write(msg + '\n')

    def copyFile(self, original, copy):
        if not os.path.exists(copy):
            shutil.copy(original, copy)
            self.logUpdate('[INFO]{} is copied to {}'.format(original, copy))
        return copy

    def thrmoFilter(
            self, directory, inputBed, ligInput, capInput="", doSaveUnfiltered=True,
    ):
        probeLen = self.probeLen1
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
            .query('p1_tm >= {} & p2_tm >= {}'.format(minProbeTm, minProbeTm))
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
        return ProbeitUtils.getSubtractedBed(inputBed, absThermoProperOutput, directory + 'primer3.neg.bed')

    def filterInputData(self):
        clustName = self.dedupDir + 'genomes.clu'
        # CLUTER POSITIVE GENOME
        posGenome = self.copyFile(original=self.inputGenome, copy=self.inputDir + 'genome.fa')
        self.deDupGenome = clustName + '_rep_seq.fasta'
        if self.needCluster:
            self.logUpdate('[INFO]deduplicate positive fasta')
            if not os.path.exists(self.deDupGenome):
                tempDir = self.dedupDir + 'temp' + os.path.sep
                msg, err = ProbeitUtils.clusterGenome(posGenome, clustName, tempDir, self.deDupGenomeIdentity)
                self.logUpdate(msg + err)
                ProbeitUtils.sortFasta(self.deDupGenome, self.dedupDir + 'sorted.fasta')
                ProbeitUtils.renameFasta(self.dedupDir + 'sorted.fasta', self.deDupGenome)
            else:
                self.copyFile(posGenome, self.deDupGenome)
        self.lookup1 = ProbeitUtils.makeLookup(self.deDupGenome, self.workDir + 'genome.lookup', simpleGenome=True)
        
        # MAKE 40 MERS FROM POSITIVE GENOME
        posProbes = self.maskingDir + 'probes.fa'
        self.logUpdate('[INFO]get {}mer probes from positive genome'.format(self.probeLen1))
        f = open(self.deDupGenome)
        w = open(posProbes, 'w')
        for title, seq in SimpleFastaParser(f):
            header = '>' + title.split()[0].strip()
            genomeLen = len(seq.strip())
            lines = [header + '_' + str(i + 1) + '\n' + seq[i:i + self.probeLen1] + '\n' for i in
                     range(genomeLen - self.probeLen1)]
            w.writelines(lines)
        f.close()
        w.close()
        
        # MAKE 40MERS FROM BOTH THE POSAITIVE GENOME AND THE NEGATIVE GENOME COORDINATE FILE(BED)
        self.logUpdate('[INFO]remove probes found in negative genome')
        negProbesCoords = self.maskingDir + 'mmseqs.search'
        if not os.path.exists(negProbesCoords):
            ProbeitUtils.searchNegative(posProbes, self.negGenome, negProbesCoords, self.maskingDir, self.ridNegIdentity)
        
        # MAKE DEDUPLICATE POSITIVE GENOME COORDINATE FILE(BED)
        deDupGenomeCoords = clustName + '_rep_seq.bed'
        with open(self.deDupGenome) as f:
            with open(deDupGenomeCoords, 'w') as w:
                for title, seq in SimpleFastaParser(f):
                    header = title.split()[0].strip()
                    w.write(header + '\t1\t' + str(len(seq.strip())) + '\n')
        
        # COPY 40MERS FROM THE NEGATIVE GENOME COORDINATE FILE(BED) TO FILTER DIR
        simpleNegProbesCoords = self.thermoFilteringDir + 'mmseqs.txt'
        with open(negProbesCoords) as f:
            with open(simpleNegProbesCoords, 'w')as w:
                for line in f:
                    c1 = line.split()[0].strip()
                    b = c1.split('_')[-1]
                    c1 = c1.split('_')[0]
                    w.write(c1 + '\t' + str(int(b) - 1) + '\t' + str(int(b) - 1 + self.probeLen1) + '\n')
        
        # MAKE NEGATIVE REMOVED COORDINATE FILE(BED)
        negRemCoords = ProbeitUtils.getSubtractedBed(
            deDupGenomeCoords, simpleNegProbesCoords, self.thermoFilteringDir + 'crosstaxa.bed'
        )
        
        self.logUpdate('[INFO]filter probes with thermodynamic features')
        probesForTheromFilter = self.thermoFilteringDir + 'probes.fa'
        with open(posProbes) as f:
            with open(probesForTheromFilter, 'w') as w:
                for title, seq in SimpleFastaParser(f):
                    header = ('>' + title).split('_')
                    w.write('_'.join(header[:-1]) + '\t' + header[-1] + '\n')
                    w.write(seq + '\n')
                    
        # TO FILTER 40MERS by THERMODYNAMIC FEATURES
        thermoImpCoords = self.thrmoFilter(self.thermoFilteringDir, deDupGenomeCoords, probesForTheromFilter)
        self.thermoBed = ProbeitUtils.getSubtractedBed(
            negRemCoords, thermoImpCoords, self.thermoFilteringDir + 'crosstaxa.primer3.bed'
        )

        # MAKE SEQ.FA
        self.window1 = self.dedupDir + 'seq.fasta'
        if not os.path.exists(self.window1):
            self.logUpdate('[INFO]compute mappability')
            ProbeitUtils.simplifyFastaHeaders(self.deDupGenome, self.window1)

    def make2ndWindow(self):
        # MAKE SEQ.FA
        self.window2 = self.inputDir2 + 'seq.fa'
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
        with open(self.minzResult1) as f2:  # 1st setcover result
            probes = dict()
            probIdx = 1
            matchedKmers = sum([line.split(';')[1].split('|') for line in f2.readlines()], [])
            for kmer in matchedKmers:
                genome = int(kmer.split(',')[0])
                probeStart = int(kmer.split(',')[1]) + 1
                probeEnd = probeStart + self.probeLen1
                probes[probIdx] = [genome, probeStart, probeEnd]
                probIdx += 1
        #
        with open(self.window2, 'w') as w:  # seq.fa
            for key in probes:
                seqIdx = probes[key][0]
                probeStart = probes[key][1]
                probeEnd = probes[key][2]
                inputSeq = seqs[seqIdx].strip()
                seqLen = len(inputSeq)
                start = probeStart - self.windowSize if probeStart > self.windowSize else 1
                end = probeStart + self.probeLen1 + self.windowSize
                end = end if end < seqLen else seqLen
                maskMap = [list(range(probes[k][1] - 1, probes[k][2] - 1)) for k in probes if probes[k][0] == seqIdx]
                maskMap = sum(maskMap, [])
                outputSeq = ''.join(['N' if pos in maskMap else inputSeq[pos] for pos in range(start - 1, end)])
                outputHeader = ':'.join([seqName[seqIdx], str(probeStart - 1), str(probeEnd - 1)])
                outputHeader = '\t'.join([outputHeader, str(start), str(end), str(seqLen)])
                w.write(outputHeader + '\n')
                w.write(outputSeq + '\n')
        self.lookup2 = ProbeitUtils.makeLookup(self.window2, self.workDir + 'probe1.lookup', simpleGenome=True)

    def trimProbes(self):
        # MAKE PROBEs UNERSTANDABLE
        self.probe1 = self.workDir + 'probe1.fa'
        self.probe2 = self.workDir + 'probe2.fa'

        with open(self.lookup1)as f:
            genomeKeys = {int(i.split()[1]):i.split()[0] for i in f}
        with open(self.lookup2)as f:
            probe1Keys = {int(i.split()[1]):i.split()[0] for i in f}

        cnt = 0
        namesForProbes1 = {}
        lines1 = []
        for h, s in SimpleFastaParser(open(self.tempProbe1)):
            keys, pos = ProbeitUtils.parseGenmapPattern(h)
            keys = [genomeKeys[k] for k in keys]
            probes = ['{}:{}:{}'.format(keys[i], pos[i], pos[i] + 40) for i in range(len(keys))]
            for i in probes:
                namesForProbes1[i] = 'p1_{}'.format(cnt)
            lines1.append('>p1_{}\t{}\n{}\n'.format(cnt, ';'.join(probes), s))
            cnt += 1
        with open(self.probe1, 'w') as w:
            w.writelines(lines1)

        cnt = 0
        lines2 = []
        for h, s in SimpleFastaParser(open(self.tempProbe2)):
            keys, pos = ProbeitUtils.parseGenmapPattern(h)
            keys = set([namesForProbes1[probe1Keys[k]] for k in keys])
            lines2.append('>p2_{}\t{}\n{}\n'.format(cnt, ';'.join([p for p in keys]), s) )
            cnt += 1
        with open(self.probe2, 'w') as w:
            w.writelines(lines2)

    def excute(self):
        self.getPars()
        self.checkArgs()
        self.makeWorkDir()
        self.logUpdate(' '.join(['{} {}'.format(i[0], i[1]).strip() for i in self.args]))
        self.filterInputData()
        # self.make1stProbe()
        self.minzResult1 = self.minzDir1 + 'result'
        self.tempProbe1 = self.workDir + '{}mer.fa'.format(self.probeLen1)
        msg = ProbeitUtils.makeProbe(
            self.indexDir1, self.mapDir1, self.mapError1, self.probeLen1, self.window1, self.deDupGenome,
            self.lookup1, self.tempProbe1, self.minzResult1, self.minzResult1 + '.bed', self.minzCovered1,
            self.probeLen1, self.minzEarlyStop1, self.minzSimScore1, self.minzRepeats1, selector='-S ' + self.thermoBed
        )
        self.logUpdate(msg)
        self.make2ndWindow()
        self.tempProbe2 = self.workDir + '{}mer.fa'.format(self.probeLen2)
        self.minzResult2 = self.minzDir2 + 'result'
        ProbeitUtils.makeProbe(
            self.indexDir2, self.mapDir2, self.mapError2, self.probeLen2, self.window2, self.window2,
            self.lookup2, self.tempProbe2, self.minzResult2, self.minzResult2 + '.bed',
            self.minzCovered2, 1, self.minzEarlyStop2, self.minzSimScore2, self.minzRepeats2
        )
        self.trimProbes()
        ProbeitUtils.sortFasta(self.probe1, self.workDir + 'sorted1.fa')
        ProbeitUtils.sortFasta(self.probe2, self.workDir + 'sorted2.fa')
        return

    @staticmethod
    def printUsage():
        print("Usage")
        print(" -h|--help NONE")
        print("\t Show usage")
        print("python probeit.py -p [POSITIVE GENOME]-n [NEGATIVE GENOME] -o [DIR]")
        print("REQUIRED")
        print(" -p|--positive FASTA file")
        print("\t The genome which MUST be covered by the probes.")
        print(" -n|--negative FASTA file")
        print("\t The genome which MUST NOT be covered by the probes.")
        print(" -o|--output DIR")
        print("\t Output directory The Directory is automatically created by Probeit.")
        print("OPTIONAL")
        print(" --window-size INT[200]")
        print("\t size of windows for 2nd probes")
        print("--not-cluster-positive NONE")
        print("\t Use it when you **DO NOT** need to cluster positive genome")
        print("--probe-len1 INT[40]")
        print("\t Length of 1st Probes")
        print("--probe-len2 INT[20]")
        print("\t Length of 2nd Probes")
        print(" --probe-error1 INT[0]")
        print("\t The number of error allowed in 1st Probes")
        print(" --probe-error2 INT[1]")
        print("\t The number of error allowed in 2nd Probes")
        print("--minimizing-covered1 INT[1]")
        print("\t The number of times each Seqs from positive genome should be covered by 1st Probes")
        print("--minimizing-covered2 INT[1]")
        print("\t The number of times each 1st Probe should be covered by 2nd Probes")
        print("--minimizing-repeats1 INT[1]")
        print("\t The number of random iterations when minimizing 1st Probes")
        print("--minimizing-repeats2 INT[10]")
        print("\t The number of random iterations when minimizing 2nd Probes")
        print("--minimizing-earlystop-criteria1 FLOAT[0.9]")
        print("\t Proportion of Seqs from positive genome covered by 1st Probes for earlystop upon minimizing of 1st Probe")
        print("--minimizing-earlystop-criteria2 FLOAT[0.99]")
        print("\t Proportion of 1st Probes covered by 2nd Probes for earlystop upon minimizing of 2nd Probe")
        quit()


class SNP:
    args = []
    shortParmas = 'hr:a:s:p:m:o:'
    longParams = [
        # usage
        'help',
        # required
        'reference=', 'strain=', 'positions=', 'mutations=', 'output=', 'annotation=',
        # optional
        'max-window=', 'window-size=', 'probe-len1=', 'probe-len2=', 'probe-error2=',
        'minimizing-covered2=', 'minimizing-earlystop-criteria2=', 'minimizing-repeats2=',
        # hidden
        'setcover-simscore=', 'mutation-search-kmer='
    ]
    refGenome = ''
    refGenomeAnnot = ''
    strGenome = ''
    posList = []
    snpList = []
    workDir = ''
    isMaxWindow = False
    windowSize = 200
    probLen1 = 40
    probeLen2 = 20
    minzCovered2 = 1
    minzEarlyStop2 = 0.99
    minzSimScore2 = 20
    minzRepeats2 = 10
    mapError2 = 1
    mutSearchKmer = 12


    def __init__(self, args):
        self.args = getopt.getopt(args, self.shortParmas, self.longParams)[0]

    @staticmethod
    def getArgList(value, isInt=False):
        if isInt:
            return [int(i.strip()) for i in value.split(',')]
        else:
            return value.split(',')

    def getPars(self):
        for opt, val in self.args:
            if opt in ('-h', '--help'):
                self.printUsage()
            try:
                # required
                self.refGenome = val if opt in ('-r', '--reference') else self.refGenome
                self.strGenome = val if opt in ('-s', '--strain') else self.strGenome
                self.workDir = val if opt in ('-o', '--output') else self.workDir
                self.posList = self.getArgList(val, isInt=True) if opt in ('-p', '--positions') else self.posList
                self.snpList = self.getArgList(val) if opt in ('-m', '--mutations') else self.snpList
                self.refGenomeAnnot = val if opt in ('-a', '--annotation') else self.refGenomeAnnot
                # optional
                self.isMaxWindow = True if opt == '--max-window' else self.isMaxWindow
                self.windowSize = int(val) if opt == '--window-size' else self.windowSize
                self.probLen1 = int(val) if opt == '--probe-len1' else self.probLen1
                self.probeLen2 = int(val) if opt == '--probe-len2' else self.probeLen2
                self.mapError2 = int(val) if opt == '--probe-error2' else self.mapError2  # genmap map error
                self.minzCovered2 = int(val) if opt == '--minimizing-covered2' else self.minzCovered2
                self.minzEarlyStop2 = float(val) if opt == '--minimizing-earlystop-criteria2' else self.minzEarlyStop2
                self.minzRepeats2 = int(val) if opt == '--minimizing-repeats2' else self.minzRepeats2
                # hidden args
                self.minzSimScore2 = int(val) if opt == '--minimizing-simscore2' else self.minzSimScore2
                self.mutSearchKmer = int(val) if opt == '--mutation-search-kmer' else self.mutSearchKmer

            except Exception as e:
                print(e)
                print("Your arguments: snp {}".format(ProbeitUtils.getUserArgs(self.args)))
                self.printUsage()
        else:
            if not self.isMaxWindow:
                self.windowSize = self.windowSize - (max(self.posList) - min(self.posList))
        return

    def checkArgs(self):
        print("Your arguments: {}".format('snp ' + ProbeitUtils.getUserArgs(self.args)))
        message = "{}You didn't input a proper argument for '{}' parameter or missed it. "
        isBadArguments = False
        if not self.refGenome:
            print(message.format('[ERROR]', '--reference'))
            isBadArguments = True
        if not self.strGenome:
            print(message.format('[ERROR]', '--strain'))
            isBadArguments = True
        if not self.posList:
            print(message.format('[ERROR]', '--positions'))
            isBadArguments = True
        if not self.snpList:
            print(message.format('[ERROR]', '--mutaions'))
            isBadArguments = True
        if not self.workDir:
            print(message.format('[ERROR]', '--output'))
            isBadArguments = True
        if isBadArguments:
            self.printUsage()
        else:
            print('You input proper arguments.' if self.refGenomeAnnot else message.format('[WARN]', '--annotation'))

    def logUpdate(self, msg):
        print(msg)
        with open(self.log, 'a') as w:
            w.write(msg+'\n')

    def makeWorkDir(self):
        self.workDir = self.workDir + os.path.sep
        self.tempDir = self.workDir + 'temp' + os.path.sep
        self.log = self.workDir + 'log.txt'
        self.indexDir = self.workDir + 'index'
        self.blastDir = self.workDir + 'blast' + os.path.sep
        self.probe1byPosDir = self.workDir + 'probe1byPos' + os.path.sep
        self.input2Dir = self.workDir + 'input2' + os.path.sep
        self.map2Dir = self.workDir + 'mapping2' + os.path.sep
        self.minz2Dir = self.workDir + 'setcover2' + os.path.sep
        ProbeitUtils.delDir(self.workDir)
        os.makedirs(self.workDir)
        os.makedirs(self.tempDir)
        os.makedirs(self.blastDir)
        os.makedirs(self.probe1byPosDir)
        os.makedirs(self.input2Dir)
        os.makedirs(self.map2Dir)
        os.makedirs(self.minz2Dir)

    def getStrKmerNearSNP(self, mutation, seqWithSNP):
        searchProbe = '{}blast.fa'.format(self.blastDir)
        with open(searchProbe, 'w') as w:
            w.write('>{}\n{}\n'.format(mutation, seqWithSNP))
        blastOutput, msg = ProbeitUtils.searchSNPs(
            self.blastDir, searchProbe, self.strGenome, 'blast.tsv', self.mutSearchKmer
        )
        self.logUpdate(msg)
        return blastOutput

    @staticmethod
    def parseBlastResult(blastResult):
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

    @staticmethod
    def parseMutation(mutation):
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

    @staticmethod
    def getOrfStartPos(annotation, orf):
        with open(annotation) as annotationFile:
            for line in annotationFile:
                if '#' not in line and 'gene=' + orf in line and line.split()[2] == 'gene':
                    return int(line.split()[3])
            else:
                return -1

    @staticmethod
    def getReferenceSeq(refGenomeFasta):
        for h, s in SimpleFastaParser(open(refGenomeFasta)):
            return s
        return ''

    @staticmethod
    def checkCodon(threeMer):
        try:
            Seq(threeMer).translate()
            return True
        except Exception:
            return False

    def makeProbesByPos(self):
        self.probesByPos = {pos: [] for pos in self.posList + [-1]}
        minPos = min(self.posList)
        maxPos = max(self.posList)
        refSeq = self.getReferenceSeq(self.refGenome)
        if not refSeq:
            self.logUpdate('[ERROR]Failure to get reference sequence from reference genome.')
            self.printUsage()
        for snp in self.snpList:
            self.logUpdate('')
            self.logUpdate('[INFO]SNP {}'.format(snp))
            mutType, orf, mutation = self.parseMutation(snp)
            if mutType not in ['aa', 'nt']:
                self.logUpdate('[ERROR]SNP {} has a wrong format.'.format(snp))
                continue
            if mutType == 'aa':
                if self.refGenomeAnnot == '':
                    self.logUpdate('[ERROR]For Amino Acid based SNPs reference annotation needed.')
                    continue
                orfStartPos = self.getOrfStartPos(self.refGenomeAnnot, orf)
                if orfStartPos == -1:
                    self.logUpdate('[ERROR]Failure to find snp {} in reference annotaion.'.format(snp))
                    continue
                aa1, aa2, mutPos = mutation[0], mutation[-1], int(mutation[1:-1])
                codonStartPos = orfStartPos + (mutPos - 1) * 3 - 1
                codonEndPos = orfStartPos + mutPos * 3 - 1
                refCodon = refSeq[codonStartPos: codonEndPos]
                if aa1 != Seq(refCodon).translate():
                    self.logUpdate('[ERROR]Failure to find SNP {} in reference genome'.format(snp))
                    continue
                seqWithSNP = refSeq[codonStartPos - (maxPos-1): codonEndPos + (self.probLen1 - minPos)]
                strainKmerNearSNP = self.getStrKmerNearSNP(mutation, seqWithSNP)  # blast.fa
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
                    df['locSNP'] = df['diffNT'].apply(lambda x: x + maxPos - 1)
                    df['SNPbyNT'] = df.apply(
                        lambda x: '{}{}{}'.format(x[5][x[6]], codonStartPos + x[6], x[4][x[6]]), axis=1
                    )
                    df['WTsequence'] = seqWithSNP
                except:
                    self.logUpdate('[ERROR]Failure to find snp {} in the strain genome or reference genome'.format(snp))
                    continue
                wtSequence, stSequence, ntSNP, locSnp, found = self.parseBlastResult(blastResult=df)
                self.logUpdate('[INFO]aa:{}:{} converted nt:{}'.format(orf, mutation, ntSNP))
                mutSeqs = ParaSeqs(ntSNP, '{}:{}'.format(orf, mutation), wtSequence, stSequence, mutLoc=locSnp)
            else:

                nt1, nt2, snpPos = mutation[0], mutation[-1], int(mutation[1:-1])
                refNT = refSeq[snpPos]
                if nt1 != refNT:
                    self.logUpdate('[ERROR]Failure to find SNP {} in reference genome'.format(snp))
                    continue
                seqWithSNP = refSeq[snpPos - (maxPos - 1):snpPos + 1 + (self.probLen1 - minPos)]
                blastOutput = self.getStrKmerNearSNP(mutation, seqWithSNP)
                df = pd.read_csv(blastOutput, sep='\t', header=None)
                df.columns = ['subGenome', 'SNPbyNT', 'match', 'STsequence']
                try:
                    df['WTsequence'] = seqWithSNP
                    df['locSNP'] = maxPos - 1
                    df = df[df.STsequence.apply(lambda x: len(x) == len(seqWithSNP))]
                    df = df[df.STsequence.apply(lambda x: x[maxPos - 1]) == nt2]
                except:
                    self.logUpdate('[ERROR]Failure to find snp {} in the strain genome or reference genome'.format(snp))
                    continue
                wtSequence, stSequence, ntSNP, locSnp, found = self.parseBlastResult(blastResult=df)
                mutSeqs = ParaSeqs(ntSNP, '', wtSequence, stSequence, mutLoc=locSnp)
            if found < 0 or not found:
                self.logUpdate('[ERROR]Failure to find SNP {} in strain genome'.format(snp))
                continue
            self.probesByPos[-1].append(mutSeqs)
            for pos in self.posList:
                wtProbe, stProbe = mutSeqs.getProbesWithPos(pos)
                paraSeq = ParaSeqs(mutSeqs.ntSnp, mutSeqs.aaSnp, wtProbe, stProbe, found=found, probLen=self.probLen1)
                self.probesByPos[pos].append(paraSeq)

    def make1stProbe(self):
        probeLines = []
        for pos in self.posList + [-1]:
            probCSV = '{}pos{}.csv'.format(self.probe1byPosDir, pos) if pos>0 else self.probe1byPosDir + 'merged.csv'
            csvWriter = open(probCSV, 'w')
            csvWriter.write('WT sequence,ST sequence,found,ntSNP,aaSNP\n')
            for probs in self.probesByPos[pos]:
                csvWriter.write(
                    '{},{},{},{},{}\n'.format(probs.wtSeq, probs.stSeq, probs.found, probs.ntSnp, probs.aaSnp))
                probeLines.append(
                    '>{}{};{}\n'.format(probs.ntSnp, '=' + probs.aaSnp if probs.aaSnp else '', pos))
                probeLines.append('{}\n'.format(probs.stSeq))
            csvWriter.close()
        self.probe1 = self.workDir + 'probe1.fa'
        with open(self.probe1, 'w') as w:
            w.writelines(probeLines)

    @staticmethod
    def makeProbMaskCoords(inputDF, outputBed):
        inputDF['seqID'] = inputDF['seqID'].apply(lambda x: x.split(';')[0])
        inputDF.to_csv(outputBed, sep='\t', header=False, index=False)
        return outputBed

    def make2ndWindow(self):
        snpPatternKmers = '{}patterns.fasta'.format(self.input2Dir)
        snpMaskingBed = '{}masked.bed'.format(self.input2Dir)
        with open(snpPatternKmers, 'w') as w:
            w.writelines(
                ['>{}{}\n{}\n'.format(p.ntSnp, '=' + p.aaSnp if p.aaSnp else '', p.stSeq) for p in self.probesByPos[-1]]
            )
        self.lookupTSV = ProbeitUtils.getPatternPosition(snpPatternKmers, self.strGenome, self.workDir + 'lookup.tsv')
        snpMaskingBed = self.makeProbMaskCoords(
            pd.read_csv(self.lookupTSV, sep='\t')[['seqID', 'start', 'end']], snpMaskingBed
        )
        # MAKE PROBEs MAKSED FASTA
        self.window2 = ProbeitUtils.getWindowFasta(
            self.strGenome, snpMaskingBed, self.input2Dir + 'masked.fa', self.input2Dir + 'window.bed',
            self.input2Dir + 'seq.fa', self.windowSize
        )
        self.lookup2 = ProbeitUtils.makeLookup(self.window2, self.input2Dir + 'name.lookup')

    def trimProbes(self):
        self.probe2 = self.workDir + 'probe2.fa'
        maskDF = pd.read_csv(self.lookupTSV, sep='\t')
        kmers = list(maskDF['patternName'])
        # PARSING TEMP 2ND FASTA
        w = open(self.probe2, 'w')
        for h, s in SimpleFastaParser(open(self.tempProbe2)):
            kmerIndex = ProbeitUtils.parseGenmapPattern(h)[0]
            soveredSNPs = ':'.join(list((set([kmers[i] for i in kmerIndex]))))
            w.write('>{}\n'.format(soveredSNPs))
            w.write(s + '\n')

    def excute(self):
        self.getPars()
        self.checkArgs()
        self.makeWorkDir()
        self.logUpdate('[INFO]Your arguments: snp ' + ProbeitUtils.getUserArgs(self.args) + '\n')
        self.logUpdate("[INFO]make 1st probes")
        self.makeProbesByPos()
        self.make1stProbe()
        self.logUpdate("[INFO]make 2nd probes")
        self.make2ndWindow()
        ProbeitUtils.delDir(self.indexDir)
        self.tempProbe2 = '{}{}mer.fasta'.format(self.workDir, self.probeLen2)
        self.minzResult2 = self.minz2Dir + 'result'
        ProbeitUtils.makeProbe(
            self.indexDir, self.map2Dir, self.mapError2, self.probeLen2, self.window2, self.window2, self.lookup2,
            self.tempProbe2, self.minzResult2, self.minzResult2 + '.bed', self.minzCovered2, 1, self.minzEarlyStop2,
            self.minzSimScore2, self.minzRepeats2
        )
        self.trimProbes()
        self.logUpdate("[INFO]done")
        ProbeitUtils.sortFasta(self.probe1, self.workDir + 'sorted1.fa')
        ProbeitUtils.sortFasta(self.probe2, self.workDir + 'sorted2.fa')
        return

    @staticmethod
    def printUsage():
        print("Usage")
        print(" -h|--help NONE")
        print('\t Show usage')
        print("python probeit.py snp ")
        print("\t -r [REF GENOME] -s [STR GENOME] -p [positions] -m [SNPs] -o [DIR] -a [REF ANNOTATION]")
        print("REQUIRED OPTIONS")
        print(" -r|--reference FASTA file")
        print("\t The wildtype genome.")
        print(" -s|--strain FASTA file")
        print("\t The strain Genome.")
        print(" -p|--positions COMMA SEPARATED INT ARRAY")
        print("\t Position List: The position of this indicates the position of the SNP on the 1st Probes")
        print(" -m|--mutations COMMA SEPARATED SNP ARRAY")
        print("\t SNP List. Both amino acid differences and nucleotide differences are allowed.")
        print(" -o|--output DIR")
        print("\t Output directory. The Directory is automatically created by Probeit.")
        print(" -a|--annotation GFF file")
        print("\t The wildtype genome annotation. Only required when using amino acid differences in the -m option.")
        print("ADDITIONAL OPTIONS")
        print(" --max-window NONE")
        print("\t When you need maximum window mode, then use this option. Default window mode is minimum window.")
        print(" --window-size INT[200]")
        print("\t size of windows for 2nd probes")
        print(" --probe-len1 INT[40]")
        print("\t Length of 1st Probes")
        print(" --probe-len2 INT[20]")
        print("\t Length of 2nd Probes")
        print(" --probe-error2 INT[1]")
        print("\t The number of error allowed in 2nd Probes")
        print(" --minimizing-covered2 INT[1]")
        print("\t The number of times each 1st Probe should be covered by 2nd Probes")
        print(" --minimizing-repeats2 INT[10]")
        print("\t The number of random iterations when minimizing 2nd Probes")
        print(" --minimizing-earlystop-criteria2 FLOAT[0.99]")
        print("\t Proportion of 1st Probes covered by 2nd Probes for earlystop upon minimizing of 2nd Probe")
        quit()

def main(args):
    probeit = Probeit(args)
    probeit.checkArgs()
    quit()


if __name__ == '__main__':
    main(sys.argv[1:])
