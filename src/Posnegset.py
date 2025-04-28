#!/usr/bin/env python
from .config import Config
from Bio.SeqIO.FastaIO import SimpleFastaParser
from Bio.Seq import Seq
import os
import getopt
from .Utils import ProbeitUtils, ThermoFilter, Kmer, BedLine


class PosNegSet:
    args = []
    shortParams = 'hp:n:o:'
    longParams = [
        # usage
        'help',
        # required
        'positive=', 'negative=', 'output=',
        # optional
        'not-make-probe2', 'remove-redundancy', 'not-thermo-filter', 'ligation-probe', 'probe1-len=', 'probe2-len=',
        'window-size=', 'threads=',
        # genmap
        'probe1-error=', 'probe2-error=',
        # setcover
        'probe1-cover=', 'probe2-cover=', 'probe1-repeat=', 'probe2-repeat=', 'probe1-earlystop=', 'probe2-earlystop=',
        # hidden
        'dedup-id=', 'rid-neg-id=', 'probe1-dist=', 'probe2-dist='
    ]
    # CLASS VARs
    inputGenome = ''
    negGenome = ''
    workDir = ''
    pLen1 = 40
    pLen2 = 20
    cluIdentity = 0.97  # for mmseqs linclust
    ridNegId = 0.90  # for mmseqs easy-search
    cmError1 = 0
    cmError2 = 1
    scCoverage1 = 1
    scCoverage2 = 1
    doRemoveRedundancy = False
    needProbe2 = True
    doThermoFilter1 = True  # thermo-filter
    doThermoFilter2 = True
    isLigationProbe = False  # thermo-filter
    windowSize = 200
    scEarlyStop1 = 0.9
    scScore1 = 11
    scRepeats1 = 1
    scEarlyStop2 = 0.99
    scScore2 = 20
    scRepeats2 = 10
    threads = 8

    def __init__(self, args):
        self.cluFile2 = None
        self.cluFile = None
        self.impKmers1 = set()
        self.impKmers2 = set()
        self.probe1Index = []
        self.lookupDict = {}
        self.negativeKmers = set()
        args = getopt.getopt(args, self.shortParams, self.longParams)[0]
        print(f"Your arguments: probeit posnegset {ProbeitUtils.getUserArgs(args)}")
        # PARSE PARAMETERS
        for opt, val in args:
            match opt:
                case '-h' | '--help': self.printUsage()
                case '-p' | '--positive': self.inputGenome = str(val)
                case '-n' | '--negative': self.negGenome = str(val)
                case '-o' | '--output': self.workDir = str(val)
                case '--threads': self.threads = int(val)
                case '--window-size': self.windowSize = int(val)
                case  '--remove-redundancy': self.doRemoveRedundancy = True
                case '--not-make-probe2':  self.needProbe2 = False
                case '--not-thermo-filter': self.doThermoFilter1 = False; self.doThermoFilter2 = False
                case '--ligation-probe': self.isLigationProbe = True
                case '--probe1-len': self.pLen1 = int(val)
                case '--probe2-len': self.pLen2 = int(val)
                case '--probe1-error': self.cmError1 = int(val)
                case '--probe2-error':  self.cmError2 = int(val)
                case '--probe1-cover': self.scCoverage1 = int(val)
                case '--probe2-cover': self.scCoverage2 = int(val)
                case '--probe1-repeat': self.scRepeats1 = int(val)
                case '--probe2-repeat': self.scRepeats2 = int(val)
                case '--probe1-earlystop':  self.scEarlyStop1 = float(val) / 100
                case '--probe2-earlystop': self.scEarlyStop2 = float(val) / 100
                # HIDDEN args
                # identity in mmseqs cluster
                case '--dedup-id': self.cluIdentity = float(val)
                # identity in mmeqs search
                case '--rid-neg-id': self.ridNegId = float(val)
                # setcover similarity
                case '--probe1-dist': self.scScore1 = int(val)
                case '--probe2-dist': self.scScore2 = int(val)
                case _: self.printUsage(opt)

        # FOR DEBUGGING
        print(f'inputGenome: {self.inputGenome}')
        print(f'negGenome: {self.negGenome}')
        print(f'workDir: {self.workDir}')
        print(f'threads: {self.threads}')
        print(f'windowSize: {self.windowSize}')
        print(f'doRemoveRedundancy: {self.doRemoveRedundancy}')
        print(f'doDesignProbe2: {self.needProbe2}')
        print(f'doThermoFilter: {self.doThermoFilter1}')
        print(f'isLigationProbe: {self.isLigationProbe}')
        print(f'pLen1: {self.pLen1}')
        print(f'pLen2: {self.pLen2}')
        print(f'cmError1: {self.cmError1}')
        print(f'cmError2: {self.cmError2}')
        print(f'scCoverage1: {self.scCoverage1}')
        print(f'scCoverage2: {self.scCoverage2}')
        print(f'scRepeats1: {self.scRepeats1}')
        print(f'scRepeats2: {self.scRepeats2}')
        print(f'scEarlyStop1: {self.scEarlyStop1}')
        print(f'scEarlyStop2: {self.scEarlyStop2}')
        print(f'cluIdentity: {self.cluIdentity}')
        print(f'ridNegId: {self.ridNegId}')
        print(f'scScore1: {self.scScore1}')
        print(f'scScore2: {self.scScore2}')

        # CHECK USER"S ARGUMENTS
        message = "{}You didn't input a proper argument for '{}' parameter or missed it."
        isBadArguments = False

        if not self.inputGenome:
            print(message.format('[ERR]', '--positive'))
            isBadArguments = True

        if not self.negGenome:
            print('[INFO] Probeit posnegset without negative genomes.')

        if not self.workDir:
            print(message.format('[ERR]', '--output'))
            isBadArguments = True

        if os.path.exists(self.workDir):
            print(f'[ERR] The directory named {self.workDir} already exists.')
            isBadArguments = True

        if isBadArguments:
            self.printUsage()

        print('You input proper arguments.')

        # DIRECTORIES
        self.workDir = ProbeitUtils.defineDirectory(self.workDir.split(os.path.sep)[0], make=True)
        self.inputDir1 = ProbeitUtils.defineDirectory('input1', make=True, root=self.workDir)
        self.maskingDir = ProbeitUtils.defineDirectory('neg_filter', make=True, root=self.workDir)
        self.thermoFiltering1Dir = ProbeitUtils.defineDirectory('thermo_filter1', make=True, root=self.workDir)
        self.cmDir1 = ProbeitUtils.defineDirectory('mapping_probe1', make=True, root=self.workDir)
        self.idxDir1 = ProbeitUtils.defineDirectory('index_probe1', make=False, root=self.workDir)
        self.scDir1 = ProbeitUtils.defineDirectory('setcover_probe1', make=True, root=self.workDir)

        if self.needProbe2:
            self.inputDir2 = ProbeitUtils.defineDirectory('input2', make=True, root=self.workDir)
            self.thermoFiltering2Dir = ProbeitUtils.defineDirectory('thermo_filter2', make=True, root=self.workDir)
            self.cmDir2 = ProbeitUtils.defineDirectory('mapping_probe2', make=True, root=self.workDir)
            self.idxDir2 = ProbeitUtils.defineDirectory('index_probe2', make=False, root=self.workDir)
            self.scDir2 = ProbeitUtils.defineDirectory('setcover_probe2', make=True, root=self.workDir)

        # FILES
        self.log = ProbeitUtils.defineFile(self.workDir, Config.log)
        self.lookup1 = ProbeitUtils.defineFile(self.workDir, 'genome.lookup')
        self.clusteredGenome = self.inputGenome
        self.window1FASTA = ProbeitUtils.defineFile(self.inputDir1, Config.window)
        self.window1PosBED = ProbeitUtils.defineFile(self.inputDir1, 'window.bed')
        self.posKmers1FASTA = ProbeitUtils.defineFile(self.inputDir1, 'kmers.fa')
        self.negRemPosBED = ProbeitUtils.defineFile(self.maskingDir, 'negRemPos.bed')
        self.thermoCalc1TSV = ProbeitUtils.defineFile(self.thermoFiltering1Dir, "thermo_uncut.tsv")
        self.scPosBed1 = ProbeitUtils.defineFile(self.scDir1, 'result.bed')
        self.tempProbe1 = ProbeitUtils.defineFile(self.workDir, 'temp1.fa')
        self.probe1 = ProbeitUtils.defineFile(self.workDir, Config.probe1)
        self.rcProbe1 = ProbeitUtils.defineFile(self.workDir, Config.rcprobe1)

        if self.needProbe2:
            self.genome2 = ProbeitUtils.defineFile(self.inputDir2, "genome.fa")
            self.window2FASTA = self.genome2 #ProbeitUtils.defineFile(self.inputDir2, Config.window)
            self.window2PosBED = ProbeitUtils.defineFile(self.inputDir2, 'window.bed')
            self.lookup2 = ProbeitUtils.defineFile(self.workDir, 'probe1.lookup')
            self.kmers2FASTA = ProbeitUtils.defineFile(self.inputDir2, 'kmers.fa')
            self.thermoCalc2TSV = ProbeitUtils.defineFile(self.thermoFiltering2Dir, 'thermo_uncut.tsv')
            self.scPosBed2 = ProbeitUtils.defineFile(self.scDir2, 'result.bed')
            self.tempProbe2 = ProbeitUtils.defineFile(self.workDir, 'temp2.fa')
            self.probe2 = ProbeitUtils.defineFile(self.workDir, Config.probe2)
            self.rcProbe2 = ProbeitUtils.defineFile(self.workDir, Config.rcprobe2)

        self.logUpdate(
            '[INFO]Your arguments: probeit posnegset ' + ' '.join(['{} {}'.format(i[0], i[1]).strip() for i in args]),
            False)

    def logUpdate(self, msg, needPrint=True):
        if needPrint:
            print(msg)

        with open(self.log, 'a') as w:
            w.write(msg + '\n')

    def negFilter(self, pLen):
        negKmerResult = ProbeitUtils.defineFile(self.maskingDir, Config.negKmers)
        negKmerPosBED = ProbeitUtils.defineFile(self.maskingDir, Config.negBed)
        # REMOVE NEGATIVE K-MERS
        self.logUpdate('[INFO]remove probes found in negative genome')
        self.logUpdate(ProbeitUtils.ridNegKmers(self.posKmers1FASTA, self.negGenome, negKmerResult, self.maskingDir, self.ridNegId, self.threads))

        # EXTRACT NEGATIVE REMOVED POSITIONS BED
        w = open(negKmerPosBED, 'w')
        for i in open(negKmerResult):
            kmers = ProbeitUtils.parseKmers(i)[1:]
            for kmer in kmers:
                w.write(BedLine(kmer.idx, kmer.sPos, kmer.sPos + pLen).getLine())

        w.close()
        self.logUpdate(ProbeitUtils.getSubtractedBed(self.window1PosBED, negKmerPosBED, self.negRemPosBED), False)

    def getNegativeKmers(self):
        genomeLength = {name: 0 for name in self.lookupDict.values()}
        for h, s in SimpleFastaParser(open(self.window1FASTA)):
            genomeLength[self.lookupDict[h]] = len(s)

        negativeKmers = set()
        if not self.negGenome:
            return negativeKmers

        prevGenomeIdx = 0
        prevEnd = 0
        for line in open(self.negRemPosBED):
            genomeIdx, currStart, currEnd = map(int, line.strip().split())
            if genomeIdx == prevGenomeIdx:
                negativeKmers.update(set([Kmer(idx=genomeIdx, sPos=i) for i in range(prevEnd, currStart)]))
                prevEnd = currEnd
                continue

            negativeKmers.update(set([Kmer(idx=prevGenomeIdx,sPos=i) for i in range(prevEnd, genomeLength[prevGenomeIdx])]))
            negativeKmers.update(set([Kmer(idx=genomeIdx, sPos=i) for i in range(currStart)]))
            prevGenomeIdx = genomeIdx
            prevEnd = currEnd

        negativeKmers.update(set([Kmer(idx=prevGenomeIdx,sPos=i) for i in range(prevEnd, genomeLength[prevGenomeIdx])]))
        return negativeKmers

    def getMaskedGenome(self):
        def _getWindowSeq(kmer, idx1, idx2):
            gIdx = kmer.idx
            sPos = kmer.sPos
            ePos = sPos + self.pLen1
            header = window1Headers[gIdx]
            seq = window1Seqs[gIdx]
            seqLen = len(seq)
            windowStartPos = 0 if sPos - self.windowSize < 0 else sPos - self.windowSize
            windowEndPos = seqLen if ePos + self.windowSize > seqLen else ePos + self.windowSize
            seq = seq[windowStartPos: windowEndPos]
            return f'>probe_{idx1}_{idx2}:{header}:{sPos}\n{seq}\n'

        window1Headers = []
        window1Seqs = []
        for h, s in SimpleFastaParser(open(self.inputGenome)):
            window1Headers.append(h.split()[0])
            window1Seqs.append(s)

        probe1Kmers = [ProbeitUtils.parseKmers(line) for line in open(self.scPosBed1)]
        for kmers in probe1Kmers:
            for kmer in kmers:
                ePos = kmer.sPos + self.pLen1
                window1Seqs[kmer.idx] = window1Seqs[kmer.idx][:kmer.sPos] + 'N' * self.pLen1 + window1Seqs[kmer.idx][ePos:]

        with open(self.genome2, 'w') as w:
            w.writelines([_getWindowSeq(kmer, i, j) for i, kmers in enumerate(probe1Kmers) for j, kmer in enumerate(kmers)])

    def makeProbe1(self):
        genomeKeys = {idx:name for name, idx in self.lookupDict.items()}
        probe1Writer = open(self.probe1, 'w')
        rcprobe1Writer = open(self.rcProbe1, 'w')
        probeIdx = 0
        for seq, kmers in {s: ProbeitUtils.parseKmers(h) for h, s in SimpleFastaParser(open(self.tempProbe1))}.items():
            rcSeq = Seq(seq).reverse_complement()
            parsedKmers = []

            for kmerIdx, kmer in enumerate(kmers):
                genomeName = genomeKeys[kmer.idx]
                parsedKmers.append(f'probe_{probeIdx}_{kmerIdx}:{genomeName}:{kmer.sPos}')

            self.probe1Index += [kmer.split(':')[0] for kmer in parsedKmers]
            probe1Writer.write(f">probe_{probeIdx}\t{'|'.join(parsedKmers)}\n{seq}\n")
            rcprobe1Writer.write(f">probe_{probeIdx}\t{'|'.join(parsedKmers)}\n{rcSeq}\n")
            probeIdx += 1

        probe1Writer.close()
        rcprobe1Writer.close()

    def makeProbe2(self):
        probe2Writer = open(self.probe2, 'w')
        rcProbe2Writer = open(self.rcProbe2, 'w')
        probeIdx = 0
        for h, seq in sorted([(h, s) for h, s in SimpleFastaParser(open(self.tempProbe2))]):
            rcSeq = Seq(seq).reverse_complement()
            keys, _ = ProbeitUtils.parseGenmapPattern(h)
            keys = [self.probe1Index[k] for k in sorted(keys)]
            probe2Writer.write('>cap_{}\t|{}\n{}\n'.format(probeIdx, '|'.join(keys), seq))
            rcProbe2Writer.write('>cap_{}\t|{}\n{}\n'.format(probeIdx, '|'.join(keys), rcSeq))
            probeIdx += 1

        probe2Writer.close()
        rcProbe2Writer.close()

    def run(self):
        # self.genome : input
        # self.genomeFASTA : clustered input
        # MAKE PROBE1
        self.logUpdate("[INFO] make 1st probes")

        # REMOVE REDUNDANCY FROM INPUT GENOME
        if self.doRemoveRedundancy:
            self.logUpdate('[INFO]deduplicate positive fasta')
            tempDir = ProbeitUtils.defineDirectory('temp', make=False, root=self.inputDir1)
            msg, self.clusteredGenome, self.cluFile = ProbeitUtils.clusterGenome(self.inputGenome, ProbeitUtils.defineFile(self.inputDir1, 'dedup'), tempDir, self.cluIdentity, threads=self.threads)
            self.logUpdate(msg)

        # MAKE WINDOW and KMERS FOR PROBE1
        ProbeitUtils.simplifyFastaHeaders(self.clusteredGenome, self.window1FASTA)
        self.lookupDict = ProbeitUtils.makeLookup(self.inputGenome, self.lookup1, True)
        ProbeitUtils.makeGenomePos(self.window1FASTA, self.window1PosBED)
        if self.negGenome or self.doThermoFilter1:
            ProbeitUtils.extractKmers(self.window1FASTA, self.posKmers1FASTA, self.pLen1)

        # DO NEG FILTER FOR PROBE1
        if self.negGenome:
            self.negFilter(self.pLen1)
            self.negativeKmers = self.getNegativeKmers()

        # DO THERMO FILTER FOR PROBE1
        if self.doThermoFilter1:
            thermoFilter1 = ThermoFilter(self.posKmers1FASTA, self.pLen1, self.thermoCalc1TSV, self.isLigationProbe)
            msg, self.impKmers1 = thermoFilter1.run()
            self.logUpdate(msg)

        # Compute Mapping
        uniqComMap = ProbeitUtils.defineFile(self.cmDir1, Config.uniqMapping)
        match self.cmError1:
            case 0:
                ProbeitUtils.simpleComputeMappability(self.window1FASTA, self.window1PosBED, uniqComMap, self.pLen1, improperKmers=self.impKmers1 | self.negativeKmers)
            case _:
                ProbeitUtils.computeMappability(self.window1FASTA, self.idxDir1, self.cmError1, self.pLen1, self.cmDir1, uniqComMap, threads=self.threads, improperKmers=self.impKmers1 | self.negativeKmers)

        if self.doRemoveRedundancy:
            ProbeitUtils.expandMapResult(uniqComMap, self.window1FASTA, self.inputGenome, self.cluFile, self.lookupDict, self.pLen1)

        msg = ProbeitUtils.makeProbe(self.tempProbe1, self.scPosBed1, self.inputGenome, self.lookup1, self.pLen1, self.scCoverage1, self.scEarlyStop1, self.scScore1,  self.scRepeats1, uniqComMap)
        self.logUpdate(msg, False)

        # COMPLETE MAKE PROBE1
        self.makeProbe1()

        if not self.needProbe2:
            print('COMPLETE!!!')
            return

        # MAKE PROBE2
        self.logUpdate("[INFO] make 2nd probes")

        # MAKE WINDOW FOR PROBE2
        self.getMaskedGenome()
        self.lookupDict = ProbeitUtils.makeLookup(self.genome2, self.lookup2, True)

        # REMOVE REDUNDANCY FROM MASKED GENOME
        if self.doRemoveRedundancy:
            msg, self.window2FASTA, self.cluFile2 = ProbeitUtils.clusterGenome(self.genome2, ProbeitUtils.defineFile(self.inputDir2, 'dedup'), ProbeitUtils.defineDirectory('temp', make=False, root=self.inputDir2), self.cluIdentity, threads=self.threads)
            self.logUpdate(msg)

        ProbeitUtils.makeGenomePos(self.window2FASTA, self.window2PosBED)

        # THERMO FILTER FOR PROBE2
        if self.doThermoFilter2:
            ProbeitUtils.extractKmers(self.window2FASTA, self.kmers2FASTA, self.pLen2)
            thermoFilter2 = ThermoFilter(self.kmers2FASTA, self.pLen2, self.thermoCalc2TSV)
            msg, self.impKmers2 = thermoFilter2.run()
            self.logUpdate(msg)

        # Compute Mapping
        uniqComMap = ProbeitUtils.defineFile(self.cmDir2, Config.uniqMapping)
        match self.cmError2:
            case 0:
                ProbeitUtils.simpleComputeMappability(self.window2FASTA, self.window2PosBED, uniqComMap, self.pLen2, improperKmers=self.impKmers2)
            case _:
                ProbeitUtils.computeMappability(self.window2FASTA, self.idxDir2, self.cmError2, self.pLen2, self.cmDir2, uniqComMap, threads=self.threads, improperKmers=self.impKmers2)

        if self.doRemoveRedundancy:
            ProbeitUtils.expandMapResult(uniqComMap, self.window2FASTA, self.genome2, self.cluFile2, self.lookupDict, self.pLen2)

        msg = ProbeitUtils.makeProbe(self.tempProbe2, self.scPosBed2, self.genome2, self.lookup2, self.pLen2, self.scCoverage2, self.scEarlyStop2, self.scScore2, self.scRepeats2, uniqComMap, overlap=True)
        self.logUpdate(msg, False)

        # COMPLETE MAKE PROBE2
        self.makeProbe2()
        print('COMPLETE!!!')
        return

    @staticmethod
    def printUsage(wrong=''):
        if wrong:
            print(f"[ERR] Unknown parameter found: {wrong}")
        print("Probeit posnegset")
        print("It generates a probe set with sequences included in the positive genome but not in the negative genome")
        print("probeit -p [POSITIVE GENOME]-n [NEGATIVE GENOME] -o [DIR]")

        print("Usage")
        print(" -h|--help NONE")
        print("\t Show usage")

        print("REQUIRED OPTIONS")
        print(" -p|--positive FASTA file")
        print("\t The genome which MUST be covered by the probes.")
        print(" -o|--output DIR")
        print("\t Output directory The Directory is automatically created by Probeit.")

        print("ADDITIONAL OPTIONS")
        print(" -n|--negative FASTA file")
        print("\t The genome which MUST NOT be covered by the probes.")
        print(" --threads INT[8]")
        print("\t number of CPU-cores used")
        print(" --window-size INT[200]")
        print("\t size of windows for 2nd probes")
        print(" --remove-redundancy NONE")
        print("\t Use it when you NEED to cluster positive genome")
        print(" --not-make-probe2 NONE")
        print("\t Use it when you DO NOT need to make 2nd probes")
        print(" --not-thermo-filter NONE")
        print("\t Use it when you DO NOT need the thermodynamic filter")
        print(" --ligation-probe NONE")
        print("\t Use it when you want to make ligation probes as probe1.")
        print(" --probe1-len INT[40]")
        print("\t Length of 1st Probes")
        print(" --probe2-len INT[20]")
        print("\t Length of 2nd Probes")

        print("ADDITIONAL OPTIONS FOR GENMAP PART: Genmap calculates mappability by summarizing k-mers.")
        print(" --probe1-error INT[0]")
        print("\t The number of error allowed in 1st Probes")
        print(" --probe2-error INT[1]")
        print("\t The number of error allowed in 2nd Probes")

        print("ADDITIONAL OPTIONS FOR SETCOVER PART: Setcover makes probe sets cover targets with the minimum number of probes.")
        print("--probe1-cover INT[1]")
        print("\t The number of times each Seqs from positive genome should be covered by 1st Probes")
        print("--probe2-cover INT[1]")
        print("\t The number of times each 1st Probe should be covered by 2nd Probes")
        print("--probe1-repeat INT[1]")
        print("\t The number of random iterations when minimizing 1st Probes")
        print("--probe2-repeat INT[10]")
        print("\t The number of random iterations when minimizing 2nd Probes")
        print("--probe1-earlystop INT[90]")
        print("\t Early stop picking new probes if X% of sequences are covered at least N(--probe1-cover) times")
        print("--probe2-earlystop INT[99]")
        print("\t Early stop picking new probes if X% of sequences are covered at least N(--probe2-cover) times")
        quit()
