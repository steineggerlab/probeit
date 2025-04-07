#!/usr/bin/env python
from .config import Config
from Bio.SeqIO.FastaIO import SimpleFastaParser
from Bio.Seq import reverse_complement
import os
import getopt
from .Utils import ProbeitUtils, ThermoFilter


class Primer:
    args = []
    shortParams = 'hp:n:o:'
    longParams = [
        # usage
        'help',
        # required
        'positive=', 'negative=', 'output=',
        # optional
        'remove-redundancy', 'not-thermo-filter', 'primer-len=', 'max-amp-len=', 'min-amp-len=', 'threads=',
        'primer-cover=',
        # genmap
        'error=',
        # hidden
        'probe1-cover=', 'probe1-repeat=',
        'dedup-id=', 'rid-neg-id=', 'probe1-dist=',
    ]
    inputGenome = ''
    negGenome = ''
    workDir = ''
    pLen1 = 20
    maxAmpLen = 200
    minAmpLen = 100
    primerCoverage = 3
    cluIdentity = 0.97  # for mmseqs linclust
    ridNegId = 0.90  # for mmseqs easy-search
    cmError1 = 0
    scCoverage1 = 30
    doRemoveRedundancy = False
    doThermoFilter1 = True  # thermo-filter
    isLigationProbe = False  # thermo-filter
    scEarlyStop1 = 0.9
    scScore1 = 11
    scRepeats1 = 1
    threads = 8
    impKmers1 = set()
    probe1Index = []
    lenSeq = {}

    def __init__(self, args):
        args = getopt.getopt(args, self.shortParams, self.longParams)[0]
        print(f"Your arguments: probeit posnegset {ProbeitUtils.getUserArgs(args)}")
        # PARSE PARAMETERS
        for opt, val in args:
            match opt:
                case '-h' | '--help': self.printUsage()
                case '-p' | '--positive': self.inputGenome = str(val)
                case '-n' | '--negative': self.negGenome = str(val)
                case '-o' | '--output': self.workDir = str(val)
                # optional args
                case '--threads': self.threads = int(val)
                case '--remove-redundancy': self.doRemoveRedundancy = True
                case '--not-thermo-filter': self.doThermoFilter1 = False
                case '--primer-len': self.pLen1 = int(val)
                case '--max-amp-len': self.maxAmpLen = int(val)
                case '--min-amp-len': self.minAmpLen = int(val)
                case '--primer-cover': self.primerCoverage = int(val)
                case '--error': self.cmError1 = int(val)
                # HIDDEN args
                case '--probe1-cover': self.scCoverage1 = int(val)
                case '--probe1-repeat': self.scRepeats1 = int(val)
                case '--probe1-earlystop': self.scEarlyStop1 = float(val) / 100
                # identity in mmseqs cluster
                case '--dedup-id': self.cluIdentity = float(val)
                # identity in mmeqs search
                case '--rid-neg-id': self.ridNegId = float(val)
                # setcover similarity
                case '--probe1-dist': self.scScore1 = int(val)
                case _: self.printUsage(opt)

        # FOR DEBUGGING
        print(f'inputGenome: {self.inputGenome}')
        print(f'negGenome: {self.negGenome}')
        print(f'workDir: {self.workDir}')
        print(f'threads: {self.threads}')
        print(f'doRemoveRedundancy: {self.doRemoveRedundancy}')
        print(f'doThermoFilter: {self.doThermoFilter1}')
        print(f'pLen1: {self.pLen1}')
        print(f'cmError1: {self.cmError1}')
        print(f'scCoverage1: {self.scCoverage1}')
        print(f'scRepeats1: {self.scRepeats1}')
        print(f'scEarlyStop1: {self.scEarlyStop1}')
        print(f'cluIdentity: {self.cluIdentity}')
        print(f'ridNegId: {self.ridNegId}')
        print(f'scScore1: {self.scScore1}')

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

        # FILES
        self.log = ProbeitUtils.defineFile(self.workDir, Config.log)
        self.lookup1 = ProbeitUtils.defineFile(self.workDir, 'genome.lookup')
        self.genomeFASTA = self.inputGenome
        self.window1FASTA = ProbeitUtils.defineFile(self.inputDir1, Config.window)
        self.window1PosBED = ProbeitUtils.defineFile(self.inputDir1, 'window.bed')
        self.posKmers1FASTA = ProbeitUtils.defineFile(self.inputDir1, 'kmers.fa')
        self.negRemPosBED = ProbeitUtils.defineFile(self.maskingDir, 'negRemPos.bed')
        self.thermoCalc1TSV = ProbeitUtils.defineFile(self.thermoFiltering1Dir, "thermo_uncut.tsv")
        self.scPosBed1 = ProbeitUtils.defineFile(self.scDir1, 'result.bed')
        self.tempProbe1 = ProbeitUtils.defineFile(self.workDir, 'temp1.fa')
        self.primerFASTA = ProbeitUtils.defineFile(self.workDir, Config.primer)
        # self.rcProbe1 = ProbeitUtils.defineFile(self.workDir, Config.rcProbe1)
        self.logUpdate(
            '[INFO]Your arguments: probeit primer ' + ' '.join(['{} {}'.format(i[0], i[1]).strip() for i in args]),
            False)

    def logUpdate(self, msg, doPrint=True):
        if doPrint:
            print(msg)
        with open(self.log, 'a') as w:
            w.write(msg + '\n')

    def negFilter(self):
        negKmerResult = ProbeitUtils.defineFile(self.maskingDir, 'search.tsv')
        negKmerPosBED = ProbeitUtils.defineFile(self.maskingDir, 'search.bed')
        # REMOVE NEGATIVE K-MERS
        self.logUpdate('[INFO]remove probes found in negative genome')
        self.logUpdate(
            ProbeitUtils.ridNegKmers(self.posKmers1FASTA, self.negGenome, negKmerResult, self.maskingDir, self.ridNegId,
                                     self.threads))

        # MAKE DEDUPLICATED POSITIVE GENOME COORDINATE FILE(BED)
        w = open(self.window1PosBED, 'w')
        idx = 0
        for _, seq in SimpleFastaParser(open(self.window1FASTA)):
            w.write(f'{idx}\t0\t{str(len(seq.strip()))}\n')
            idx += 1
        w.close()

        # EXTRACT NEGATIVE REMOVED POSITIONS BED
        w = open(negKmerPosBED, 'w')
        for line in open(negKmerResult):
            kmers = ProbeitUtils.parseKmers(line)[1:]
            for kmer in kmers:

                w.write(f'{kmer.idx}\t{kmer.sPos}\t{kmer.sPos + self.pLen1}\n')
        w.close()
        self.logUpdate(ProbeitUtils.getSubtractedBed(self.window1PosBED, negKmerPosBED, self.negRemPosBED), False)

    def easyComMap(self, uniqComMap):
        w = open(uniqComMap, 'w')
        self.lenSeq.clear()
        idx = 0
        for h, sPos in SimpleFastaParser(open(self.window1FASTA)):
            self.lenSeq[str(idx)] = len(sPos)
            idx += 1
        negativeKmers = set()
        prevGenome = 0
        prevEnd = 0
        for line in open(self.negRemPosBED):
            genomeIdx, sPos, ePos = map(int, line.strip().split())
            if genomeIdx == prevGenome:
                negativeKmers.update(set([(genomeIdx,i) for i in range(prevEnd, sPos)]))
                prevEnd = ePos
                continue

            negativeKmers.update(set([(prevGenome,i) for i in range(prevEnd, self.lenSeq[prevGenome])]))
            negativeKmers.update(set([(genomeIdx,i) for i in range(0, sPos)]))
            prevGenome = genomeIdx
            prevEnd = ePos

        negativeKmers.update(set([f'{prevGenome},{i}' for i in range(prevEnd, self.lenSeq[prevGenome])]))
        for h, sPos in SimpleFastaParser(open(self.posKmers1FASTA)):
            kmers = ProbeitUtils.parseKmers(h)
            kmerSet = set(kmers)
            isThermoImproper = kmerSet & self.impKmers1
            isUnregularNtFound = set(sPos) != Config.nucleotideSet & set(sPos)
            isNegative = {kmers[0]} & negativeKmers
            if isNegative or isThermoImproper or isUnregularNtFound:
                continue
            w.write(h + '\n')
        w.close()

    def makeWindow1(self):
        ProbeitUtils.sortFasta(self.genomeFASTA)
        ProbeitUtils.simplifyFastaHeaders(self.genomeFASTA, self.window1FASTA)
        ProbeitUtils.makeLookup(self.window1FASTA, self.lookup1, self.window1PosBED)

    def makePrimers(self):
        def kmerParser(rawKmer):
            g, p = map(int, rawKmer.split(','))
            return g, p

        def isCovered(l, key, cov):
            cnt = l.count(key)
            return cnt >= cov

        uniqComMap = ProbeitUtils.defineFile(self.cmDir1, 'uniq.genmap.csv')
        self.easyComMap(uniqComMap)
        msg = ProbeitUtils.makeProbe(
            self.tempProbe1, self.scPosBed1, self.window1FASTA, self.lookup1, self.pLen1,
            self.scCoverage1, self.scEarlyStop1, self.scScore1, self.scRepeats1, uniqComMap,
        )
        self.logUpdate(msg, False)

        # PRIMER
        numGenomes = int([i.split()[0] for i in open(self.lookup1)][-1]) + 1
        genomeAndPairs = {i: [] for i in range(numGenomes)}
        kmers = dict()
        seqs = dict()
        idx = 0
        for h, s in SimpleFastaParser(open(self.tempProbe1)):
            for kmer in h.split()[1].split('|'):
                kmers[kmerParser(kmer)] = idx
                seqs[idx] = s
            idx += 1

        for genomeIdx in range(numGenomes):
            primerCandidates = sorted([(g, p) for g, p in kmers.keys() if g == genomeIdx])
            for primerIdx, primerCandidate1 in enumerate(primerCandidates):
                sPos = primerCandidate1[1]
                for primerCandidate2 in primerCandidates[primerIdx + 1:]:
                    ePos = primerCandidate2[1]
                    ampLen = ePos - sPos + self.pLen1
                    if self.minAmpLen <= ampLen <= self.maxAmpLen:
                        genomeAndPairs[genomeIdx].append((kmers[primerCandidate1], kmers[primerCandidate2]))
        # primers
        with open(self.lookup1) as f:
            genomeKeys = {int(i.split()[0]): i.split()[1] for i in f}

        pairsAndCoveredGenomes = dict()
        for genome, pairs in genomeAndPairs.items():
            for pair in pairs:
                pairsAndCoveredGenomes.setdefault(pair, [])
                pairsAndCoveredGenomes[pair].append(genome)

        pairsCoverage = dict()
        for k, v in pairsAndCoveredGenomes.items():
            pairsCoverage.setdefault(len(v), [])
            pairsCoverage[len(v)].append(k)
        coveredGenomes = []
        primerPairsSet = set()

        for _, pair in sorted(pairsCoverage.items(), reverse=True):
            for g in pairsAndCoveredGenomes[pair]:
                if not isCovered(coveredGenomes, g, self.primerCoverage):
                    coveredGenomes.append(g)
                    primerPairsSet.add(pair)

                if len(coveredGenomes) < numGenomes:
                    continue

                break

        idx = 0
        w = open(self.primerFASTA, 'w')
        for i, j in primerPairsSet:
            genomes = [genome for genome, primerList in genomeAndPairs.items() if (i, j) in primerList]
            genomes = [genomeKeys[i] for i in genomes]
            genomes = ';'.join(genomes)
            w.write(f'>primer_pair{idx}\t{genomes}\n{seqs[i]}\n{reverse_complement(seqs[j])}\n')
            idx += 1
        w.close()

    def run(self):
        # MAKE PROBE1
        self.logUpdate("[INFO] make 1st probes")

        # REMOVE REDUNDANCY FROM INPUT GENOME
        if self.doRemoveRedundancy:
            self.logUpdate('[INFO]deduplicate positive fasta')
            clusteredGenome = ProbeitUtils.defineFile(self.inputDir1, 'dedup')
            tempDir = ProbeitUtils.defineDirectory('temp', make=False, root=self.inputDir1)
            msg, self.genomeFASTA = ProbeitUtils.clusterGenome(self.inputGenome, clusteredGenome, tempDir,
                                                               self.cluIdentity, threads=self.threads)
            self.logUpdate(msg)

        # MAKE WINDOW and KMERS FOR PROBE1
        self.makeWindow1()
        ProbeitUtils.extractKmers(self.window1FASTA, self.posKmers1FASTA, self.pLen1)

        # DO NEG FILTER FOR PROBE1
        if self.negGenome:
            self.negFilter()

        # DO THERMO FILTER FOR PROBE1
        if self.doThermoFilter1:
            thermoFilter1 = ThermoFilter(self.posKmers1FASTA, self.pLen1, self.thermoCalc1TSV, self.isLigationProbe)
            msg, self.impKmers1 = thermoFilter1.run()
            self.logUpdate(msg)

        # COMPLETE MAKE PROBE1
        self.makePrimers()

        print('COMPLETE!!!')
        return

    @staticmethod
    def printUsage(wrong=''):
        if wrong:
            print(f"[ERR] Unknown parameter found: {wrong}")
        print("Probeit primer")
        print("It generates primer pairs with sequences included in the positive genome but not in the negative genome")
        print("probeit primer -p [POSITIVE GENOME]-n [NEGATIVE GENOME] -o [DIR]")

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
        print(" --remove-redundancy NONE")
        print("\t Use it when you NEED to cluster positive genome")
        print(" --not-thermo-filter NONE")
        print("\t Use it when you DO NOT need the thermodynamic filter")
        print(" --primer-len INT[20]")
        print("\t Length of 1st Probes")
        print(" --max-amp-len INT[200]")
        print("\t Maximum length of primers")
        print(" --min-amp-len INT[100]")
        print("\t Minimum length of primers")
        print("--primer-cover INT[3]")
        print("\t The number of times each Seqs from positive genome should be covered by primer pairs")
        print("ADDITIONAL OPTIONS FOR GENMAP PART: Genmap calculates mappability by summarizing k-mers.")
        print(" --error INT[0]")
        print("\t The number of error allowed in primers")

        # print("ADDITIONAL OPTIONS FOR SETCOVER PART: Setcover makes probe sets cover targets with the minimum number of probes.")
        # print("--probe1-cover INT[1]")
        # print("\t The number of times each Seqs from positive genome should be covered by 1st Probes")
        # print("--probe1-repeat INT[1]")
        # print("\t The number of random iterations when minimizing 1st Probes")
        # print("--probe1-earlystop INT[90]")
        # print("\t Early stop picking new probes if X% of sequences are covered at least N(--probe1-cover) times")

        quit()