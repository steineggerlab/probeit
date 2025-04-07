#!/usr/bin/env python
from .config import Config
from Bio.SeqIO.FastaIO import SimpleFastaParser
from Bio.Seq import Seq
import pandas as pd
import os
import getopt
from .Utils import ParaSeqs, ProbeitUtils, ThermoFilter

class SNP:
    args = []
    shortParams = 'hr:a:s:p:m:o:'
    longParams = [
        # usage
        'help',
        # required
        'reference=', 'strain=', 'positions=', 'mutations=', 'output=', 'annotation=',
        # optional
        'not-make-probe2', 'not-thermo-filter', 'max-window', 'window-size=', 'probe1-len=', 'probe2-len=', 'threads=',
        # genmap
        'probe-error2=',
        # setcover
        'probe2-cover=', 'probe2-earlystop=', 'probe2-repeat=',
        # hidden
        'probe2-dist=', 'search-kmer='
    ]
    refGenome = ''
    refGenomeAnnot = ''
    strGenome = ''
    posList = []
    snpList = []
    workDir = ''
    needProbe2 = True
    isMaxWindow = False
    windowSize = 200
    pLen1 = 40
    pLen2 = 20
    scCoverage2 = 1
    scEarlyStop2 = 0.99
    scScore2 = 20
    scRepeats2 = 10
    cmError2 = 1
    searchKmer = 12
    threads = 8
    doThermoFilter2 = True
    impKmers2 = set()

    def __init__(self, args):
        args = getopt.getopt(args, self.shortParams, self.longParams)[0]
        print(f"Your arguments: probeit snp {ProbeitUtils.getUserArgs(args)}")
        # PARSE PARAMETERS
        # isReference = True
        for opt, val in args:
            match opt:
                case '-h' | '--help': self.printUsage()
                # required
                case '-r' | '--reference': self.refGenome = val
                case '-s' | '--strain': self.strGenome = val
                case '-o' | '--output': self.workDir = val
                case '-p' | '--positions':   self.posList = ProbeitUtils.getArgList(val, isInt=True)
                case '-m' | '--mutations': self.snpList = ProbeitUtils.getArgList(val)
                case '-a' | '--annotation': self.refGenomeAnnot = val
                # optional
                case '--not-make-probe2': self.needProbe2 = False
                case '--not-thermo-filter': self.doThermoFilter2 = False
                case '--threads': self.threads = int(val)
                case '--max-window': self.isMaxWindow = True
                case '--window-size': self.windowSize = int(val)
                case '--probe1-len': self.pLen1 = int(val)
                case '--probe2-len': self.pLen2 = int(val)
                case '--probe2-error': self.cmError2 = int(val)
                case '--probe2-cover': self.scCoverage2 = int(val)
                case '--probe2-earlystop': self.scEarlyStop2 = float(val) / 100
                case '--probe2-repeat': self.scRepeats2 = int(val)
                # hidden args
                case '--probe2-dist': self.scScore2 = int(val)
                case '--search-kmer': self.searchKmer = int(val)
                case _: self.printUsage(opt)

        self.snpList = sorted(list(set(self.snpList)))

        # FOR DEBUGGING
        print(f'refGenome: {self.refGenome}')
        print(f'strGenome: {self.strGenome}')
        print(f'workDir: {self.workDir}')
        print(f'posList: {self.posList}')
        print(f'snpList: {self.snpList}')
        print(f'refGenomeAnnot: {self.refGenomeAnnot}')
        print(f'needProbe2: {self.needProbe2}')
        print(f'threads: {self.threads}')
        print(f'isMaxWindow: {self.isMaxWindow}')
        print(f'windowSize: {self.windowSize}')
        print(f'probLen1: {self.pLen1}')
        print(f'probeLen2: {self.pLen2}')
        print(f'cmError2: {self.cmError2}')
        print(f'scCoverage2: {self.scCoverage2}')
        print(f'scEarlyStop2: {self.scEarlyStop2}')
        print(f'scRepeats2: {self.scRepeats2}')
        print(f'scScore2: {self.scScore2}')
        print(f'searchKmer: {self.searchKmer}')

        validPosList = []
        for p in self.posList:
            if 0 < p <= self.pLen1:
                validPosList.append(int(p))
                continue
            print('[ERROR] {} is not proper for position list.')
        self.posList = validPosList

        # print("Your arguments: {}".format('snp ' + ProbeitUtils.getUserArgs(args)))
        # CHECK USER"S ARGUMENTS
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
            print(message.format('[ERROR]', '--mutations'))
            isBadArguments = True
        if not self.workDir:
            print(message.format('[ERROR]', '--output'))
            isBadArguments = True

        if os.path.exists(self.workDir):
            print(f'[ERR] The directory named {self.workDir} already exists.')
            isBadArguments = True

        if isBadArguments:
            self.printUsage()
        print('You input proper arguments.' if self.refGenomeAnnot else message.format('[WARN]', '--annotation'))

        if not self.isMaxWindow:
            self.windowSize = self.windowSize - (max(self.posList) - min(self.posList))
        # DIRECTORIES
        self.workDir = ProbeitUtils.defineDirectory(self.workDir.split(os.path.sep)[0], make=True)
        self.searchDir = ProbeitUtils.defineDirectory('search', make=True, root=self.workDir)
        self.inputDir1 = ProbeitUtils.defineDirectory('input1', make=True, root=self.workDir)
        if self.needProbe2:
            self.inputDir2 = ProbeitUtils.defineDirectory('input2', make=True, root=self.workDir)
            self.thermoFilterDir2 = ProbeitUtils.defineDirectory('thermo_filter2', make=True, root=self.workDir)
            self.idxDir2 = ProbeitUtils.defineDirectory('index', make=False, root=self.workDir)
            self.cmDir2 = ProbeitUtils.defineDirectory('mapping2', make=True, root=self.workDir)
            self.scDir2 = ProbeitUtils.defineDirectory('setcover2', make=True, root=self.workDir)
        # FILES
        self.log = ProbeitUtils.defineFile(self.workDir, 'log.txt')
        self.posProbeCSVs = {pos: ProbeitUtils.defineFile(self.inputDir1, f'pos{pos}.csv') for pos in self.posList}
        self.posProbeCSVs[-1] = ProbeitUtils.defineFile(self.inputDir1, 'merged.csv')
        self.probe1 = ProbeitUtils.defineFile(self.workDir, Config.probe1)
        self.rcProbe1 = ProbeitUtils.defineFile(self.workDir, Config.rcprobe1)
        if self.needProbe2:
            self.window1FASTA = ProbeitUtils.defineFile(self.inputDir1, Config.window)
            self.lookup = ProbeitUtils.defineFile(self.workDir, 'probe1.lookup')
            self.window1PosBED = ProbeitUtils.defineFile(self.inputDir1, 'window1.bed')
            self.window1TSV = ProbeitUtils.defineFile(self.inputDir1, 'window1.tsv')
            self.maskedGenomeFASTA = ProbeitUtils.defineFile(self.inputDir2, 'masked.fa')
            self.window2PosBED = ProbeitUtils.defineFile(self.inputDir2, 'window2.bed')
            self.window2FASTA = ProbeitUtils.defineFile(self.inputDir2, Config.window)
            self.kmers2FASTA = ProbeitUtils.defineFile(self.thermoFilterDir2, 'kmers.fa')
            self.thermoCalcTSV2 = ProbeitUtils.defineFile(self.thermoFilterDir2, 'thermo_uncut.tsv')
            self.tempProbe2 = ProbeitUtils.defineFile(self.workDir, 'temp2.fa')
            self.probe2 = ProbeitUtils.defineFile(self.workDir, Config.probe2)
            self.rcprobe2 = ProbeitUtils.defineFile(self.workDir, Config.rcprobe2)
        self.logUpdate('[INFO]Your arguments: probeit snp ' + ProbeitUtils.getUserArgs(args), False)
        # VARIABLES
        self.probesByPos = {pos: [] for pos in self.posList + [-1]}

    def logUpdate(self, msg, doPrint=True):
        if doPrint:
            print(msg)

        with open(self.log, 'a') as w:
            w.write(msg + '\n')

    def _searchSnpFromStrGenome(self, mutation, seqWithSNP):
        snpKmers = ProbeitUtils.defineFile(self.searchDir, f'search_{mutation}.tsv')
        searchProbe = ProbeitUtils.defineFile(self.searchDir, f'search_{mutation}.fa')
        with open(searchProbe, 'w') as w:
            w.write('>{}\n{}\n'.format(mutation, seqWithSNP))
        self.logUpdate(ProbeitUtils.searchSNPs(self.searchDir, searchProbe, self.strGenome, snpKmers, self.searchKmer, threads=self.threads))
        return snpKmers

    @staticmethod
    def _parseSearchResult(searchResult):
        #TODO
        found = -1
        groupedDf = searchResult.groupby(['WTsequence', 'STsequence', 'locSNP', 'SNPbyNT'])
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
        n = len(parsedMutation)
        mutType = parsedMutation[0]
        match mutType:
            case 'nt':
               if n != 2:
                   return None, None, None
               snp = parsedMutation[1]
               orf = None

            case 'aa':
                if n != 3:
                    return None, None, None
                _, orf, snp = tuple(parsedMutation)
                orf = orf.replace('orf', 'ORF').replace('Orf', 'ORF')

            case _:
                return None, None, None

        return mutType, orf, snp

    @staticmethod
    def getOrfStartPos(annotation, orf):
        annotationFile = open(annotation).readlines()
        for line in annotationFile:
            if '#' not in line and 'gene=' + orf in line and line.split()[2] == 'gene':
                return int(line.split()[3])

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

    def _trimSearchResult(self, result, wtSeq, mutationType, maxPos, snp, startPos=-1, wtCodon=''):
        result = pd.read_csv(result, sep='\t', header=None)
        match mutationType:
            case 'aa':
                aa2 = snp[-1]
                print('AA')
                result.columns = ['subGenome', 'SNPbyAA', 'match', 'STsequence']
                result = result[result.STsequence.apply(lambda x: len(x) == len(wtSeq))]
                if result.empty:
                    return result

                result['STcodon'] = result.STsequence.apply(lambda x: x[maxPos - 1:maxPos + 2])
                result = result[result.STcodon.apply(lambda x: self.checkCodon(x))]
                if result.empty:
                    return result
                result = result[result.STcodon.apply(lambda x: Seq(x).translate() == aa2)]

                if result.empty:
                    return result

                result['WTcodon'] = wtCodon
                result['diffNT'] = result.apply(lambda x: [i for i in range(len(x[4])) if x[4][i] != x[5][i]], axis=1)
                result['diffNT'] = result.diffNT.apply(lambda x: x[0] if len(x) == 1 else -1)
                result['locSNP'] = result['diffNT'].apply(lambda x: x + maxPos - 1)
                result['SNPbyNT'] = result.apply(lambda x: '{}{}{}'.format(x[5][x[6]], startPos + x[6], x[4][x[6]]), axis=1)
                result['WTsequence'] = wtSeq
                return result

            case 'nt':
                print('NT')
                nt2 = snp[-1]
                result.columns = ['subGenome', 'SNPbyNT', 'match', 'STsequence']
                result['WTsequence'] = wtSeq
                result['locSNP'] = maxPos - 1
                result = result[result.STsequence.apply(lambda x: len(x) == len(wtSeq))]
                if result.empty:
                    return result

                result = result[result.STsequence.apply(lambda x: x[maxPos - 1] == nt2)]
                return result

            case _:
                return pd.DataFrame()

    def findSNPsAA(self, orf, mutation, snp, refSeq, maxPos, minPos):
        # DO NOT EXIST REFERENCE GENOME ANNOTATION
        if self.refGenomeAnnot == '':
            self.logUpdate('[warn]For Amino Acid based SNPs reference annotation needed.')
            return -1, None

        # NOT A VALID SNP
        orfStartPos = self.getOrfStartPos(self.refGenomeAnnot, orf)
        if orfStartPos == -1:
            self.logUpdate('[warn]Failure to find snp {} in reference annotaion.'.format(snp))
            return -1, None

        aa1, aa2, mutPos = mutation[0], mutation[-1], int(mutation[1:-1])
        codonStartPos = orfStartPos + (mutPos - 1) * 3 - 1
        codonEndPos = orfStartPos + mutPos * 3 - 1
        refCodon = refSeq[codonStartPos: codonEndPos]
        if aa1 != Seq(refCodon).translate():
            self.logUpdate('[warm]Failure to find SNP {} in reference genome'.format(snp))
            return -1, None

        # SEARCH
        seqWithSNP = refSeq[codonStartPos - (maxPos - 1): codonEndPos + (self.pLen1 - minPos)]
        searchResult = self._searchSnpFromStrGenome(mutation, seqWithSNP)
        trimmedResult = self._trimSearchResult(searchResult, seqWithSNP, 'aa', maxPos, snp, codonStartPos, refCodon)
        if trimmedResult.empty:
            self.logUpdate(f'[WARN] Problems occured searching snp {snp} in the strain genome or reference genome')
            return -1, None

        # RETURN RESULT
        wtSequence, stSequence, ntSNP, locSnp, found = self._parseSearchResult(searchResult=trimmedResult)
        self.logUpdate('[INFO]aa:{}:{} converted to nt:{}'.format(orf, mutation, ntSNP))
        aaOrf = '{}:{}'.format(orf, mutation)
        return found, ParaSeqs(ntSNP, aaOrf, wtSequence, stSequence, mutLoc=locSnp, probLen=self.pLen1)

    def findSNPsNT(self, mutation, snp, refSeq, maxPos, minPos):
        nt1, nt2, snpPos = mutation[0], mutation[-1], int(mutation[1:-1])
        refNT = refSeq[snpPos]
        # NOT A VALID SNP
        if nt1 != refNT:
            self.logUpdate('[warn]Failure to find SNP {} in reference genome'.format(snp))
            return -1, None

        #  SEARCH
        seqWithSNP = refSeq[snpPos - (maxPos - 1):snpPos + 1 + (self.pLen1 - minPos)]
        searchResult = self._searchSnpFromStrGenome(mutation, seqWithSNP)
        trimmedResult = self._trimSearchResult(searchResult, seqWithSNP, 'nt', maxPos, snp)
        if trimmedResult.empty:
            self.logUpdate(f'[WARN] Problems occured searching snp {snp} in the strain genome or reference genome')
            return -1, None

        # RETURN RESULT
        wtSequence, stSequence, ntSNP, locSnp, found = self._parseSearchResult(searchResult=trimmedResult)
        return found, ParaSeqs(ntSNP, '', wtSequence, stSequence, mutLoc=locSnp, probLen=self.pLen1)

    def findSNPs(self):
        minPos = min(self.posList)
        maxPos = max(self.posList)
        refSeq = self.getReferenceSeq(self.refGenome)
        if not refSeq:
            self.logUpdate('[warn]Failure to get reference sequence from reference genome.')
            self.printUsage()
        for snp in self.snpList:
            self.logUpdate('[INFO]SNP {}'.format(snp))
            mutType, orf, mutation = self.parseMutation(snp)
            found, mutSeqs = self.findSNPsAA(orf, mutation, snp, refSeq, maxPos, minPos) if mutType == 'aa' else self.findSNPsNT(mutation, snp, refSeq, maxPos, minPos)
            if found <= 0 or not found:
                self.logUpdate('[WARN] Failure to find SNP {} in strain genome'.format(snp))
                continue
            self.probesByPos[-1].append(mutSeqs)
            for pos in self.posList:
                wtProbe, stProbe = mutSeqs.getProbesWithPos(pos)
                paraSeq = ParaSeqs(mutSeqs.ntSnp, mutSeqs.aaSnp, wtProbe, stProbe, found=found, probLen=self.pLen1)
                self.probesByPos[pos].append(paraSeq)

    @staticmethod
    def makeMaskedPosBED(inputDF, outputBed):
        inputDF['seqID'] = inputDF['seqID'].apply(lambda x: x.split(';')[0])
        inputDF.to_csv(outputBed, sep='\t', header=False, index=False)

    def makeProbe1(self):
        probeLines = []
        rcpPobeLines = []
        kmerIdx = 0
        for pos in self.posList + [-1]:
            probeCSV = self.posProbeCSVs[pos]
            csvWriter = open(probeCSV, 'w')
            csvWriter.write('WT sequence,ST sequence,found,ntSNP,aaSNP\n')
            probeIdx = 0

            for p in self.probesByPos[pos]:
                csvWriter.write(f'{p.wtSeq},{p.stSeq},{p.found},{p.ntSnp},{p.aaSnp}\n')
                if pos != -1:
                    seq = p.stSeq
                    rcSeq = Seq(seq).reverse_complement()
                    probeLines.append(f'>probe{probeIdx}_{kmerIdx};{p.ntSnp}{"=" + p.aaSnp if p.aaSnp else ""};pos={pos}\n{seq}\n')
                    rcpPobeLines.append(f'>probe{probeIdx}_{kmerIdx};{p.ntSnp}{"=" + p.aaSnp if p.aaSnp else ""};pos={pos}\n{rcSeq}\n')
                probeIdx += 1

            kmerIdx += 1
            csvWriter.close()

        if not probeLines:
            self.logUpdate('[ERROR] Cannot find any SNP in strain genomes')
            self.printUsage()

        with open(self.probe1, 'w') as w:
            w.writelines(sorted(probeLines))

        with open(self.rcProbe1, 'w') as w:
            w.writelines(sorted(rcpPobeLines))

    def makeWindows(self):
        # MAKE WINDOW1
        window1 = [f'>probe{i}_{p.ntSnp}{"=" + p.aaSnp if p.aaSnp else ""}\n{p.stSeq}\n' for i, p in enumerate(self.probesByPos[-1])]
        with open(self.window1FASTA, 'w') as w:
            w.writelines(window1)

        # MAKE WINDOW 2
        msg = ProbeitUtils.getPatternPosition(self.window1FASTA, self.strGenome, self.window1TSV)
        self.logUpdate(msg, False)
        self.makeMaskedPosBED(pd.read_csv(self.window1TSV, sep='\t')[['seqID', 'start', 'end']], self.window1PosBED)
        ProbeitUtils.getWindowFasta(self.strGenome, self.window1PosBED, self.maskedGenomeFASTA, self.window2PosBED, self.window2FASTA, self.windowSize)
        ProbeitUtils.makeLookup(self.window2FASTA, self.lookup, self.window2PosBED)

    def makeProbe2(self):
        # COMPUTE MAPPABILITY
        uniqComMap = ProbeitUtils.defineFile(self.cmDir2, 'uniq.genmap.csv')
        message = '[INFO] compute mappability\n'
        match self.cmError2:
            case 0:
                ProbeitUtils.simpleComputeMappability(self.window2FASTA, self.lookup, self.window2PosBED, uniqComMap, self.pLen2, improperKmers=self.impKmers2)
            case _:
                message += ProbeitUtils.computeMappability(self.window2FASTA, self.idxDir2, self.cmError2, self.pLen2, self.cmDir2, uniqComMap, threads=self.threads, improperKmers=self.impKmers2)

        scBed = ProbeitUtils.defineFile(self.scDir2, 'result.bed')
        msg = ProbeitUtils.makeProbe(self.tempProbe2, scBed, self.window2FASTA, self.lookup, self.pLen2, self.scCoverage2, self.scEarlyStop2, self.scScore2, self.scRepeats2, uniqComMap, overlap=True)
        self.logUpdate(msg, False)
        maskDF = pd.read_csv(self.window1TSV, sep='\t')
        kmers = list(maskDF['patternName'])
        probe2Lines = []
        rcprobe2Lines = []

        for h, seq in SimpleFastaParser(open(self.tempProbe2)):
            rcSeq = Seq(seq).reverse_complement()
            kmerIndex = ProbeitUtils.parseGenmapPattern(h)[0]
            coveredSNPs = list((set([kmers[i] for i in kmerIndex])))
            probe2Lines.append(f'{":".join(coveredSNPs)}\n{seq}\n')
            rcprobe2Lines.append(f'{":".join(coveredSNPs)}\n{rcSeq}\n')

        with open(self.probe2, 'w') as w:
            for i, line in enumerate(sorted(probe2Lines)):
                w.write(f'>cap{i}\t{line}')

        with open(self.rcprobe2, 'w') as w:
            for i, line in enumerate(sorted(rcprobe2Lines)):
                w.write(f'>cap{i}\t{line}')

    def run(self):
        # START MAKING 1st PROBEs
        self.logUpdate("[INFO]make 1st probes")

        # FIND SNPs FROM STRAIN GENOME
        self.findSNPs()

        # MAKE PROBE1
        self.makeProbe1()

        # COMPLETE MAKING PROBE1
        if not self.needProbe2:
            print('COMPLETE!!!')
            return

        # START MAKING 2ND PROBEs
        self.logUpdate("[INFO]make 2nd probes")

        # MAKE WINDOW FOR PROBE2
        self.makeWindows()

        # THERMO FILTER
        if self.doThermoFilter2:
            ProbeitUtils.extractKmers(self.window2FASTA, self.kmers2FASTA, self.pLen2)
            thermoFilter2 = ThermoFilter(self.kmers2FASTA, self.pLen2, self.thermoCalcTSV2)
            msg, self.impKmers2 = thermoFilter2.run()
            self.logUpdate(msg)

        # MAKE PROBE2
        self.makeProbe2()

        # COMPLETE MAKING PROBE2
        print('COMPLETE!!!')
        return

    @staticmethod
    def printUsage(wrong=''):
        if wrong:
            print(f"[ERR] Unknown parameter found: {wrong}")
        print("Probeit snp")
        print("It generates a probe set which detect input amino acid SNPs from strain genome.")
        print("probeit snp -r [REF GENOME] -s [STR GENOME] -p [positions] -m [SNPs] -o [DIR] -a [REF ANNOTATION]")

        print("Usage")
        print(" -h|--help NONE")
        print('\t Show usage')

        print("REQUIRED OPTIONS")
        print(" -r|--reference FASTA file")
        print("\t The wild type genome.")
        print(" -s|--strain FASTA file")
        print("\t The strain Genome.")
        print(" -p|--positions COMMA SEPARATED INT ARRAY")
        print("\t Position List: The position of this indicates the position of the SNP on the 1st Probes")
        print(" -m|--mutations COMMA SEPARATED SNP ARRAY")
        print("\t SNP List. Both amino acid differences and nucleotide differences are allowed.")
        print(" -o|--output DIR")
        print("\t Output directory. The Directory is automatically created by Probeit.")
        print(" -a|--annotation GFF file")
        print("\t The wild type genome annotation. Only required when using amino acid differences in the -m option.")

        print("ADDITIONAL OPTIONS")
        print(" --not-make-probe2 NONE")
        print("\t Use it when you DO NOT need to make 2nd probes")
        print(" --threads INT[8]")
        print("\t number of CPU-cores used")
        print(" --max-window NONE")
        print("\t When you need maximum window mode, then use this option. Default window mode is minimum window.")
        print(" --window-size INT[200]")
        print("\t size of windows for 2nd probes")
        print(" --probe1-len INT[40]")
        print("\t Length of 1st Probes")
        print(" --probe2-len INT[20]")
        print("\t Length of 2nd Probes")
        print(" --not-make-probe2 NONE")
        print("\t Use it when you DO NOT need to make 2nd probes")

        print("OPTIONS FOR GENMAP PART: Genmap calculates mappability by summarizing k-mers.")
        print(" --probe2-error INT[1]")
        print("\t The number of error allowed in 2nd Probes")

        print("OPTIONS FOR SETCOVER PART: Setcover makes probesets cover targets with the minimum number of probes.")
        print(" --probe2-cover INT[1]")
        print("\t The number of times each 1st Probe should be covered by 2nd Probes")
        print(" --probe2-repeat INT[10]")
        print("\t The number of random iterations when minimizing 2nd Probes")
        print(" --probe2-earlystop INT[99]")
        print("\t Early stop picking new probes if X% of sequences are covered at least N(--probe2-cover) times.")
        quit()