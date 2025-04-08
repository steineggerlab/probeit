#!/usr/bin/env python
from Bio.SeqIO.FastaIO import SimpleFastaParser
from Bio.Seq import Seq, reverse_complement
from pandas import read_csv, merge, DataFrame
from primer3 import calc_tm, calc_hairpin_tm, calc_homodimer_tm
from numpy import vectorize
from os import path, system, makedirs
from shutil import rmtree
from subprocess import Popen, PIPE
from re import compile, findall
from .config import Config

class Kmer:
    def __init__(self, line="", idx=-1, sPos=-1):
        if line:
            self.idx, self.sPos = map(int, line.split(','))
            return
        self.idx = idx
        self.sPos = sPos

    def getIdx(self):
        return self.idx

    def getStartPos(self):
        return self.sPos

    def __eq__(self, other):
        if self.__class__ != other.__class__:
            return False

        if self.idx != other.idx:
            return False

        if self.sPos != other.sPos:
            return False

        return True

    def __ne__(self, other):
        if self.__class__ != other.__class__:
            return True

        if self.idx != other.idx:
            return True

        if self.sPos != other.sPos:
            return True

        return False

    def __lt__(self, other):
        if self.idx >= other.idx:
            return False
        if self.sPos >= other:
            return False
        return True

    def __le__(self, other):
        if self.idx > other.idx:
            return False
        if self.sPos > other:
            return False
        return True

    def __gt__(self, other):
        if self.idx <= other.idx:
            return False
        if self.sPos <= other:
            return False
        return True

    def __ge__(self, other):
        if self.idx < other.idx:
            return False
        if self.sPos < other:
            return False
        return True

    def __hash__(self):
        return hash((self.idx, self.sPos))

    def __str__(self):
        return f'{self.idx},{self.sPos}'

    def getStr(self):
        return f'{self.idx},{self.sPos}'

class BedLine:
    def __init__(self, idx, sPos, ePos):
        self.idx = idx
        self.sPos = sPos
        self.ePos = ePos

    def getLine(self):
        return f'{self.idx}\t{self.sPos}\t{self.ePos}\n'

# To implement increment operator like i++ in C.
class Incrementer:
    def __init__(self, start=-1):
        self.value = start
        
    def get(self):
        self.value += 1
        return self.value

    def reset(self, val=-1):
        self.value = val

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
        if self.__class__ != other.__class__:
            return False
        if self.ntSnp != other.ntSnp:
            return False
        if self.wtSeq != other.wtSeq:
            return False
        if self.stSeq != other.stSeq:
            return False
        return True

    def getProbesWithPos(self, pos):
        start = self.mutLoc - pos + 1
        end = self.mutLoc - pos + self.probLen + 1
        return self.wtSeq[start:end], self.stSeq[start:end]

class ProbeitUtils:
    # RUN COMMANDLINE
    @staticmethod
    def runCommand(command, verbose=False):
        print('[CLI] ' + command)
        if verbose:
            commandList = command.split()
            sp = Popen(commandList, stdout=PIPE, stderr=PIPE)
            stdout, stderr = sp.communicate()
            return stdout.decode('UTF-8'), stderr.decode('UTF-8')

        system(command)

    # ARGUMENTS PARSING
    @classmethod
    def getUserArgs(cls, args):
        return ' '.join(['{} {}'.format(i[0], i[1]) if len(i) == 2 else i for i in args])

    # FILES AND DIRECTORIES RELATED
    @staticmethod
    def defineFile(directory, fileName):
        return directory + fileName

    @classmethod
    def defineDirectory(cls, dirName, root='', make=True):
        directory = f'{root}{dirName}{path.sep}'
        if make:
            makedirs(directory)
        return directory

    @staticmethod
    def getFileName(file, keepFilenameExtension=True):
        file = file.split(path.sep)[-1]
        if keepFilenameExtension:
            return file

        return '.'.join(file.split('.')[:-1])

    @classmethod
    def delDir(cls, directory):
        if path.isdir(directory):
            print(f'[INFO] The directory named {directory} is removed not.')
            rmtree(directory)
            return

        print(f'[INFO] The directory named {directory} does not exist.')

    # FASTA FILES RELATED
    @classmethod
    def getSubseqFasta(cls, coordinateBed, inputFasta, outputFasta):
        command1 = "seqkit subseq --quiet --bed {} {} > {}".format(coordinateBed, inputFasta, outputFasta)
        cls.runCommand(command1)
        return '[CLI] ' + command1 + '\n'

    @classmethod
    def getWindowFasta(cls, genomeFasta, maskingBed, maskedGenomeFasta, windowBed, windowFasta, windowSize):
        command1 = "bedtools maskfasta -fi {} -bed {} -fo {}".format(genomeFasta, maskingBed, maskedGenomeFasta)
        cls.runCommand(command1)
        inputDF = read_csv(maskingBed, sep='\t', header=None)
        inputDF[1] = inputDF[1].apply(lambda x: x - windowSize)
        inputDF[2] = inputDF[2].apply(lambda x: x + windowSize)
        inputDF.to_csv(windowBed, sep='\t', header=False, index=False)

        command2 = "bedtools getfasta -fi {} -bed {} > {}".format(maskedGenomeFasta, windowBed, windowFasta)
        cls.runCommand(command2)

    @staticmethod
    def simplifyFastaHeaders(inputFasta, outputFasta):
        def _getSimpleHeader(header):
            return header.split()[0]

        with open(outputFasta, 'w') as w:
            w.writelines([f'>{_getSimpleHeader(h)}\n{seq}\n' for h, seq in SimpleFastaParser(open(inputFasta))])


    @staticmethod
    def extractKmers(genomeFasta, kmersFasta, pLen):
        def writeKmer(seqAndKmers):
            seq = seqAndKmers[0]
            kmers = seqAndKmers[1]
            repKmer = kmers[0]
            strKmers = [kmer.getStr() for kmer in kmers]
            return f'>{repKmer.getStr()};{"|".join(strKmers)}\n{seq}\n'

        kmers = {}
        genomeIdx = Incrementer()
        for _, s in SimpleFastaParser(open(genomeFasta)):
            l = len(s)
            for pos in range(l - pLen + 1):
                kmer = Kmer(idx=genomeIdx.get(), sPos=pos)
                seq = s[pos:pos + pLen].upper()
                kmers.setdefault(seq, [])
                kmers[seq].append(kmer)

        with open(kmersFasta, 'w') as w:
            w.writelines(map(writeKmer, kmers.items()))

    @classmethod
    def getPatternPosition(cls, patternFasta, genomeFasta, positionsTSV):
        command = "seqkit locate -f {} {} > {}".format(patternFasta, genomeFasta, positionsTSV)
        cls.runCommand(command)
        df = read_csv(positionsTSV, sep='\t')
        df.sort_values(by=list(df.columns), ignore_index=True).to_csv(positionsTSV, sep='\t', index=False)
        return '[CLI] {}\n'.format(command)

    @classmethod
    def clusterGenome(cls, inputFile, outputFile, directory, seqIdentity, threads):
        cmd = 'mmseqs easy-linclust {} {} {} -v 3 --kmer-per-seq-scale 0.5 --kmer-per-seq 1000 --min-seq-id {} --cov-mode 1 -c 0.95 --remove-tmp-files 0 --threads {}'
        cmd = cmd.format(inputFile, outputFile, directory, seqIdentity, threads)
        stdout, stderr = cls.runCommand(cmd, verbose=True)
        msg = stdout + stderr
        outputFile += '_rep_seq.fasta'
        return msg, outputFile

    @classmethod
    def ridNegKmers(cls, posKmers, negative, output, outputDir, seqInProbe, thread):
        tempDir = cls.defineDirectory("tmp", root=outputDir)
        command = f"mmseqs easy-search {posKmers} {negative} {output} {tempDir} -v 3 --spaced-kmer-mode 0 -k 13 --mask 0 --format-output query,target,pident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue -c 0.9 --min-seq-id {seqInProbe} --cov-mode 2 --alignment-mode 4 --search-type 3 --threads {thread} --remove-tmp-files 0"
        stdOut, stdErr = cls.runCommand(command, verbose=True)
        print(stdOut)
        print(stdErr)
        return f'[CLI] {command}\n' + stdOut + stdErr

    # COMPUTE MAPPABILITY
    @classmethod
    def simpleComputeMappability(cls, genome, lookup, positionBED, outputCSV, pLen, improperKmers):
        positiveKmers = set()
        w = open(outputCSV, 'w')
        for line in open(positionBED):
            gIdx, sPos, ePos = map(int, line.split())

            for i in range(sPos, ePos):
                positiveKmers.add(Kmer(idx=gIdx, sPos=i))

        seqAndKmers = {}
        idx = Incrementer()
        for h, seq in SimpleFastaParser(open(genome)):
            l = len(seq)
            for i in range(l - pLen + 1):
                kmer = Kmer(idx=idx.get(), sPos=i)
                seq = seq[i:i + pLen].upper()
                if set(seq) - Config.nucleotideSet:
                    continue
                seqAndKmers.setdefault(seq, [])
                seqAndKmers[seq].append(kmer)

        for seq, kmers in seqAndKmers.items():
            repKmer = kmers[0]
            isPositive = repKmer in positiveKmers
            isImproper = seq in improperKmers
            if isPositive and not isImproper:
                kmerStr = '|'.join([kmer.getStr() for kmer in kmers])
                w.write(f'{repKmer};{kmerStr}\n')
        w.close()

    @classmethod
    def computeMappability(cls, genome, indexDir, error, length, outputDir, outputCSV, threads, improperKmers):
        w = open(outputCSV, 'w')
        cls.delDir(indexDir)
        command1 = f"genmap index -F {genome} -I {indexDir}  > /dev/null"
        command2 = f"genmap map --no-reverse-complement -E {error} --csv -K {length} -t -b --frequency-large -I {indexDir} -O {outputDir} -T {threads}  > /dev/null"
        cls.runCommand(command1)
        cls.runCommand(command2)
        inputCSV = cls.defineFile(outputDir, f'{cls.getFileName(genome, False)}.genmap.csv')
        for line in open(inputCSV):
            parsedKmers = findall(r'[0-9]+,[0-9]+',line)
            if not parsedKmers:
                continue
            repKmer = Kmer(line=parsedKmers[0])
            firstKmer = Kmer(line=parsedKmers[1])
            if repKmer == firstKmer and repKmer not in improperKmers:
                w.write(line)
        w.close()
        return '[CLI] {}\n[CLI] {}\n'.format(command1, command2)

    # BED FILES RELATED
    @classmethod
    def getSubtractedBed(cls, positiveBed, negativeBed, bedFileName):
        command = "bedtools subtract -a {} -b {} > {}".format(positiveBed, negativeBed, bedFileName)
        cls.runCommand(command)
        return '[CLI] ' + command

    # LOOKUP FILES RELATED
    @classmethod
    def makeLookup(cls, windowFasta, lookup, genomePos=''):
        i = Incrementer()
        with open(lookup, 'w') as w:
            w.writelines([f'{i.get()}\t{header.strip()}\n' for header, _ in SimpleFastaParser(open(windowFasta))])

        if not genomePos:
            return

        i.reset()
        with open(genomePos, 'w') as w:
            w.writelines([BedLine(i.get(), 0, len(seq)).getLine() for _, seq in SimpleFastaParser(open(windowFasta))])

    # SNPs RELATED
    @classmethod
    def searchSNPs(cls, workDir, inputFasta, strGenomeFasta, result, kmer, threads):
        searchDir = cls.defineDirectory('search', root=workDir, make=True)
        tempDir = cls.defineDirectory('temp', root=searchDir, make=False)
        searchDB = cls.defineFile(searchDir, 'searchDB')
        strainDB = cls.defineFile(searchDir, 'strainDB')
        aln = cls.defineFile(searchDir, 'mmseqs.aln')
        cmd1 = f'mmseqs createdb {inputFasta} {searchDB}'
        cmd2 = f'mmseqs createdb {strGenomeFasta} {strainDB}'
        cmd3 = f'mmseqs search {searchDB} {strainDB} {aln} {tempDir} --search-type 3 -k {kmer} --threads {threads}'
        cmd4 = f'mmseqs convertalis {searchDB} {strainDB} {aln} {result} --format-output target,query,tseq,tstart,tend --search-type 3 --threads {threads}'
        out1, err1 = cls.runCommand(cmd1, verbose=True)
        out2, err2 = cls.runCommand(cmd2, verbose=True)
        out3, err3 = cls.runCommand(cmd3, verbose=True)
        out4, err4 = cls.runCommand(cmd4, verbose=True)
        df = read_csv(result, sep='\t', header=None)
        df.columns = ['substr', 'snp', 'strseq', 'start', 'end']
        df['aln'] = df.apply(lambda x: x[2][int(x[3] - 1):int(x[4])], axis=1)
        df['len'] = df.aln.apply(lambda x: len(x) - 1)
        df = df[['substr', 'snp', 'len', 'aln']]
        df.to_csv(result, header=False, index=False, sep='\t')
        print(err1 + err2 + err3 + err4)
        ProbeitUtils.delDir(searchDir)
        return out1 + out2 + out3 + out4

    # PARSE GENMAP KMER PATTERNS
    @staticmethod
    def parseKmers(line):
        return [Kmer(kmerHeader) for kmerHeader in findall(r'[0-9]+,[0-9]+', line)]

    # @staticmethod
    # def parseKmers2(line):
    #     return  [tuple(map(int, kmer.split(","))) for kmer in findall(r'[0-9]+,[0-9]+', line)]

    @staticmethod
    def parseGenmapPattern(header):
        p1 = compile('[0-9]+,')
        p2 = compile(',[0-9]+')
        return [int(i[:-1]) for i in p1.findall(header)], [int(i[1:]) for i in p2.findall(header)]

    # ARGUMENTS PARSING
    @staticmethod
    def getArgList(value, isInt=False):
        if isInt:
            return [int(i) for i in value.split(',')]

        return value.split(',')

    # OTHERS
    @classmethod
    def setCover(cls, coverage, length, eStop, dist, reps, mapCSV, genome, lookup, setcoverResultBed, probeLen):
        def writeSetcoverResult(line):
            repKmer = ProbeitUtils.parseKmers(line)[0]
            matchedKmers = line.split(';')[1].strip()
            idx = repKmer.idx
            pos = repKmer.sPos
            return f"{genomeAndIdx[idx]}\t{pos}\t{pos + probeLen}\t{matchedKmers}\n"

        filePath = path.sep.join(path.realpath(__file__).split(path.sep)[:-1])
        setcoverPath = '{}{}{}{}{}'.format(filePath, path.sep, 'setcover', path.sep, 'setcover')
        setcoverPath = setcoverPath if path.exists(setcoverPath) else 'setcover'
        command = " -c {} -l {} -p {} -d {} -i {} {} {}".format(coverage, length, eStop, dist, reps, mapCSV, genome)
        stdOut, stdErr = cls.runCommand(setcoverPath + command, verbose=True)
        genomeAndIdx = dict()
        print(stdOut)
        for line in open(lookup):
            idx = int(line.split()[0].strip())
            genome = line.split()[1].strip()
            genomeAndIdx[idx] = genome

        with open(setcoverResultBed, 'w') as w:
            w.writelines(map(writeSetcoverResult, stdOut.strip().split('\n')))

        return f'{setcoverPath}{command}\n{stdOut}{stdErr}'

    @classmethod
    def makeProbe(cls, output, scResultBed, window, lookup, probeLen, scCoverage, scEarlyStop, scScore, scRepeats, uniqComMap, overlap=False, ):
        scLen = 1 if overlap else probeLen
        # SETCOVER
        message = "[INFO] minimize probe set\n"
        message += cls.setCover(scCoverage, scLen, scEarlyStop, scScore, scRepeats, uniqComMap, window, lookup, scResultBed, probeLen)
        # MAKE PROBEs
        message += cls.getSubseqFasta(scResultBed, window, output)
        return message


class ThermoFilter:
    def __init__(self, kmersFASTA, pLen, outputTSV, isLigational=False):
        self.kmersFASTA = kmersFASTA
        self.pLen = pLen
        self.outputTSV = outputTSV
        self.isLigational = isLigational
        self.minGC = Config.getMinGC()
        self.maxGC = Config.getMaxGC()
        self.maxHomoDimerTm = Config.getMaxhomoDimerTm()
        self.maxHairpinTm = Config.getMaxhairpinTm()
        self.minProbeTm = Config.getMinProbeTm()
        maxRepeat = Config.getMaxRepeat()
        self.problematicSeqs = ['A' * maxRepeat, 'T' * maxRepeat, 'C' * maxRepeat, 'G' * maxRepeat]
        self.kmerList = []
        self.query = ['p1_intrinsic_probs', f'p1_tm<{self.minProbeTm}', f'p1_GC_perc<{self.minGC}',
                      f'p1_GC_perc>{self.maxGC}', f'p1_homodimer_tm>{self.maxHomoDimerTm}',
                      f'p1_hairpin_tm>{self.maxHairpinTm}']
        self.sorter = ['p1_hairpin_tm', 'p1_homodimer_tm']
        self.ascending = [True, True, True, True] if self.isLigational else [True, True]
        if self.isLigational:
            self.query += ['p2_intrinsic_probs', f'p2_tm<{self.minProbeTm}', f'p2_GC_perc<{self.minGC}',
                           f'p2_GC_perc>{self.maxGC}', f'p2_homodimer_tm>{self.maxHomoDimerTm}',
                           f'p2_hairpin_tm>{self.maxHairpinTm}']
            self.sorter += ['p2_hairpin_tm', 'p2_homodimer_tm']

    def makeLogMessage(self, numKmers, numImpKmers):
        msg = '[INFO]filter probes with thermodynamic features\n'
        msg += f"\tMinimum Tm: {self.minProbeTm}\n"
        msg += f"\tMinimum GC percentage: {self.minGC}\n"
        msg += f"\tMaximum GC percentage: {self.maxGC}\n"
        msg += f"\tHomodimer maximum Tm: {self.maxHomoDimerTm}\n"
        msg += f"\tHairpin maximum Tm: {self.maxHairpinTm}\n"
        msg += f"[INFO] {numKmers} ligation probe candidate sets inputted.\n"
        msg += f'[INFO] Improper Ligation probe candidates are removed: {numImpKmers}\n'
        return msg

    @staticmethod
    def getContentGC(oligo):
        gcCount = oligo.count('G') + oligo.count('C')
        return gcCount / len(oligo)

    def hasLowComplexity(self, candidate_oligo):
        for s in self.problematicSeqs:
            if s in candidate_oligo:
                return True

        return False

    @staticmethod
    def getImpKmersSet(impKmers):
        return set(sum(impKmers['id'].apply(lambda x: ProbeitUtils.parseKmers(x)[1:]), []))

    def getKmersDF(self):
        with open(self.kmersFASTA) as f:
            identity = []
            posStartList = []
            posEndList = []
            kmers = []
            rc = []
            p1 = []
            p2 = []
            for title, kmer in SimpleFastaParser(f):
                split_name = title.split('\t', 2)
                identity.append(split_name[0])
                pStart = int(split_name[1]) if len(split_name) > 1 else 0
                pEnd = pStart + self.pLen
                posStartList.append(pStart)
                posEndList.append(pEnd)
                kmers.append(kmer)
                rcKmer = str(reverse_complement(Seq(kmer)))
                rc.append(rcKmer)
                if not self.isLigational:
                    p1.append(rcKmer)
                    continue
                mid_pos = round(self.pLen / 2)
                p1.append(rcKmer[0: mid_pos])
                p2.append(rcKmer[mid_pos: self.pLen])

        self.kmerList = p1 + p2
        cols = ['id', 'chromStart', 'chromEnd', 'genome_segment', 'rc', 'p1']
        outputDF = DataFrame(list(zip(identity, posStartList, posEndList, kmers, rc, p1)), columns=cols)
        if self.isLigational:
            outputDF['p2'] = p2

        return outputDF

    def getThermoFeaturesDF(self):
        thermoFeatures = DataFrame(list(set(self.kmerList)), columns=['p'])
        thermoFeatures = thermoFeatures.assign(ultimate_base=thermoFeatures['p'].str[-1])
        thermoFeatures = thermoFeatures.assign(penultimate_base=thermoFeatures['p'].str[-2])
        thermoFeatures = thermoFeatures.assign(tm=vectorize(calc_tm)(thermoFeatures['p']))
        thermoFeatures['hairpin_tm'] = list(map(calc_hairpin_tm, thermoFeatures['p']))
        thermoFeatures = thermoFeatures.assign(homodimer_tm=vectorize(calc_homodimer_tm)(thermoFeatures['p']))
        thermoFeatures = thermoFeatures.assign(intrinsic_probs=vectorize(self.hasLowComplexity)(thermoFeatures['p']))
        thermoFeatures = thermoFeatures.assign(GC_perc=vectorize(self.getContentGC)(thermoFeatures['p']))
        return thermoFeatures

    def getThermoFeaturedKmersDF(self, kmersDF, thermoFeaturesDF):
        joinedDF = merge(kmersDF, thermoFeaturesDF.add_prefix('p1_'), how='left', left_on='p1', right_on='p1_p')
        if self.isLigational:
            joinedDF = merge(joinedDF, thermoFeaturesDF.add_prefix('p2_'), how='left', left_on='p2', right_on='p2_p')
        return joinedDF

    @staticmethod
    def dfFilter(df, query, sorter, ascending):
        filtered = df.query(' | '.join(query))
        return filtered.sort_values(by=sorter, ascending=ascending)

    def run(self):
        kmersDF = self.getKmersDF()
        thermoFeaturesDF = self.getThermoFeaturesDF()
        thermoFeaturedKmersDF = self.getThermoFeaturedKmersDF(kmersDF, thermoFeaturesDF)
        thermoFeaturedKmersDF.to_csv(self.outputTSV, sep='\t')
        impKmers = self.dfFilter(thermoFeaturedKmersDF, self.query, self.sorter, self.ascending)
        return self.makeLogMessage(len(kmersDF), len(impKmers)), self.getImpKmersSet(impKmers)

