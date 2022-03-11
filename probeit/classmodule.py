#!/usr/bin/env python
from operator import ne
from .config import Config
from sys import stderr
from defusedxml import DTDForbidden
from pandas.core.accessor import DirNamesMixin
from Bio.SeqIO.FastaIO import SimpleFastaParser
from Bio.Seq import Seq
from Bio.Seq import reverse_complement
import pandas as pd
import primer3
import numpy as np
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
        start = self.mutLoc - pos + 1
        end = self.mutLoc - pos + self.probLen + 1
        return self.wtSeq[start:end], self.stSeq[start:end] 


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
    def delDir(cls, directory):
        if os.path.isdir(directory):
            print(f'[INFO] The directory named {directory} is removed not.')
            shutil.rmtree(directory)
        else:
            print(f'[INFO] The directory named {directory} does not exist.')

    @classmethod
    def sortFasta(cls, inputFasta):
        fastaList = []
        fastaList = [(h,s) for h, s in SimpleFastaParser(open(inputFasta))]
        fastaList = sorted(['>{}\n{}\n'.format(i[0], i[1]) for i in fastaList])
        with open(inputFasta, 'w') as w:
            w.writelines(fastaList)

    # TO CALL SEQKIT MODULES
    @classmethod
    def renameFasta(cls, oldFasta, newFasta):
        command = "seqkit rename {} > {} -w 0".format(oldFasta, newFasta)
        cls.runCommand(command)

    @classmethod
    def getSubseqFasta(cls, coordinateBed, inputFasta, outputFasta):
        command1 = "seqkit subseq --quiet --bed {} {} > {}".format(coordinateBed, inputFasta, outputFasta)
        cls.runCommand(command1)
        return '[CLI] '+command1+'\n'

    @classmethod
    def getPatternPosition(cls, patternFasta, genomeFasta, positonsTSV):
        command = "seqkit locate -f {} {} > {}".format(patternFasta, genomeFasta, positonsTSV)
        cls.runCommand(command)
        df = pd.read_csv(positonsTSV, sep='\t')
        df.sort_values(by=list(df.columns), ignore_index=True).to_csv(positonsTSV, sep='\t', index=False)
        return '[CLI] {}\n'.format(command)

    @staticmethod
    def getFileName(file, keepFilenameExtension=True):
        file = file.split(os.path.sep)[-1]
        if keepFilenameExtension:
            return file 
        else:
            return '.'.join(file.split('.')[:-1])

    # COMPUTE MAPPABILITY
    @classmethod
    def simpleComputeMappability(cls, genome, lookup, positionBED, outputCSV, pLen):
        positive = set()
        genomeNames = [i.split()[1].strip() for i in open(lookup)]
        for line in open(positionBED):
            g, s, e = line.strip().split()
            s = int(s)
            e = int(e)
            g = genomeNames.index(g)
            for i in range(s,e):
                positive.add(f'{g},{i}')
        kmers = {}
        idx = 0
        alphabets = {'A', 'C', 'G', 'T'}
        for h,s in SimpleFastaParser(open(genome)):
            l = len(s)
            for i in range(l-pLen+1):
                name = f'{idx},{i}'
                seq = s[i:i+pLen].upper()
                if set(seq) == alphabets&set(seq):
                    kmers.setdefault(seq,[])
                    kmers[seq].append(name)
            idx += 1
        w = open(outputCSV, 'w')
        for s, names in kmers.items():
            rep = names[0]
            isPositive = {rep}&positive!=set()
            if isPositive:
                w.write(f'{rep};{"|".join(names)}\n')
        w.close()
        return

    @classmethod
    def computeMappability(cls, genome, indexDir, error, len, outputDir, outputCSV, positiveCoords='', threads=8,  thermoImpKmers=set()):
        selector = f'-S {positiveCoords}' if positiveCoords else ''
        cls.delDir(indexDir)
        command0 = " > /dev/null"
        command1 = "genmap index -F {} -I {}".format(genome, indexDir)
        cls.runCommand(command1 + command0)
        command2 = "genmap map --no-reverse-complement -E {} {} --csv -K {} -t -b --frequency-large -I {} -O {} -T {}"
        command2 = command2.format(error, selector, len, indexDir, outputDir, threads)
        cls.runCommand(command2 + command0)
        inputCSV = cls.defineFile(outputDir, f'{cls.getFileName(genome, False)}.genmap.csv')
        w = open(outputCSV, 'w')
        for line in open(inputCSV):
            [header, first] = line.split('|')[0].strip().split(';')
            if header == first and not {header}&thermoImpKmers:
                w.write(line)
        w.close()
        return '[CLI] {}{}\n[CLI] {}{}\n'.format(command1, command0, command2, command0)

    # TO CALL BEDTOOLS MODULES
    @classmethod
    def getSubtractedBed(cls, positiveBed, negativeBed, bedFileName):
        command = "bedtools subtract -a {} -b {} > {}".format(positiveBed, negativeBed, bedFileName)
        cls.runCommand(command)
        return '[CLI] '+command

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

    @classmethod
    def makeLookup(cls, windowFasta, lookup, genomePos=''):
        headers = [header.strip() for header, _ in SimpleFastaParser(open(windowFasta))] 
        w = open(lookup, 'w')
        for i, header in enumerate(headers):
            w.write(f'{i}\t{header}\n')
        if genomePos: 
            lengths = [len(seq) for _, seq in SimpleFastaParser(open(windowFasta))]
            w = open(genomePos, 'w')
            for i, length in enumerate(lengths):
                w.write(f'{i}\t0\t{length}\n')
    
        

    @classmethod
    def searchSNPs(cls, workDir, inputFasta, strGenomeFasta, result, kmer=12, threads=8):
        searchDir = cls.defineDirectory('search', root=workDir, make=True)
        tempDir = cls.defineDirectory('temp', root=searchDir, make=False)
        searchdb = cls.defineFile(searchDir, 'searchDB' )
        strdb =  cls.defineFile(searchDir, 'strainDB')
        aln = cls.defineFile(searchDir , 'mmseqs.aln')
        cmd0 = ' --threads {}'.format(threads)
        cmd1 = 'mmseqs createdb {} {}'
        cmd2 = 'mmseqs createdb {} {}'
        cmd3 = 'mmseqs search {} {} {} {} --search-type 3 -k {}'
        cmd4 = 'mmseqs convertalis {} {} {} {} --format-output target,query,tseq,tstart,tend --search-type 3'
        out1, err1 = cls.runCommand(cmd1.format(inputFasta, searchdb), verbose=True)
        out2, err2 = cls.runCommand(cmd2.format(strGenomeFasta, strdb), verbose=True)
        out3, err3 = cls.runCommand(cmd3.format(searchdb, strdb, aln, tempDir, kmer, threads) + cmd0, verbose=True)
        out4, err4 = cls.runCommand(cmd4.format(searchdb, strdb, aln, result, threads) + cmd0, verbose=True)
        df = pd.read_csv(result, sep='\t', header=None)
        df.columns = ['substr', 'snp', 'strseq', 'start', 'end']
        df['aln'] = df.apply(lambda x: x[2][int(x[3]-1):int(x[4])], axis=1)
        df['len'] = df.aln.apply(lambda x: len(x)-1)
        df = df[['substr', 'snp', 'len', 'aln']]
        df.to_csv(result, header=False, index=False, sep='\t')
        print(err1 + err2 + err3 + err4)
        ProbeitUtils.delDir(searchDir)
        return out1 + out2 + out3 + out4

    @classmethod
    def clusterGenome(cls, input, output, directory, seqIdentity, threads=8):
        command = (
            f'mmseqs easy-linclust {input} {output} {directory} -v 3 --kmer-per-seq-scale 0.5 --kmer-per-seq 1000 ' +
            f'--min-seq-id {seqIdentity} --cov-mode 1 -c 0.95 --remove-tmp-files 0 --threads {threads}' 
        )
        stdout, stderr = cls.runCommand(command, verbose=True)
        return output+'_rep_seq.fasta',  stdout + stderr

    @classmethod
    def ridNegKmers(cls, posKmers, negative, output, outputDir, seqInProbe, thread=8):
        tempDir = cls.defineDirectory("tmp", root=outputDir)
        command = (
            f'mmseqs easy-search {posKmers} {negative} {output} {tempDir} -v 3 --spaced-kmer-mode 0 -k 13 --mask 0 ' +
            f'--format-output query,target,pident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue ' +
            f'-c 0.9 --min-seq-id {seqInProbe} --cov-mode 2 --alignment-mode 4 --search-type 3 --threads {thread}' + ' --remove-tmp-files 0'
        )
        stdOut, stdErr = cls.runCommand(command, verbose=True)
        print(stdOut)
        print(stdErr)
        return f'[CLI] {command}\n' + stdOut + stdErr

    @staticmethod
    def parseKmers(line):
        return re.findall(r'[0-9]+,[0-9]+', line)

    @classmethod
    def setCover(cls, coverage, length, eStop, dist, reps, mapCSV, genome, lookup, setcoverResultBed, probeLen):
        filePath = os.path.sep.join(os.path.realpath(__file__).split(os.path.sep)[:-1])
        setcoverPath = '{}{}{}{}{}'.format(filePath, os.path.sep, 'setcover', os.path.sep, 'setcover')
        setcoverPath = setcoverPath   if os.path.exists(setcoverPath) else 'setcover'
        command = " -c {} -l {} -p {} -d {} -i {} {} {}".format(coverage, length, eStop, dist, reps, mapCSV, genome)
        stdOut, stdErr = cls.runCommand(setcoverPath + command, verbose=True)
        genomeAndIdx = dict()
        print(stdOut)
        for line in open(lookup):
            idx = int(line.split()[0].strip())
            genome = line.split()[1].strip()
            genomeAndIdx[idx] = genome
        with open(setcoverResultBed, 'w') as w:
            for line in stdOut.strip().split('\n'):
                matchedKmers = line.split(';')[1].strip()
                idx = int(line.split(';')[0].split(',')[0])
                pos = line.split(';')[0].split(',')[1]
                w.write('\t'.join([genomeAndIdx[idx], pos, str(int(pos) + probeLen), matchedKmers]) + '\n')
        return f'{setcoverPath}{command}\n{stdOut}{stdErr}'    

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
            cls, output, scResultBed, window, lookup, probeLen, 
            scCoverage, scEarlyStop, scScore, scRepeats, uniqComMap, 
            overlap=False, 
    ):
        scLen = 1 if overlap else probeLen 
        # SETCOVER
        message = "[INFO] minimize probe set\n"
        message += ProbeitUtils.setCover(scCoverage, scLen, scEarlyStop, scScore, scRepeats, uniqComMap, window, lookup, scResultBed, probeLen)
        # MAKE PROBEs
        message += ProbeitUtils.getSubseqFasta(scResultBed, window, output)
        return message

    @staticmethod
    def parseGenmapPattern(header):
        p1 = re.compile('[0-9]+,')
        p2 = re.compile(',[0-9]+')
        return [int(i[:-1]) for i in p1.findall(header)], [int(i[1:]) for i in p2.findall(header)]

    @staticmethod
    def defineFile(directory, fileName):
        return directory + fileName

    @classmethod
    def defineDirectory(cls, dirName, root='', make=True):
        directory = f'{root}{dirName}{os.path.sep}'
        if make:
            os.makedirs(directory)
        return directory


class Probeit:
    args = []
    subWork = None

    def __init__(self, args):
        self.args = args

    def do(self):
        if self.args == [] or self.args[0] == '-h' or self.args[0] == '--help':
            self.printUsage()
            return
        elif self.args[0] == 'posnegset':
            print('CURRENT: ', os.getcwd())
            subWork = PosNegSet(self.args[1:])
            subWork.run()
            return
        elif self.args[0] == 'snp':
            print('CURRENT: ', os.getcwd())
            subWork = SNP(self.args[1:])
            subWork.run()
            return
        else:
            self.printUsage()
            return

    @staticmethod
    def printUsage():
        print("PROBEIT V2.2")
        print('probeit <workflow> [<args>]')
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
        'not-make-probe2', 'remove-reduncancy', 'not-thermo-filter', 'ligation-probe', 'probe1-len=', 'probe2-len=', 'window-size=', 'threads=',
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
    ridNegId = 0.90   # for mmseqs easy-search
    cmError1 = 0
    cmError2 = 1
    scCoverage1 = 1
    scCoverage2 = 1
    doRemoveRedundancy = False
    doDesignProbe2 = True
    doThermoFilter = True # thermofilter
    isLigationProbe = False # thermofilter
    windowSize = 200
    scEarlyStop1 = 0.9
    scScore1 = 11
    scRepeats1 = 1
    scEarlyStop2 = 0.99
    scScore2 = 20
    scRepeats2 = 10
    threads = 8

    def __init__(self, args):
        args = getopt.getopt(args, self.shortParams, self.longParams)[0]
        print(f"Your arguments: posnegset {ProbeitUtils.getUserArgs(args)}")
        # PARSE PARAMETERS
        for opt, val in args:
            if opt in ('-h', '--help'):
                self.printUsage()
            try:
                self.inputGenome = str(val) if opt in ('-p', '--positive') else self.inputGenome
                self.negGenome = str(val) if opt in ('-n', '--negative') else self.negGenome
                self.workDir = str(val) if opt in ('-o', '--output') else self.workDir
                # optional args
                self.threads = int(val) if opt == '--threads' else self.threads
                self.windowSize = int(val) if opt == '--window-size' else self.windowSize
                self.doRemoveRedundancy = True if opt == '--remove-reduncancy' else self.doRemoveRedundancy
                self.doDesignProbe2 = False if opt == '--not-make-probe2' else self.doDesignProbe2
                self.doThermoFilter = False if opt == '--not-thermo-filter' else self.doThermoFilter
                self.isLigationProbe = True if opt == '--ligation-probe' else self.isLigationProbe
                self.pLen1 = int(val) if opt == '--probe1-len' else self.pLen1
                self.pLen2 = int(val) if opt == '--probe2-len' else self.pLen2
                self.cmError1 = int(val) if opt == '--probe1-error' else self.cmError1
                self.cmError2 = int(val) if opt == '--probe2-error' else self.cmError2
                self.scCoverage1 = int(val) if opt == '--probe1-cover' else self.scCoverage1
                self.scCoverage2 = int(val) if opt == '--probe2-cover' else self.scCoverage2
                self.scRepeats1 = int(val) if opt == '--probe1-repeat' else self.scRepeats1
                self.scRepeats2 = int(val) if opt == '--probe2-repeat' else self.scRepeats2
                self.scEarlyStop1 = float(val)/100 if opt == '--probe1-earlystop' else self.scEarlyStop1
                self.scEarlyStop2 = float(val)/100 if opt == '--probe2-earlystop' else self.scEarlyStop2
                # HIDDEN args
                # identity in mmseqs cluster
                self.cluIdentity = float(val) if opt == '--dedup-id' else self.cluIdentity
                # identity in mmeqs seqrch
                self.ridNegId = float(val) if opt == '--rid-neg-id' else self.ridNegId
                # setcover similarity
                self.scScore1 = int(val) if opt == '--probe1-dist' else self.scScore1
                self.scScore2 = int(val) if opt == '--probe2-dist' else self.scScore2
            except Exception as e:
                self.printUsage()
        
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
            print(f'[ERR] The directory named {self.workDir} already exsits.')
            isBadArguments = True
        if isBadArguments:
            self.printUsage()
        print('You input proper arguments.')
        
        # DIRECTORIES
        self.workDir = ProbeitUtils.defineDirectory(self.workDir.split(os.path.sep)[0], make=True)
        self.inputDir1 = ProbeitUtils.defineDirectory('input1', make=True, root=self.workDir)
        self.maskingDir = ProbeitUtils.defineDirectory('neg_filter', make=True, root=self.workDir)
        self.thermoFilteringDir = ProbeitUtils.defineDirectory('thermo_filter', make=True, root=self.workDir)
        self.cmDir1 = ProbeitUtils.defineDirectory('mapping_probe1', make=True, root=self.workDir)
        self.idxDir1 = ProbeitUtils.defineDirectory('index_probe1', make=False, root=self.workDir)
        self.scDir1 = ProbeitUtils.defineDirectory('setcover_probe1', make=True, root=self.workDir)
        if self.doDesignProbe2:
            self.inputDir2 = ProbeitUtils.defineDirectory('input2', make=True, root=self.workDir)
            self.cmDir2 = ProbeitUtils.defineDirectory('mapping_probe2', make=True, root=self.workDir)
            self.idxDir2 = ProbeitUtils.defineDirectory('index_probe2', make=False, root=self.workDir)
            self.scDir2 = ProbeitUtils.defineDirectory('setcover_probe2', make=True, root=self.workDir)

        # FILES 
        self.log = ProbeitUtils.defineFile(self.workDir, Config.log)
        self.lookup1 = ProbeitUtils.defineFile(self.workDir, 'genome.lookup')
        self.genome  = self.inputGenome 
        self.window1FASTA = ProbeitUtils.defineFile(self.inputDir1, Config.window)
        self.window1PosBED = ProbeitUtils.defineFile(self.inputDir1, 'window.bed') 
        self.posKmers = ProbeitUtils.defineFile(self.inputDir1, 'kmers.fa')
        self.negRemPosBED = ProbeitUtils.defineFile(self.maskingDir, 'negRemPos.bed')
        self.thermoImpPosBED = ProbeitUtils.defineFile(self.thermoFilteringDir, "thermo_neg.bed")
        self.thermoCalcTSV = ProbeitUtils.defineFile(self.thermoFilteringDir, "thermo_uncut.tsv")
        self.scPosBed1 = ProbeitUtils.defineFile(self.scDir1, 'result.bed')       
        self.tempProbe1 = ProbeitUtils.defineFile(self.workDir, 'temp1.fa')
        self.probe1 = ProbeitUtils.defineFile(self.workDir, Config.probe1)
        if self.doDesignProbe2:
            self.window2FASTA = ProbeitUtils.defineFile(self.inputDir2, Config.window)
            self.window2PosBED = ProbeitUtils.defineFile(self.inputDir1, 'window.bed') 
            self.lookup2 = ProbeitUtils.defineFile(self.workDir, 'probe1.lookup')
            self.tempProbe2 = ProbeitUtils.defineFile(self.workDir, 'temp2.fa')
            self.probe2 = ProbeitUtils.defineFile(self.workDir, Config.probe2)
        self.logUpdate('[INFO]Your arguments: ' + ' '.join(['{} {}'.format(i[0], i[1]).strip() for i in args]), False)
        
    def logUpdate(self, msg, doPrint=True):
        if doPrint:
            print(msg)
        with open(self.log, 'a') as w:
            w.write(msg + '\n')

    def thermoFilter(self):
        self.logUpdate('[INFO]filter probes with thermodynamic features') 
        probeLen = self.pLen1
        minGC = Config.getMinGC()
        maxGC = Config.getMaxGC()
        maxHomoDimerTm = Config.getMaxhomoDimerTm()
        maxHairpinTm = Config.getMaxhairpinTm()
        minProbeTm = Config.getMinProbeTm()
        problematicSeqs = ['A' * 5, 'T' * 5, 'C' * 5, 'G' * 5]

        # method for calculating gc ratio
        def getContentGC(oligo):
            gcs = oligo.count('G') + oligo.count('C')
            return gcs / len(oligo)

        # method for calculating complexity
        def hasLowComplexity(candidate_oligo):
            for s in problematicSeqs:
                if s in candidate_oligo:
                    return True
            return False
        self.logUpdate("filter probes based on primer3")
        # PRINTOUT CUT-OFF VALUES
        self.logUpdate("Minimum Tm: " + str(minProbeTm))
        self.logUpdate("Minimum GC percentage: " + str(minGC))
        self.logUpdate("Maximum GC percentage: " + str(maxGC))
        self.logUpdate("Homodimer maximum Tm: " + str(maxHomoDimerTm))
        self.logUpdate("Hairpin maximum Tm: " + str(maxHairpinTm))
        
        #  TO MAKE LIGATION PROBE CANDIDATES DF
        with open(self.posKmers) as f:
            identity = []
            posStart = []
            posEnd = []
            seqs = []
            rc = []
            p1 = []
            p2 = []
            for title, kmer in SimpleFastaParser(f):
                split_name = title.split('\t', 2)
                identity.append(split_name[0])
                if len(split_name) > 1:
                    this_posStart = int(split_name[1])
                else:
                    this_posStart = 0
                this_posEnd = this_posStart + probeLen
                posStart.append(this_posStart)
                posEnd.append(this_posEnd)
                seqs.append(kmer)
                rcKmer = str(reverse_complement(Seq(kmer)))
                rc.append(rcKmer)
                if self.isLigationProbe:
                    mid_pos = round(probeLen / 2)
                    p1.append(rcKmer[0:mid_pos])
                    p2.append(rcKmer[mid_pos:probeLen])
                else:
                    p1.append(rcKmer)
                    # p2.append(rcKmer)
        ligCands = pd.DataFrame(
            list(zip(identity, posStart, posEnd, seqs, rc, p1, p2)),
            columns=['id', 'chromStart', 'chromEnd', 'genome_segment', 'rc', 'p1', 'p2'],
        )
        self.logUpdate(str(len(ligCands)) + " ligation probe candidate sets inputted")
        
        # TO MAKE THERMO FEATURES DF
        thermos = pd.DataFrame(list(set(p1 + p2)), columns=['p'])
        self.logUpdate("There were " + str(len(thermos)) + " unique probes")
        thermos = thermos.assign(ultimate_base=thermos['p'].str[-1])
        thermos = thermos.assign(penultimate_base=thermos['p'].str[-2])
        self.logUpdate("Ultimate and penultimate bases assigned")
        thermos = thermos.assign(tm=np.vectorize(primer3.calcTm)(thermos['p']))
        thermos['hairpin_tm'] = list(map(primer3.calcHairpinTm, thermos['p']))
        thermos = thermos.assign(homodimer_tm=np.vectorize(primer3.calcHomodimerTm)(thermos['p']))
        thermos = thermos.assign(intrinsic_probs=np.vectorize(hasLowComplexity)(thermos['p']))
        self.logUpdate("Thermodynamic calculations complete")
        thermos = thermos.assign(GC_perc=np.vectorize(getContentGC)(thermos['p']))
        self.logUpdate("GC perc assigned")
       
        # TO MAKE ligProbsAndThermos
        ligCandsAndThermos = pd.merge(
            ligCands, thermos.add_prefix('p1_'), how='left', left_on='p1', right_on='p1_p'
        )
        ligCandsAndThermos = (
            pd.merge(ligCandsAndThermos, thermos.add_prefix('p2_'), how='left', left_on='p2', right_on='p2_p')
        )
        
        # FILTERING
        sorder = ['p1_hairpin_tm', 'p2_hairpin_tm', 'p1_homodimer_tm', 'p2_homodimer_tm']
        ascending = [True, True, True, True]
        query = ['p1_intrinsic_probs', f'p1_tm<{minProbeTm}', f'p1_GC_perc<{minGC}', f'p1_GC_perc>{maxGC}', f'p1_homodimer_tm>{maxHomoDimerTm}', f'p1_hairpin_tm>{maxHairpinTm}']
        query2 = ['p2_intrinsic_probs', f'p2_tm<{minProbeTm}', f'p2_GC_perc<{minGC}', f'p2_GC_perc>{maxGC}', f'p2_homodimer_tm>{maxHomoDimerTm}', f'p2_hairpin_tm>{maxHairpinTm}']
        query = query + query2 if self.isLigationProbe else query
        impLigProbs = ligCandsAndThermos.query(' | '.join(query))  
        self.logUpdate('Improper Ligation probe candidates are removed: ' + str(len(impLigProbs)))
        impLigProbs = impLigProbs.sort_values(by=sorder, ascending=ascending)
        impLigProbs.rename(columns={'id': 'chrom'})[['chrom', 'chromStart', 'chromEnd']].to_csv(self.thermoImpPosBED, index=False, sep='\t')
        ligCandsAndThermos.to_csv((self.thermoCalcTSV), sep='\t')

    def getImproperKmers(self):
        if not self.doThermoFilter:
            return set()
        temp = set()
        for line in open(self.thermoImpPosBED):
            for kmer in ProbeitUtils.parseKmers(line)[1:]:
                temp.add(kmer)
        return temp

    def removeRedundancy(self):
        self.logUpdate('[INFO]deduplicate positive fasta')
        clusteredGenome = ProbeitUtils.defineFile(self.inputDir1, 'dedup') 
        tempDir = ProbeitUtils.defineDirectory('temp', make=False, root=self.inputDir1)
        self.genome, msg = ProbeitUtils.clusterGenome(self.inputGenome, clusteredGenome, tempDir, self.cluIdentity, threads=self.threads)
        self.logUpdate(msg)

    def makeWindow1(self):
        ProbeitUtils.simplifyFastaHeaders(self.genome, self.window1FASTA)
        ProbeitUtils.makeLookup(self.window1FASTA, self.lookup1, self.window1PosBED)
        
    def extractKmers(self):
        self.logUpdate('[INFO]get {}mer probes from positive genome'.format(self.pLen1))
        f = open(self.window1FASTA)
        kmers = {}
        genomeIdx = 0
        self.logUpdate('[INFO] extracts all possible kmers from the input genome')
        for _, s in SimpleFastaParser(open(self.window1FASTA)):
            l = len(s)
            for pos in range(l-self.pLen1+1):
                name = f'{genomeIdx},{pos}'
                seq = s[pos:pos+self.pLen1].upper()
                kmers.setdefault(seq,[])
                kmers[seq].append(name)
            genomeIdx += 1
        w = open(self.posKmers, 'w')
        for s, names in kmers.items():
            rep = names[0]
            w.write(f'>{rep};{"|".join(names)}\n{s}\n')
        w.close()
        f.close()

    def negFilter(self):
        negKmerResult = ProbeitUtils.defineFile(self.maskingDir, 'search.tsv')
        negKmerPosBED = ProbeitUtils.defineFile(self.maskingDir, 'search.bed') 
        # REMOVE NEGATIVE K-MERS
        self.logUpdate('[INFO]remove probes found in negative genome')
        self.logUpdate(ProbeitUtils.ridNegKmers(self.posKmers, self.negGenome, negKmerResult, self.maskingDir, self.ridNegId, self.threads))
        
        # MAKE DEDUPLICATED POSITIVE GENOME COORDINATE FILE(BED)
        with open(self.window1FASTA) as f:
            with open(self.window1PosBED, 'w') as w:
                idx = 0
                for _, seq in SimpleFastaParser(f):
                    w.write(f'{idx}\t0\t{str(len(seq.strip()))}\n')
                    idx += 1
        
        # EXTRACT NEGATIVE REMOVED POSITIONs BED
        w = open(negKmerPosBED, 'w')
        for i in open(negKmerResult):
            kmers = ProbeitUtils.parseKmers(i)[1:]
            for kmer in kmers:
                g, s = [int(i) for i in kmer.split(',')]
                w.write(f'{g}\t{s}\t{s+40}\n')       
        w.close()
        self.logUpdate(ProbeitUtils.getSubtractedBed(self.window1PosBED, negKmerPosBED, self.negRemPosBED))
    
    def makeWindow2(self):
        def _getScKmer(scKmer, window):
            temp = scKmer.split(',')
            g = int(temp[0])
            s = int(temp[1])
            e = s + self.pLen1
            window[g] = window[g][0], window[g][1][:s]+'N'*40+window[g][1][e:]
            return g, s, e
        
        def _getWinSeq(kmer, lenList, header):
            g = kmer[0] 
            s = kmer[1]
            e = kmer[2]
            gLen = lenList[g]
            ws = s - self.windowSize if s>self.windowSize else 0
            we = e + self.windowSize 
            we = gLen if we>gLen else we
            genome, seq = window1Seqs[g]
            seq = seq[ws:we] 
            return f'>{header}:{genome}:{s}\n{seq}\n'
        
        window1Seqs =[(h, s) for h, s in SimpleFastaParser(open(self.window1FASTA))]
        window1Lens = [len(i[1]) for i in window1Seqs]
        probe1List = [[_getScKmer(kmer, window1Seqs) for kmer in ProbeitUtils.parseKmers(line)] for line in open(self.scPosBed1)]
        w = open(self.window2FASTA, 'w')
        for i, kmers in enumerate(probe1List):
            for j, kmer in enumerate(kmers):
                w.write(_getWinSeq(kmer, window1Lens, f'p1_{i}'))
        w.close()
        ProbeitUtils.makeLookup(self.window2FASTA, self.lookup2, self.window2PosBED)

    def easyComMap(self, uniqComMap):
        impKmers = self.getImproperKmers()
        w = open(uniqComMap, 'w')
        alphabets = {'A', 'C', 'G', 'T'}
        lenSeq = dict()
        idx = 0
        for h, s in SimpleFastaParser(open(self.window1FASTA)) :
            lenSeq[str(idx)] = len(s)
            idx += 1     
        negativeKmers = set()
        prevGenome = '0'
        prevEnd = 0
        for line in open(self.negRemPosBED):
            g, s, e = line.strip().split()
            s = int(s)
            e = int(e)
            if g == prevGenome:
                for i in range(prevEnd,s):
                    negativeKmers.add(f'{g},{i}')
                prevEnd = e
            else:
                for i in range(prevEnd, lenSeq[prevGenome]):
                    negativeKmers.add(f'{prevGenome},{i}')
                prevGenome = g
                prevEnd = e
                for i in range(0,s):
                    negativeKmers.add(f'{g},{i}')
        else:
            for i in range(prevEnd, lenSeq[prevGenome]):
                negativeKmers.add(f'{prevGenome},{i}')

        for h, s in SimpleFastaParser(open(self.posKmers)):
            kmers = set(ProbeitUtils.parseKmers(h))
            isThermoImproper = kmers&impKmers
            notOnlyATGC = set(s) != alphabets&set(s)
            # TODO
            # isNegative = kmers&negativeKmers
            isNegative = set([ProbeitUtils.parseKmers(h)[0]])&negativeKmers
            if isNegative or isThermoImproper or notOnlyATGC:
                continue
            w.write(h+'\n')
        w.close()

    def makeProbe1(self):
        uniqComMap = ProbeitUtils.defineFile(self.cmDir1, 'uniq.genmap.csv')
        self.easyComMap(uniqComMap)
        msg = ProbeitUtils.makeProbe(
            self.tempProbe1, self.scPosBed1, self.window1FASTA, self.lookup1, self.pLen1, 
            self.scCoverage1, self.scEarlyStop1, self.scScore1, self.scRepeats1, uniqComMap,  
        )
        self.logUpdate(msg, False)
        with open(self.lookup1)as f:
            genomeKeys = {int(i.split()[0]):i.split()[1] for i in f}
        probeIdx = 0
        # probeNames1 = {}
        w = open(self.probe1, 'w')
        for h, s in sorted([(h,s) for h,s in SimpleFastaParser(open(self.tempProbe1))]):
            keys, pos = ProbeitUtils.parseGenmapPattern(h)
            keys = [genomeKeys[k] for k in keys]
            probes = ['{}:{}:{}'.format(keys[i], pos[i], pos[i] + self.pLen1) for i in range(len(keys))]
            w.write('>p1_{}\t{}\n{}\n'.format(probeIdx, ';'.join(probes), s))
            probeIdx += 1
        w.close()

    def makeProbe2(self):
        uniqComMap = ProbeitUtils.defineFile(self.cmDir2, 'uniq.genmap.csv')
        if self.cmError2 == 0:
            ProbeitUtils.simpleComputeMappability(self.window2FASTA, self.lookup2, self.window2PosBED, uniqComMap, self.pLen2)
        else:
            ProbeitUtils.computeMappability(self.window2FASTA, self.idxDir2, self.cmError2, self.pLen2, self.cmDir2, uniqComMap, threads=self.threads)
        scBed = ProbeitUtils.defineFile(self.scDir2, 'result.bed')
        msg = ProbeitUtils.makeProbe(
            self.tempProbe2, scBed, self.window2FASTA, self.lookup2, self.pLen2,
            self.scCoverage2, self.scEarlyStop2, self.scScore2, self.scRepeats2,  uniqComMap, 
            overlap=True
        )
        self.logUpdate(msg, False)
        with open(self.lookup2)as f:
            probe1Keys = {int(i.split()[0]):i.split()[1] for i in f}
        probeIdx = 0
        w = open(self.probe2, 'w')
        for h, s in sorted([(h,s) for h, s in SimpleFastaParser(open(self.tempProbe2))]):
            keys, _ = ProbeitUtils.parseGenmapPattern(h)           
            keys = set([probe1Keys[k] for k in keys])
            w.write('>p2_{}\t{}\n{}\n'.format(probeIdx, ';'.join([p for p in keys]), s) )
            probeIdx += 1
        w.close()
    
    def run(self):
        # MAKE 1st probe
        self.logUpdate("[INFO] make 1st probes")
        if self.doRemoveRedundancy:
            self.removeRedundancy()
        ProbeitUtils.sortFasta(self.genome)
        self.makeWindow1()
        self.extractKmers()
        if self.negGenome:
            self.negFilter()
        if self.doThermoFilter:
            self.thermoFilter()
        self.makeProbe1()
        # MAKE 2nd PROBES
        if self.doDesignProbe2:
            self.logUpdate("[INFO] make 2nd probes")
            self.makeWindow2()
            self.makeProbe2()
        print('COMPLETE!!!')
        return

    @staticmethod
    def printUsage():
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
        print(" --remove-reduncancy NONE")
        print("\t Use it when you NEED to cluster positive genome")
        print(" --not-make-probe2 NONE")
        print("\t Use it when you DO NOT need to make 2nd probes")
        print(" --not-thermo-filter NONE")
        print("\t Use it when you DO NOT need the thermodynamic filter")
        print(" --ligation-probe NONE")
        print("\t Use it when you want to make ligational probes as probe1.")        
        print(" --probe1-len INT[40]")
        print("\t Length of 1st Probes")
        print(" --probe2-len INT[20]")
        print("\t Length of 2nd Probes")

        print("ADDITIONAL OPTIONS FOR GENMAP PART: Genmap calculates mappability by summarizing k-mers.")
        print(" --probe1-error INT[0]")
        print("\t The number of error allowed in 1st Probes")
        print(" --probe2-error INT[1]")
        print("\t The number of error allowed in 2nd Probes")

        print("ADDITIONAL OPTIONS FOR SETCOVER PART: Setcover makes probesets cover targets with the minimum number of probes.")
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


class SNP:
    args = []
    shortParams = 'hr:a:s:p:m:o:'
    longParams = [
        # usage
        'help',
        # required
        'reference=', 'strain=', 'positions=', 'mutations=', 'output=', 'annotation=',
        # optional
        'not-make-probe2', 'max-window', 'window-size=', 'probe1-len=', 'probe2-len=', 'threads=',
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
    probLen1 = 40
    probeLen2 = 20
    scCoverage2 = 1
    scEarlyStop2 = 0.99
    scScore2 = 20
    scRepeats2 = 10
    cmError2 = 1
    searchKmer = 12
    threads = 8


    def __init__(self, args):
        args = getopt.getopt(args, self.shortParams, self.longParams)[0]
        print(f"Your arguments: posnegset {ProbeitUtils.getUserArgs(args)}")
        # PARSE PARAMETERS
        for opt, val in args:
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
                self.needProbe2 = False if opt == '--not-make-probe2' else self.needProbe2
                self.threads = int(val) if opt == '--threads' else self.threads
                self.isMaxWindow = True if opt == '--max-window' else self.isMaxWindow
                self.windowSize = int(val) if opt == '--window-size' else self.windowSize
                self.probLen1 = int(val) if opt == '--probe1-len' else self.probLen1
                self.probeLen2 = int(val) if opt == '--probe2-len' else self.probeLen2
                self.cmError2 = int(val) if opt == '--probe2-error' else self.cmError2
                self.scCoverage2 = int(val) if opt == '--probe2-cover' else self.scCoverage2
                self.scEarlyStop2 = float(val)/100 if opt == '--probe2-earlystop' else self.scEarlyStop2
                self.scRepeats2 = int(val) if opt == '--probe2-repeat' else self.scRepeats2
                # hidden args
                self.scScore2 = int(val) if opt == '--probe2-dist' else self.scScore2
                self.searchKmer = int(val) if opt == '--search-kmer' else self.searchKmer
            except Exception as e:
                print(e)
                print("Your arguments: snp {}".format(ProbeitUtils.getUserArgs(args)))
                self.printUsage() 
        self.snpList = list(set(self.snpList))
        validPosList = []
        for p in self.posList:
            if 0 < p <= self.probLen1:
                validPosList.append(int(p))
            else:
                print('[ERROR] {} is not proper for position list.')
        self.posList = validPosList
        
        print("Your arguments: {}".format('snp ' + ProbeitUtils.getUserArgs(args)))
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
            print(message.format('[ERROR]', '--mutaions'))
            isBadArguments = True
        if not self.workDir:
            print(message.format('[ERROR]', '--output'))
            isBadArguments = True
        if os.path.exists(self.workDir):
            print(f'[ERR] The directory named {self.workDir} already exsits.')
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
            self.idxDir2 = ProbeitUtils.defineDirectory('index', make=False, root=self.workDir)
            self.cmDir2 = ProbeitUtils.defineDirectory('mapping2', make=True, root=self.workDir)
            self.scDir2 = ProbeitUtils.defineDirectory('setcover2', make=True, root=self.workDir)
        # FILES
        self.log = ProbeitUtils.defineFile(self.workDir, 'log.txt')
        self.posProbeCSVs = {pos: ProbeitUtils.defineFile(self.inputDir1, f'pos{pos}.csv') for pos in self.posList}
        self.posProbeCSVs[-1] = ProbeitUtils.defineFile(self.inputDir1, 'merged.csv')
        self.probe1 = ProbeitUtils.defineFile(self.workDir, Config.probe1)
        if self.needProbe2:
            self.window1FASTA = ProbeitUtils.defineFile(self.inputDir1, Config.window)
            self.lookup = ProbeitUtils.defineFile(self.workDir, 'probe1.lookup')
            self.window1PosBED = ProbeitUtils.defineFile(self.inputDir1, 'window1.bed')
            self.window1TSV = ProbeitUtils.defineFile(self.inputDir1, 'window1.tsv')
            self.maskedGenomeFASTA= ProbeitUtils.defineFile(self.inputDir2, 'masked.fa')
            self.window2PosBED = ProbeitUtils.defineFile(self.inputDir2, 'window2.bed')
            self.window2FASTA = ProbeitUtils.defineFile(self.inputDir2, Config.window)
            self.tempProbe2 = ProbeitUtils.defineFile(self.workDir, 'temp2.fa')
            self.probe2 = ProbeitUtils.defineFile(self.workDir, Config.probe2)
        self.logUpdate('[INFO]Your arguments: snp ' + ProbeitUtils.getUserArgs(args), False)
        # VARIABLES
        self.probesByPos = {pos: [] for pos in self.posList + [-1]}        

    @staticmethod
    def getArgList(value, isInt=False):
        if isInt:
            return [int(i.strip()) for i in value.split(',')]
        else:
            return value.split(',')

    def logUpdate(self, msg, doPrint=True):
        if doPrint:
            print(msg)
        with open(self.log, 'a') as w:
            w.write(msg+'\n')

    def _searchSnpFromStrGenome(self, mutation, seqWithSNP):
        snpKmers = ProbeitUtils.defineFile(self.searchDir, f'search_{mutation}.tsv')
        searchProbe = ProbeitUtils.defineFile(self.searchDir, f'search_{mutation}.fa')
        with open(searchProbe, 'w') as w:
            w.write('>{}\n{}\n'.format(mutation, seqWithSNP))
        self.logUpdate(
            ProbeitUtils.searchSNPs(self.searchDir, searchProbe, self.strGenome, snpKmers, self.searchKmer, threads=self.threads)
        )
        return snpKmers

    @staticmethod
    def _parseSearcgResult(searchResult):
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

    def _trimSearchResult(self, result, wtSeq, muttype, maxPos, snp, startPos=-1, wtCodon=''):
        result = pd.read_csv(result, sep='\t', header=None)
        if muttype == 'aa':
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
        elif muttype == 'nt':
            print('NT')
            nt2 = snp[-1]
            result.columns = ['subGenome', 'SNPbyNT', 'match', 'STsequence']
            result['WTsequence'] = wtSeq
            result['locSNP'] = maxPos - 1
            result = result[result.STsequence.apply(lambda x: len(x) == len(wtSeq))]
            if result.empty:
                return result
            result = result[result.STsequence.apply(lambda x: x[maxPos - 1] == nt2)]  
            if result.empty:
                return result
            return result
        else:
            return pd.DataFrame()

    def findSNPs(self): 
        def _findSNPsAA(orf, mutation):
            # DO NOT EXIST REFERENCE GENOME ANNOTATION
            if self.refGenomeAnnot == '':
                self.logUpdate('[warn]For Amino Acid based SNPs reference annotation needed.')
                return  -1, None
            
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
            seqWithSNP = refSeq[codonStartPos - (maxPos-1): codonEndPos + (self.probLen1 - minPos)]
            searchResult = self._searchSnpFromStrGenome(mutation, seqWithSNP)
            trimmedResult = self._trimSearchResult(searchResult, seqWithSNP, mutType, maxPos, snp, codonStartPos, refCodon)
            if trimmedResult.empty:
                self.logUpdate(f'[WARN] Problems occured searching snp {snp} in the strain genome or reference genome')
                return -1, None
            
            # RETURN RESULT
            wtSequence, stSequence, ntSNP, locSnp, found = self._parseSearcgResult(searchResult=trimmedResult)
            self.logUpdate('[INFO]aa:{}:{} converted to nt:{}'.format(orf, mutation, ntSNP))
            aaOrf = '{}:{}'.format(orf, mutation)
            return found, ParaSeqs(ntSNP, aaOrf, wtSequence, stSequence, mutLoc=locSnp, probLen=self.probLen1)

        def _findSNPsNT(orf, mutation):
            nt1, nt2, snpPos = mutation[0], mutation[-1], int(mutation[1:-1])
            refNT = refSeq[snpPos]
            # NOT A VALID SNP
            if nt1 != refNT:
                self.logUpdate('[warn]Failure to find SNP {} in reference genome'.format(snp))
                return -1, None
            
            #  SEARCH
            seqWithSNP = refSeq[snpPos - (maxPos - 1):snpPos + 1 + (self.probLen1 - minPos)]
            searchResult = self._searchSnpFromStrGenome(mutation, seqWithSNP)
            trimmedResult = self._trimSearchResult(searchResult, seqWithSNP, mutType, maxPos, snp)
            if trimmedResult.empty:
                self.logUpdate(f'[WARN] Problems occured searching snp {snp} in the strain genome or reference genome')
                return -1, None
            
            # RETURN RESULT
            wtSequence, stSequence, ntSNP, locSnp, found = self._parseSearcgResult(searchResult=trimmedResult)
            return found, ParaSeqs(ntSNP, '', wtSequence, stSequence, mutLoc=locSnp, probLen=self.probLen1)

        minPos = min(self.posList)
        maxPos = max(self.posList)
        refSeq = self.getReferenceSeq(self.refGenome)
        if not refSeq:
            self.logUpdate('[warn]Failure to get reference sequence from reference genome.')
            self.printUsage()
        for snp in self.snpList:
            self.logUpdate('[INFO]SNP {}'.format(snp))
            mutType, orf, mutation = self.parseMutation(snp)
            found, mutSeqs = _findSNPsAA(orf, mutation) if mutType == 'aa' else _findSNPsNT(orf, mutation)
            if found <= 0 or not found:
                self.logUpdate('[WARN] Failure to find SNP {} in strain genome'.format(snp))
                continue
            self.probesByPos[-1].append(mutSeqs)
            for pos in self.posList:
                wtProbe, stProbe = mutSeqs.getProbesWithPos(pos)
                paraSeq = ParaSeqs(mutSeqs.ntSnp, mutSeqs.aaSnp, wtProbe, stProbe, found=found, probLen=self.probLen1)
                self.probesByPos[pos].append(paraSeq)

    @staticmethod
    def makeMaskedPosBED(inputDF, outputBed):
        inputDF['seqID'] = inputDF['seqID'].apply(lambda x: x.split(';')[0])
        inputDF.to_csv(outputBed, sep='\t', header=False, index=False)

    def makeProbe1(self):
        probeLines = []
        for pos in self.posList + [-1]:
            probeCSV = self.posProbeCSVs[pos]
            csvWriter = open(probeCSV, 'w')
            csvWriter.write('WT sequence,ST sequence,found,ntSNP,aaSNP\n')
            for p in self.probesByPos[pos]:
                csvWriter.write(f'{p.wtSeq},{p.stSeq},{p.found},{p.ntSnp},{p.aaSnp}\n')
                if pos != -1:
                    probeLines.append(f'>{p.ntSnp}{"=" + p.aaSnp if p.aaSnp else ""};{pos}\n{p.stSeq}\n')
            csvWriter.close()
        if not probeLines:
            self.logUpdate('[ERROR] Cannot find any SNP in strain genomes')
            self.printUsage()
        with open(self.probe1, 'w') as w:
            w.writelines(sorted(probeLines))
        
    # TODO
    def makeWindows(self):
        # MAKE WINDOW1
        window1 = [f'>snp1_{i}_{p.ntSnp}{ "=" + p.aaSnp if p.aaSnp else ""}\n{p.stSeq}\n' for i, p in enumerate(self.probesByPos[-1])]           
        with open(self.window1FASTA, 'w') as w:
            w.writelines(window1)
        
        # MAKE WINDOW 2
        msg = ProbeitUtils.getPatternPosition(self.window1FASTA, self.strGenome, self.window1TSV)
        self.logUpdate(msg, False)
        self.makeMaskedPosBED(pd.read_csv(self.window1TSV, sep='\t')[['seqID', 'start', 'end']], self.window1PosBED)
        ProbeitUtils.getWindowFasta(
            self.strGenome, self.window1PosBED, self.maskedGenomeFASTA, self.window2PosBED, self.window2FASTA, self.windowSize
        )
        ProbeitUtils.makeLookup(self.window2FASTA, self.lookup, self.window2PosBED)

    def makeProbe2(self):
        # COMPUTE MAPPABILITY
        uniqComMap = ProbeitUtils.defineFile(self.cmDir2, 'uniq.genmap.csv')
        message = '[INFO] compute mappability\n'
        if self.cmError2 == 0:
            ProbeitUtils.simpleComputeMappability(self.window2FASTA, self.lookup, self.window2PosBED, uniqComMap, self.probeLen2)
        else:
            message += ProbeitUtils.computeMappability(self.window2FASTA, self.idxDir2, self.cmError2, self.probeLen2, self.cmDir2, uniqComMap, threads=self.threads)
        scBed = ProbeitUtils.defineFile(self.scDir2, 'result.bed')
        msg = ProbeitUtils.makeProbe(self.tempProbe2, scBed, self.window2FASTA, self.lookup, self.probeLen2, self.scCoverage2, self.scEarlyStop2, self.scScore2, self.scRepeats2, uniqComMap, overlap=True)
        self.logUpdate(msg, False)
        maskDF = pd.read_csv(self.window1TSV, sep='\t')
        kmers = list(maskDF['patternName'])  
        probe2lines = []
        for h, s in SimpleFastaParser(open(self.tempProbe2)):
            kmerIndex = ProbeitUtils.parseGenmapPattern(h)[0]
            coveredSNPs = list((set([kmers[i] for i in kmerIndex])))
            probe2lines.append(f'>{":".join(coveredSNPs)}\n{s}\n')
        with open(self.probe2, 'w') as w:
            w.writelines(sorted(probe2lines))
    
    def run(self):
        # MAKE 1st PROBEs
        self.logUpdate("[INFO]make 1st probes")
        self.findSNPs()
        self.makeProbe1()
        if self.needProbe2:
            # MAKE 2ND PROBEs
            self.logUpdate("[INFO]make 2nd probes")
            self.makeWindows()
        self.makeProbe2()
        print('COMPLETE!!!')
        return

    @staticmethod
    def printUsage():
        print("Probeit snp")
        print("It generates a probe set which detect input amino acid SNPs from strain genome.")
        print("probeit snp -r [REF GENOME] -s [STR GENOME] -p [positions] -m [SNPs] -o [DIR] -a [REF ANNOTATION]")

        print("Usage")
        print(" -h|--help NONE")
        print('\t Show usage')

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

