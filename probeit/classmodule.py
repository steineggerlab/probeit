#!/usr/bin/env python
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
    def delDir(cls, dir):
        if os.path.isdir(dir):
            print(f'[INFO] The directory named {dir} is removed not.')
            shutil.rmtree(dir)
        else:
            print(f'[INFO] The directory named {dir} does not exist.')

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
    def computeMappability(cls, genome, indexDir, error, len, outputDir, selCoords='', thread=8):
        selector = f'-S {selCoords}' if selCoords else ''
        cls.delDir(indexDir)
        command0 = " > /dev/null"
        command1 = "genmap index -F {} -I {}".format(genome, indexDir)
        cls.runCommand(command1 + command0)
        command2 = "genmap map --no-reverse-complement -E {} {} --csv -K {} -t -b --frequency-large -I {} -O {} -T {}"
        command2 = command2.format(error, selector, len, indexDir, outputDir, thread)
        cls.runCommand(command2 + command0)
        msg = '[CLI] {}{}\n[CLI] {}{}\n'.format(command1, command0, command2, command0)
        return cls.defineFile(outputDir, f'{cls.getFileName(genome, False)}.genmap.csv'), msg

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
    def makeLookup(cls, windowFasta, lookup):
        with open(windowFasta) as f:
            headers = [title.strip() for title, _ in SimpleFastaParser(f)]
        lookupLines = [f'{i}\t{header}\n' for i, header in enumerate(headers)]
        with open(lookup, 'w') as w:
            w.writelines(lookupLines)
        return lookup

    @classmethod
    def searchSNPs(cls, workDir, inputFasta, strGenomeFasta, result, kmer=12, thread=8):
        searchDir = cls.defineDirectory('search', root=workDir, make=True)
        tempDir = cls.defineDirectory('temp', root=searchDir, make=False)
        searchdb = cls.defineFile(searchDir, 'searchDB' )
        strdb =  cls.defineFile(searchDir, 'strainDB')
        aln = cls.defineFile(searchDir , 'mmseqs.aln')
        cmd0 = ' --threads {}'.format(thread)
        cmd1 = 'mmseqs createdb {} {}'
        cmd2 = 'mmseqs createdb {} {}'
        cmd3 = 'mmseqs search {} {} {} {} --search-type 3 -k {}'
        cmd4 = 'mmseqs convertalis {} {} {} {} --format-output target,query,tseq,tstart,tend --search-type 3'
        out1, err1 = cls.runCommand(cmd1.format(inputFasta, searchdb), verbose=True)
        out2, err2 = cls.runCommand(cmd2.format(strGenomeFasta, strdb), verbose=True)
        out3, err3 = cls.runCommand(cmd3.format(searchdb, strdb, aln, tempDir, kmer, thread) + cmd0, verbose=True)
        out4, err4 = cls.runCommand(cmd4.format(searchdb, strdb, aln, result, thread) + cmd0, verbose=True)
        df = pd.read_csv(result, sep='\t', header=None)
        df.columns = ['substr', 'snp', 'strseq', 'start', 'end']
        df['aln'] = df.apply(lambda x: x[2][int(x[3]-1):int(x[4])], axis=1)
        df['len'] = df.aln.apply(lambda x: len(x)-1)
        df = df[['substr', 'snp', 'len', 'aln']]
        df.to_csv(result, header=False, index=False, sep='\t')
        print(err1 + err2 + err3 + err4)
        ProbeitUtils.delDir(searchDir)
        return out1 + out2 + out3 + out4

    @staticmethod
    def copyFile(original, copy):
            shutil.copy(original, copy)

    @classmethod
    def clusterGenome(cls, input, output, dir, seqIdentity, cluster=True, thread=8):
        if not cluster:
            output = output+'_rep_seq.fasta'
            cls.copyFile(input, output)
            return output, '', ''
        command = (
            f'mmseqs easy-linclust {input} {output} {dir} -v 3 --kmer-per-seq-scale 0.5 --kmer-per-seq 1000 ' +
            f'--min-seq-id {seqIdentity} --cov-mode 1 -c 0.95 --remove-tmp-files 0 --threads {thread}' 
        )
        stdout, stderr = cls.runCommand(command, verbose=True)
        return output+'_rep_seq.fasta',  stdout, stderr

    @classmethod
    def ridNegKmers(cls, output, negative, maskOutput, outputDir, seqInProbe, thread=8):
        tempDir = cls.defineDirectory("tmp", root=outputDir)
        command = (
            f'mmseqs easy-search {output} {negative} {maskOutput} {tempDir} -v 3 --spaced-kmer-mode 0 -k 13 --mask 0 ' +
            f'--format-output query,target,pident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue ' +
            f'-c 0.9 --min-seq-id {seqInProbe} --cov-mode 2 --alignment-mode 4 --search-type 3 --threads {thread}'
        ) 
        stdOut, stdErr = cls.runCommand(command, verbose=True)
        print(stdOut)
        print(stdErr)
        return stdOut, stdErr

    @classmethod
    def parseKmers(cls,line):
        return re.findall(r'[0-9]+,[0-9]+', line)

    @classmethod
    def simplifyCMResult(cls, inputCSV, outputCSV, thermoImproperProbes=[]):        
        print(inputCSV)
        print('[INFO] Removing duplication from result of coumpute mappability.')
        w = open(outputCSV, 'w')
        with open(inputCSV) as f:
            found = set()
            for line in f:
                kmers = cls.parseKmers(line)
                if not kmers:
                    continue
                head = kmers[0]
                tail = '|'.join(kmers[1:])
                if not set(kmers)&found and head not in thermoImproperProbes:
                    w.write(f'{head};{tail}\n')
                found.add(kmers[0])
        w.close()

    @classmethod
    def setCover(cls, coverage, length, eStop, dist, reps, mapCSV, genome):
        filePath = os.path.sep.join(os.path.realpath(__file__).split(os.path.sep)[:-1])
        setcoverPath = '{}{}{}{}{}'.format(filePath, os.path.sep, 'setcover', os.path.sep, 'setcover')
        setcoverPath = setcoverPath   if os.path.exists(setcoverPath) else 'setcover'
        command = " -c {} -l {} -p {} -d {} -i {} {} {}".format(coverage, length, eStop, dist, reps, mapCSV, genome)
        stdOut, stdErr = cls.runCommand(setcoverPath + command, verbose=True)
        return stdOut, stdErr

    @classmethod
    def makeMinzProbeBed(cls, lookup, setcoverResult, setcoverResultBed, probeLen):
        genomeAndIdx = dict()
        with open(lookup) as f1:
            for line in f1:
                idx = int(line.split()[0].strip())
                genome = line.split()[1].strip()
                genomeAndIdx[idx] = genome
        with open(setcoverResult) as f2:
            with open(setcoverResultBed, 'w') as w:
                for line in f2:
                    matchedKmers = line.split(';')[1].strip()
                    idx = int(line.split(';')[0].split(',')[0])
                    pos = line.split(';')[0].split(',')[1]
                    w.write('\t'.join([genomeAndIdx[idx], pos, str(int(pos) + probeLen), matchedKmers]) + '\n')

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
            cls, output, window, lookup, indexDir, cmDir, scDir, probeLen, 
            cmError, scCoverage, scEarlyStop, scScore, scRepeats, 
            selCoords='', overlap=False, thread=8, improperKmers=[]
    ):
        scLen = 1 if overlap else probeLen 
        minzResult = cls.defineFile(scDir, 'result')
        minzBed = cls.defineFile(scDir, 'result.bed')
        uniqComMap = cls.defineFile(cmDir, 'uniq.genmap.csv')
        # COMPUTEMAPPABILITY
        message = ''
        message += '[INFO]compute mappability\n'
        comMap, msg = ProbeitUtils.computeMappability(window, indexDir, cmError, probeLen, cmDir, selCoords=selCoords, thread=thread)
        ProbeitUtils
        message += msg
        message += '[INFO]deduplicate genmap result\n'
        ProbeitUtils.simplifyCMResult(comMap, uniqComMap, improperKmers)
        message += "[INFO]minimize probe set\n"
        message += f"[CLI] setcover -c {scCoverage} -l {scLen} -p {scEarlyStop} -d {scScore} -i {scRepeats} {uniqComMap} {window}\n"
        # SETCOVER
        msg, err = ProbeitUtils.setCover(scCoverage, scLen, scEarlyStop, scScore, scRepeats, uniqComMap, window)
        message += msg
        message += err
        with open(minzResult, 'w') as w:
            w.write(msg)
        # MAKE PROBEs
        ProbeitUtils.makeMinzProbeBed(lookup, minzResult, minzBed, probeLen)
        message += ProbeitUtils.getSubseqFasta(minzBed, window, output)
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
            self.subWork = PosNegSet(self.args[1:])
            self.subWork.execute()
            return
        elif self.args[0] == 'snp':
            print('CURRENT: ', os.getcwd())
            self.subWork = SNP(self.args[1:])
            self.subWork.execute()
            return
        else:
            self.printUsage()
            return

    @staticmethod
    def printUsage():
        print("PROBEIT")
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
        'not-make-probe2', 'not-cluster', 'not-thermo-filter', 
        # 
        'probe1-len=', 'probe2-len=', 'window-size=', 'threads=',
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
    needCluster = True
    needProbe2 = True
    needThermoFilter = True
    windowSize = 200
    scEarlyStop1 = 0.9
    scScore1 = 11
    scRepeats1 = 1
    scEarlyStop2 = 0.99
    scScore2 = 20
    scRepeats2 = 10
    threads = 8

    # FILE NAMES
    PROBE1 = 'probe1.fa'
    PROBE2 = 'probe2.fa'
    TEMP1 = 'temp1.fa'
    TEMP2 = 'temp2.fa'
    # CLUSTER = ''
    # NEGATIVE_FILTER = ''
    # THERMO_IMPROPER = ''
    # THERMO_INFO = ''
    # GENMAP RELATED

    # SETCOVER RELATED
    SETCOVER_RESULT = 'result'
    SETCOVER_COORDS = 'result.bed'



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
                self.windowSize - int(val) if opt == '--window-size' else self.windowSize
                self.needCluster = False if opt == '--not-cluster' else self.needCluster
                self.needProbe2 = False if opt == '--not-make-probe2' else self.needProbe2
                self.needThermoFilter = False if opt == '--not-thermo-filter' else self.needThermoFilter
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
        self.inputDir = ProbeitUtils.defineDirectory('input', make=True, root=self.workDir)
        self.dedupDir = ProbeitUtils.defineDirectory('cluster', make=True, root=self.workDir)
        self.maskingDir = ProbeitUtils.defineDirectory('mmseqs', make=True, root=self.workDir)
        self.thermoFilteringDir = ProbeitUtils.defineDirectory('filter', make=True, root=self.workDir)
        self.cmDir1 = ProbeitUtils.defineDirectory('mapping_probe1', make=True, root=self.workDir)
        self.idxDir1 = ProbeitUtils.defineDirectory('index_probe1', make=False, root=self.workDir)
        self.scDir1 = ProbeitUtils.defineDirectory('setcover_probe1', make=True, root=self.workDir)
        if self.needProbe2:
            self.inputDir2 = ProbeitUtils.defineDirectory('input_probe2', make=True, root=self.workDir)
            self.cmDir2 = ProbeitUtils.defineDirectory('mapping_probe2', make=True, root=self.workDir)
            self.idxDir2 = ProbeitUtils.defineDirectory('index_probe2', make=False, root=self.workDir)
            self.scDir2 = ProbeitUtils.defineDirectory('setcover_probe2', make=True, root=self.workDir)

        # FILES 
        self.log = ProbeitUtils.defineFile(self.workDir, 'log.txt')
        self.lookup1 = ProbeitUtils.defineFile(self.workDir, 'genome.lookup')
        self.posKmers = ProbeitUtils.defineFile(self.maskingDir, 'probes.fa')
        self.deDupGenome  = ProbeitUtils.defineFile(self.dedupDir, 'dedup')  
        self.window1 = ProbeitUtils.defineFile(self.dedupDir, 'seq.fasta')
        self.negRemCoords = ProbeitUtils.defineFile(self.cmDir1, 'crosstaxa.bed')
        self.thermoImpCoords = ProbeitUtils.defineFile(self.thermoFilteringDir, "thermo_neg.bed")
        self.thermoCalcTSV = ProbeitUtils.defineFile(self.thermoFilteringDir, "thermo_uncut.tsv")
        self.kmersForTheromFilter = ProbeitUtils.defineFile(self.thermoFilteringDir, 'probes.fa')
        self.scResult1 = ProbeitUtils.defineFile(self.scDir1, self.SETCOVER_RESULT)
        self.setcoverCoords1 = ProbeitUtils.defineFile(self.scDir1, self.SETCOVER_COORDS)       
        self.tempProbe1 = ProbeitUtils.defineFile(self.workDir, self.TEMP1)
        self.probe1 = ProbeitUtils.defineFile(self.workDir, 'probe1.fa')
        if self.needProbe2:
            self.window2 = ProbeitUtils.defineFile(self.inputDir2,'seq.fa')
            self.lookup2 = ProbeitUtils.defineFile(self.workDir, 'probe1.lookup')
            self.tempProbe2 = ProbeitUtils.defineFile(self.workDir, self.TEMP2)
            self.probe2 = ProbeitUtils.defineFile(self.workDir, 'probe2.fa')
        self.logUpdate('[INFO]Your arguments: ' + ' '.join(['{} {}'.format(i[0], i[1]).strip() for i in args]), False)
        
    def logUpdate(self, msg, doPrint=True):
        if doPrint:
            print(msg)
        with open(self.log, 'a') as w:
            w.write(msg + '\n')

    def copyFile(self, original, copy):
        shutil.copy(original, copy)
        self.logUpdate('[INFO]{} is copied to {}'.format(original, copy))
        return copy

    def thermoFilter(self):
        probeLen = self.pLen1
        minGC = 0.30
        maxGC = 0.70
        homoDimerTmCutoff = 60
        hairpinTmMax = 60
        minProbeTm = 40
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
        self.logUpdate("Homodimer maximum Tm: " + str(homoDimerTmCutoff))
        self.logUpdate("Hairpin maximum Tm: " + str(hairpinTmMax))
        
        #  TO MAKE LIGATION PROBE CANDIDATES DF
        with open(self.kmersForTheromFilter) as f:
            identity = []
            posStart = []
            posEnd = []
            seqs = []
            rc = []
            p1 = []
            p2 = []
            for title, sequence in SimpleFastaParser(f):
                split_name = title.split('\t', 2)
                identity.append(split_name[0])
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
        impLigProbs = (
            ligCandsAndThermos.query(
                # intrinsic sequence
                'p1_intrinsic_probs | p2_intrinsic_probs | ' +
                # TM
                f'p1_tm<{minProbeTm} | p2_tm<{minProbeTm} | ' +
                # GC content
                f'p1_GC_perc<{minGC} | p2_GC_perc<{minGC} | p1_GC_perc>{maxGC} | p2_GC_perc>{maxGC} | ' +
                # homodimer TM
                f'p1_homodimer_tm>{homoDimerTmCutoff} | p2_homodimer_tm>{homoDimerTmCutoff} | ' +
                # hairpin TM
                f'p1_hairpin_tm>{hairpinTmMax} | p2_hairpin_tm>{hairpinTmMax}'
            ) 
        ) 
        self.logUpdate('Improper Ligation probe candidates are removed: ' + str(len(impLigProbs)))
        impLigProbs = impLigProbs.sort_values(by=sorder, ascending=ascending)
        impLigProbs.rename(columns={'id': 'chrom'})[['chrom', 'chromStart', 'chromEnd']].to_csv(self.thermoImpCoords, index=False, sep='\t')
        ligCandsAndThermos.to_csv((self.thermoCalcTSV), sep='\t')

    def getImproperKmers(self):
        if not self.needThermoFilter:
            return []
        genomeKeys = {i.split()[1]: int(i.split()[0]) for i in open(self.lookup1)}
        df = pd.read_csv(self.thermoImpCoords , sep='\t')
        df.chrom = df.chrom.apply(lambda x: genomeKeys[x])
        return list(df.apply(lambda x: f'{x[0]},{x[1]}', axis=1))

    def remRedundancy(self):
        self.logUpdate('[INFO]deduplicate positive fasta')
        tempDir = ProbeitUtils.defineDirectory('temp', make=False, root=self.dedupDir)
        self.deDupGenome, msg, err = ProbeitUtils.clusterGenome(self.inputGenome, self.deDupGenome, tempDir, self.cluIdentity, cluster=self.needCluster, thread=self.threads)
        self.logUpdate(msg + err)
        ProbeitUtils.sortFasta(self.deDupGenome)

    def makeWindow1(self):
        ProbeitUtils.simplifyFastaHeaders(self.deDupGenome, self.window1)
        ProbeitUtils.makeLookup(self.window1, self.lookup1)
        
    def extractKmersFromInputGenome(self):
        self.logUpdate('[INFO]get {}mer probes from positive genome'.format(self.pLen1))
        f = open(self.window1)
        w = open(self.posKmers, 'w')
        for title, seq in SimpleFastaParser(f):
            header = title.strip()
            genomeLen = len(seq.strip())
            lines = [f'>{header}_{i}\n{seq[i:i + self.pLen1]}\n' for i in range(genomeLen - self.pLen1+1)]
            w.writelines(lines)
        f.close()
        w.close()

    def remNegative(self):
        negKmerResult = ProbeitUtils.defineFile(self.maskingDir, 'mmseqs.search')
        negKmerCoords = ProbeitUtils.defineFile(self.thermoFilteringDir, 'mmseqs.txt')
        deDupGenomeCoords = ProbeitUtils.defineFile(self.dedupDir, 'dedup.bed') 
        # REMOVE NEGATIVE K-MERS
        self.logUpdate('[INFO]remove probes found in negative genome')
        ProbeitUtils.ridNegKmers(self.posKmers, self.negGenome, negKmerResult, self.maskingDir, self.ridNegId, self.threads)
        # MAKE DEDUPLICATED POSITIVE GENOME COORDINATE FILE(BED)
        with open(self.window1) as f:
            with open(deDupGenomeCoords, 'w') as w:
                for title, seq in SimpleFastaParser(f):
                    header = title.strip()
                    w.write(header + '\t0\t' + str(len(seq.strip())) + '\n')
        if self.negGenome: 
            with open(negKmerResult) as f:
                with open(negKmerCoords, 'w') as w:
                    for line in f:
                        header = line.split()[0].strip()
                        startPos = int(header.split('_')[-1])
                        header = '_'.join(header.split('_')[:-1])
                        w.write(header + '\t' + str(int(startPos)) + '\t' + str(int(startPos) + self.pLen1) + '\n')
                        w.write(f'{header}\t{startPos}\t{startPos+self.pLen1}\n')
            self.logUpdate(ProbeitUtils.getSubtractedBed(deDupGenomeCoords, negKmerCoords, self.negRemCoords))

    # TODO change name
    def getThermoImproperKmers(self):
        self.logUpdate('[INFO]filter probes with thermodynamic features') 
        with open(self.kmersForTheromFilter, 'w') as w:
            for title, seq in SimpleFastaParser(open(self.posKmers)):
                header = ('>' + title).split('_')
                w.write('_'.join(header[:-1]) + '\t' + header[-1] + '\n')
                w.write(seq + '\n')
        if self.needThermoFilter:
            self.thermoFilter()
        
    def makeWindow2(self):
        def _getScKmer(scKmer, window):
            temp = scKmer.split(',')
            g = int(temp[0])
            s = int(temp[1])
            e = s + self.pLen1
            window[g] = window[g][0], window[g][1][:s]+'N'*40+window[g][1][e:]
            return g, s, e
        def _getWinSeq(kmer, lenList):
            g = kmer[0] 
            s = kmer[1]
            e = kmer[2]
            gLen = lenList[g]
            ws = s - self.windowSize if s>self.windowSize else 0
            we = e + self.windowSize 
            we = gLen if we>gLen else we
            genome, seq = window1List[g]
            seq = seq[ws:we] 
            return f'>{genome}:{s}:{e}\t{ws}\t{we}\t{gLen}\n{seq}\n'
        window1List =[(h, s) for h, s in SimpleFastaParser(open(self.window1))]
        genomeLens = [len(i[1]) for i in window1List]
        scResultList = [_getScKmer(kmer, window1List) for l in open(self.scResult1) for kmer in ProbeitUtils.parseKmers(l)[1:]]
        window2Lines = [_getWinSeq(kmer, genomeLens) for kmer in scResultList]
        with open(self.window2, 'w') as w:
            w.writelines(window2Lines)
        ProbeitUtils.makeLookup(self.window2, self.lookup2)

    def trimProbes(self):
        # MAKE PROBEs UNERSTANDABLE
        with open(self.lookup1)as f:
            genomeKeys = {int(i.split()[0]):i.split()[1] for i in f}
        cnt = 0
        namesForProbes1 = {}
        lines1 = []
        for h, s in SimpleFastaParser(open(self.tempProbe1)):
            keys, pos = ProbeitUtils.parseGenmapPattern(h)
            keys = [genomeKeys[k] for k in keys]
            probes = ['{}:{}:{}'.format(keys[i], pos[i], pos[i] + self.pLen1) for i in range(len(keys))]
            for i in probes:
                namesForProbes1[i] = 'p1_{}'.format(cnt)
            lines1.append('>p1_{}\t{}\n{}\n'.format(cnt, ';'.join(probes), s))
            cnt += 1
        with open(self.probe1, 'w') as w:
            w.writelines(sorted(lines1))
        if not self.needProbe2:
            return
        with open(self.lookup2)as f:
            probe1Keys = {int(i.split()[0]):i.split()[1] for i in f}
        cnt = 0
        lines2 = []
        for h, s in SimpleFastaParser(open(self.tempProbe2)):
            keys, pos = ProbeitUtils.parseGenmapPattern(h)
            keys = set([namesForProbes1[probe1Keys[k]] for k in keys])
            lines2.append('>p2_{}\t{}\n{}\n'.format(cnt, ';'.join([p for p in keys]), s) )
            cnt += 1
        with open(self.probe2, 'w') as w:
            w.writelines(sorted(lines2))

    def execute(self):
        # MAKE 1st probe
        self.logUpdate("[INFO] make 1st probes")
        # self.filterInputGenome()
        self.remRedundancy()
        self.makeWindow1()
        self.extractKmersFromInputGenome()
        if self.negGenome:
            self.remNegative()
        self.getThermoImproperKmers()
        posCoords = self.negRemCoords if self.negGenome else '' 
        improperKmers = self.getImproperKmers()
        msg = ProbeitUtils.makeProbe(
            self.tempProbe1, self.window1, self.lookup1, self.idxDir1, self.cmDir1, self.scDir1, self.pLen1, 
            self.cmError1, self.scCoverage1, self.scEarlyStop1, self.scScore1, self.scRepeats1,
            selCoords=posCoords, thread=self.threads, improperKmers=improperKmers
        )
        self.logUpdate(msg, False)
        if not self.needProbe2:
            self.trimProbes()
            return
        self.logUpdate("[INFO] make 2nd probes")
        # MAKE 2nd WINDOW
        self.makeWindow2()
        # MAKE 2nd PROBES
        msg = ProbeitUtils.makeProbe(
            self.tempProbe2, self.window2, self.lookup2, self.idxDir2, self.cmDir2, self.scDir2, self.pLen2,
            self.cmError2, self.scCoverage2, self.scEarlyStop2, self.scScore2, self.scRepeats2, 
            thread=self.threads, overlap=True
        )
        self.logUpdate(msg, False)
        self.trimProbes()
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
        print(" --not-cluster NONE")
        print("\t Use it when you DO NOT need to cluster positive genome")
        print(" --not-make-probe2 NONE")
        print("\t Use it when you DO NOT need to make 2nd probes")
        print(" --not-thermo-filter NONE")
        print("\t Use it when you DO NOT need the thermodynamic filter")
        print(" --probe1-len INT[40]")
        print("\t Length of 1st Probes")
        print(" --probe2-len INT[20]")
        print("\t Length of 2nd Probes")

        print("OPTIONS FOR GENMAP PART: Genmap calculates mappability by summarizing k-mers.")
        print(" --probe1-error INT[0]")
        print("\t The number of error allowed in 1st Probes")
        print(" --probe2-error INT[1]")
        print("\t The number of error allowed in 2nd Probes")

        print("OPTIONS FOR SETCOVER PART: Setcover makes probesets cover targets with the minimum number of probes.")
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
                validPosList.append(p)
            else:
                print('[ERROR] {} is not proper for position list.')
        self.posList = validPosList
        if not self.isMaxWindow:
            self.windowSize = self.windowSize - (max(self.posList) - min(self.posList))
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
        # DIRECTORIES
        self.workDir = ProbeitUtils.defineDirectory(self.workDir.split(os.path.sep)[0], make=True)
        self.tempDir = ProbeitUtils.defineDirectory('temp', make=True, root=self.workDir)
        self.searchDir = ProbeitUtils.defineDirectory('blast', make=True, root=self.workDir)
        self.probe1byPosDir = ProbeitUtils.defineDirectory('probe1byPos', make=True, root=self.workDir)
        if self.needProbe2:
            self.inputDir2 = ProbeitUtils.defineDirectory('input2', make=True, root=self.workDir)
            self.idxDir2 = ProbeitUtils.defineDirectory('index', make=False, root=self.workDir)
            self.cmDir2 = ProbeitUtils.defineDirectory('mapping2', make=True, root=self.workDir)
            self.scDir2 = ProbeitUtils.defineDirectory('setcover2', make=True, root=self.workDir)
        # FILES
        self.log = ProbeitUtils.defineFile(self.workDir, 'log.txt')
        self.posProbeCSVs = {pos: ProbeitUtils.defineFile(self.probe1byPosDir, f'pos{pos}.csv') for pos in self.posList}
        self.posProbeCSVs[-1] = ProbeitUtils.defineFile(self.probe1byPosDir, 'merged.csv')
        self.probe1 = ProbeitUtils.defineFile(self.workDir, 'probe1.fa')
        if self.needProbe2:
            self.window1 = ProbeitUtils.defineFile(self.inputDir2, 'masking.fasta')
            self.window1Coords = ProbeitUtils.defineFile(self.inputDir2, 'masked.fasta')
            self.maskedGenome= ProbeitUtils.defineFile(self.inputDir2, 'masked.fa')
            self.window2Coords = ProbeitUtils.defineFile(self.inputDir2, 'window.bed')
            self.lookup2 = ProbeitUtils.defineFile(self.workDir, 'name.lookup')
            self.window2 = ProbeitUtils.defineFile(self.inputDir2, 'seq.fa')
            self.lookupTSV = ProbeitUtils.defineFile(self.workDir, 'lookup.tsv')
            self.tempProbe2 = ProbeitUtils.defineFile(self.workDir, 'temp2.fa')
            self.probe2 = ProbeitUtils.defineFile(self.workDir, 'probe2.fa')
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
            ProbeitUtils.searchSNPs(self.searchDir, searchProbe, self.strGenome, snpKmers, self.searchKmer, thread=self.threads)
        )
        return snpKmers

    @staticmethod
    def _parseSearcgResult(searchResult):
        print(searchResult)
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
        minPos = min(self.posList)
        maxPos = max(self.posList)
        refSeq = self.getReferenceSeq(self.refGenome)
        if not refSeq:
            self.logUpdate('[warn]Failure to get reference sequence from reference genome.')
            self.printUsage()
        for snp in self.snpList:
            self.logUpdate('[INFO]SNP {}'.format(snp))
            mutType, orf, mutation = self.parseMutation(snp)
            if mutType == 'aa':
                if self.refGenomeAnnot == '':
                    self.logUpdate('[warn]For Amino Acid based SNPs reference annotation needed.')
                    continue
                orfStartPos = self.getOrfStartPos(self.refGenomeAnnot, orf)
                if orfStartPos == -1:
                    self.logUpdate('[warn]Failure to find snp {} in reference annotaion.'.format(snp))
                    continue
                aa1, aa2, mutPos = mutation[0], mutation[-1], int(mutation[1:-1])
                codonStartPos = orfStartPos + (mutPos - 1) * 3 - 1
                codonEndPos = orfStartPos + mutPos * 3 - 1
                
                refCodon = refSeq[codonStartPos: codonEndPos]
                if aa1 != Seq(refCodon).translate():
                    self.logUpdate('[warm]Failure to find SNP {} in reference genome'.format(snp))
                    continue
                
                seqWithSNP = refSeq[codonStartPos - (maxPos-1): codonEndPos + (self.probLen1 - minPos)]
                searchResult = self._searchSnpFromStrGenome(mutation, seqWithSNP)
                trimmedResult = self._trimSearchResult(searchResult, seqWithSNP, mutType, maxPos, snp, codonStartPos, refCodon)
                if trimmedResult.empty:
                    self.logUpdate(f'[WARN] Problems occured searching snp {snp} in the strain genome or reference genome')
                    continue
                wtSequence, stSequence, ntSNP, locSnp, found = self._parseSearcgResult(searchResult=trimmedResult)
                self.logUpdate('[INFO]aa:{}:{} converted to nt:{}'.format(orf, mutation, ntSNP))
                aaOrf = '{}:{}'.format(orf, mutation)
                mutSeqs = ParaSeqs(ntSNP, aaOrf, wtSequence, stSequence, mutLoc=locSnp, probLen=self.probLen1)
            elif mutType == 'nt':
                nt1, nt2, snpPos = mutation[0], mutation[-1], int(mutation[1:-1])
                refNT = refSeq[snpPos]
                if nt1 != refNT:
                    self.logUpdate('[warn]Failure to find SNP {} in reference genome'.format(snp))
                    continue
                seqWithSNP = refSeq[snpPos - (maxPos - 1):snpPos + 1 + (self.probLen1 - minPos)]
                searchResult = self._searchSnpFromStrGenome(mutation, seqWithSNP)
                trimmedResult = self._trimSearchResult(searchResult, seqWithSNP, mutType, maxPos, snp)
                if trimmedResult.empty:
                    self.logUpdate(f'[WARN] Problems occured searching snp {snp} in the strain genome or reference genome')
                    continue
                wtSequence, stSequence, ntSNP, locSnp, found = self._parseSearcgResult(searchResult=trimmedResult)
                mutSeqs = ParaSeqs(ntSNP, '', wtSequence, stSequence, mutLoc=locSnp, probLen=self.probLen1)
            print('DONE')
            if found < 0 or not found:
                self.logUpdate('[WARN]Failure to find SNP {} in strain genome'.format(snp))
                continue
            self.probesByPos[-1].append(mutSeqs)
            for pos in self.posList:
                wtProbe, stProbe = mutSeqs.getProbesWithPos(pos)
                paraSeq = ParaSeqs(mutSeqs.ntSnp, mutSeqs.aaSnp, wtProbe, stProbe, found=found, probLen=self.probLen1)
                self.probesByPos[pos].append(paraSeq)
            else:
                self.logUpdate('[warn]SNP {} has a wrong format.'.format(snp))
                continue
                

    def makeProbe1(self):
        probeLines = []
        for pos in self.posList + [-1]:
            probCSV = self.posProbeCSVs[pos]
            csvWriter = open(probCSV, 'w')
            csvWriter.write('WT sequence,ST sequence,found,ntSNP,aaSNP\n')
            for p in self.probesByPos[pos]:
                csvWriter.write(f'{p.wtSeq},{p.stSeq},{p.found},{p.ntSnp},{p.aaSnp}\n')
                if pos != -1:
                    probeLines.append(f'>{p.ntSnp}{"=" + p.aaSnp if p.aaSnp else ""};{pos}\n{p.stSeq}\n')
            csvWriter.close()
        with open(self.probe1, 'w') as w:
            w.writelines(sorted(probeLines))

    @staticmethod
    def makeMaskedCoords(inputDF, outputBed):
        inputDF['seqID'] = inputDF['seqID'].apply(lambda x: x.split(';')[0])
        inputDF.to_csv(outputBed, sep='\t', header=False, index=False)

    def makeWindows(self):
        print( self.probesByPos[-1])
        window1 = [f'>{p.ntSnp}{ "=" + p.aaSnp if p.aaSnp else ""}\n{p.stSeq}\n' for p in self.probesByPos[-1]]
        if not window1:
            self.logUpdate('[ERROR] Cannot find any SNP in strain genomes')
            self.printUsage()
        with open(self.window1, 'w') as w:
            w.writelines(window1)
        msg = ProbeitUtils.getPatternPosition(self.window1, self.strGenome, self.lookupTSV)
        self.logUpdate(msg, False)
        self.makeMaskedCoords(pd.read_csv(self.lookupTSV, sep='\t')[['seqID', 'start', 'end']], self.window1Coords)
        ProbeitUtils.getWindowFasta(
            self.strGenome, self.window1Coords, self.maskedGenome, self.window2Coords, self.window2, self.windowSize
        )
        ProbeitUtils.makeLookup(self.window2, self.lookup2)

    def trimProbes(self):
        maskDF = pd.read_csv(self.lookupTSV, sep='\t')
        kmers = list(maskDF['patternName'])
        probe2lines = []
        for h, s in SimpleFastaParser(open(self.tempProbe2)):
            kmerIndex = ProbeitUtils.parseGenmapPattern(h)[0]
            coveredSNPs = list((set([kmers[i] for i in kmerIndex])))
            probe2lines.append(f'>{":".join(coveredSNPs)}\n{s}\n')
        with open(self.probe2, 'w') as w:
            w.writelines(sorted(probe2lines))

    def execute(self):
        # MAKE 1st PROBEs
        self.logUpdate("[INFO]make 1st probes")
        self.findSNPs()
        self.makeProbe1()
        if not self.needProbe2:
            return
        # MAKE 2ND PROBEs
        self.logUpdate("[INFO]make 2nd probes")
        self.makeWindows()
        msg = ProbeitUtils.makeProbe(
            self.tempProbe2, self.window2, self.lookup2, self.idxDir2, self.cmDir2, self.scDir2, self.probeLen2, 
            self.cmError2, self.scCoverage2, self.scEarlyStop2, self.scScore2, self.scRepeats2, 
            thread=self.threads, overlap=True
        )
        self.logUpdate(msg, False)
        self.trimProbes()
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

