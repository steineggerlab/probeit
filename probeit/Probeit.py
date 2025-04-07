#!/usr/bin/env python
from .config import Config
from Bio.SeqIO.FastaIO import SimpleFastaParser
from Bio.Seq import Seq, reverse_complement
import pandas as pd
import primer3
import numpy as np
import os
import shutil
import getopt
import subprocess
import re
import paraSeq
class ProbeitUtils:    
    # RUN COMMANDLINE
    @staticmethod
    def runCommand(command, verbose=False):
        print('[CLI] '+command)
        if verbose:
            commandList = command.split()
            sp = subprocess.Popen(commandList, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            stdout, stderr = sp.communicate()
            return stdout.decode('UTF-8'), stderr.decode('UTF-8')

        os.system(command)

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
        directory = f'{root}{dirName}{os.path.sep}'
        if make:
            os.makedirs(directory)
        return directory

    @staticmethod
    def getFileName(file, keepFilenameExtension=True):
        file = file.split(os.path.sep)[-1]
        if keepFilenameExtension:
            return file 

        return '.'.join(file.split('.')[:-1])

    @classmethod
    def delDir(cls, directory):
        if os.path.isdir(directory):
            print(f'[INFO] The directory named {directory} is removed not.')
            shutil.rmtree(directory)
            return

        print(f'[INFO] The directory named {directory} does not exist.')

    # FASTA FILES RELATED
    @classmethod
    def sortFasta(cls, inputFasta):
        fastaList = [(h,s) for h, s in SimpleFastaParser(open(inputFasta))]
        fastaList = sorted(['>{}\n{}\n'.format(i[0], i[1]) for i in fastaList])
        with open(inputFasta, 'w') as w:
            w.writelines(fastaList)

    @classmethod
    def getSubseqFasta(cls, coordinateBed, inputFasta, outputFasta):
        command1 = "seqkit subseq --quiet --bed {} {} > {}".format(coordinateBed, inputFasta, outputFasta)
        cls.runCommand(command1)
        return '[CLI] '+command1+'\n'

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

    @staticmethod
    def simplifyFastaHeaders(inputFasta, outputFasta):
        f = open(inputFasta)
        w = open(outputFasta, 'w')
        for title, seq in SimpleFastaParser(f):
            w.write(('>' + title).split()[0].strip() + '\n')
            w.write(seq + '\n')

        f.close()
        w.close()

    @staticmethod
    def extractKmers(genomeFasta, kmersFasta, pLen):
            kmers = {}
            genomeIdx = 0
            w = open(kmersFasta, 'w')
            for _, s in SimpleFastaParser(open(genomeFasta)):
                l = len(s)
                for pos in range(l-pLen+1):
                    name = f'{genomeIdx},{pos}'
                    seq = s[pos:pos+pLen].upper()
                    kmers.setdefault(seq,[])
                    kmers[seq].append(name)
                genomeIdx += 1
            for s, names in kmers.items():
                rep = names[0]
                w.write(f'>{rep};{"|".join(names)}\n{s}\n')
            w.close()

    @classmethod
    def getPatternPosition(cls, patternFasta, genomeFasta, positionsTSV):
        command = "seqkit locate -f {} {} > {}".format(patternFasta, genomeFasta, positionsTSV)
        cls.runCommand(command)
        df = pd.read_csv(positionsTSV, sep='\t')
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
        command = (
            f'mmseqs easy-search {posKmers} {negative} {output} {tempDir} -v 3 --spaced-kmer-mode 0 -k 13 --mask 0 ' +
            f'--format-output query,target,pident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue ' +
            f'-c 0.9 --min-seq-id {seqInProbe} --cov-mode 2 --alignment-mode 4 --search-type 3 --threads {thread}' +
            ' --remove-tmp-files 0'
        )
        stdOut, stdErr = cls.runCommand(command, verbose=True)
        print(stdOut)
        print(stdErr)
        return f'[CLI] {command}\n' + stdOut + stdErr

    # COMPUTE MAPPABILITY
    @classmethod
    def simpleComputeMappability(cls, genome, lookup, positionBED, outputCSV, pLen, improperKmers):
        positiveKmers = set()
        genomeNames = [i.split()[1].strip() for i in open(lookup)]
        w = open(outputCSV, 'w')
        for line in open(positionBED):
            g, kmerSeq, e = line.strip().split()
            kmerSeq = int(kmerSeq)
            e = int(e)
            g = genomeNames.index(g)
            for i in range(kmerSeq,e):
                positiveKmers.add(f'{g},{i}')
        kmerAndNames = {}
        idx = 0
        nucleotides = {'A', 'C', 'G', 'T'}
        for h, kmerSeq in SimpleFastaParser(open(genome)):
            l = len(kmerSeq)
            for i in range(l - pLen + 1):
                kmerName = f'{idx},{i}'
                seq = kmerSeq[i:i+pLen].upper()
                if set(seq) != nucleotides & set(seq):
                    continue

                kmerAndNames.setdefault(seq,[])
                kmerAndNames[seq].append(kmerName)
            idx += 1

        for kmerSeq, kmerNames in kmerAndNames.items():
            repName = kmerNames[0]
            isPositive = {repName} & positiveKmers != set()
            isImproper = kmerSeq in improperKmers
            if isPositive and not isImproper:
                w.write(f'{repName};{"|".join(kmerNames)}\n')
        w.close()

    @classmethod
    def computeMappability(cls, genome, indexDir, error, length, outputDir, outputCSV, threads, improperKmers):
        w = open(outputCSV, 'w')
        cls.delDir(indexDir)
        command0 = " > /dev/null"
        command1 = "genmap index -F {} -I {}".format(genome, indexDir)
        cls.runCommand(command1 + command0)
        command2 = "genmap map --no-reverse-complement -E {} --csv -K {} -t -b --frequency-large -I {} -O {} -T {}"
        command2 = command2.format(error, length, indexDir, outputDir, threads)
        cls.runCommand(command2 + command0)
        inputCSV = cls.defineFile(outputDir, f'{cls.getFileName(genome, False)}.genmap.csv')

        for line in open(inputCSV):
            [repName, firstName] = line.split('|')[0].strip().split(';')
            if repName == firstName and not {repName} & improperKmers:
                w.write(line)
        w.close()
        return '[CLI] {}{}\n[CLI] {}{}\n'.format(command1, command0, command2, command0)

    # BED FILES RELATED
    @classmethod
    def getSubtractedBed(cls, positiveBed, negativeBed, bedFileName):
        command = "bedtools subtract -a {} -b {} > {}".format(positiveBed, negativeBed, bedFileName)
        cls.runCommand(command)
        return '[CLI] '+command

    # LOOKUP FILES RELATED
    @classmethod
    def makeLookup(cls, windowFasta, lookup, genomePos=''):
        w = open(lookup, 'w')
        headers = [header.strip() for header, _ in SimpleFastaParser(open(windowFasta))]
        for i, header in enumerate(headers):
            w.write(f'{i}\t{header}\n')
        w.close()
        if not genomePos:
            return

        w = open(genomePos, 'w')
        lengths = [len(seq) for _, seq in SimpleFastaParser(open(windowFasta))]
        for i, length in enumerate(lengths):
            w.write(f'{i}\t0\t{length}\n')
        w.close()
        
    # SNPs RELATED
    @classmethod
    def searchSNPs(cls, workDir, inputFasta, strGenomeFasta, result, kmer, threads):
        searchDir = cls.defineDirectory('search', root=workDir, make=True)
        tempDir = cls.defineDirectory('temp', root=searchDir, make=False)
        searchDB = cls.defineFile(searchDir, 'searchDB' )
        strainDB =  cls.defineFile(searchDir, 'strainDB')
        aln = cls.defineFile(searchDir , 'mmseqs.aln')
        cmd0 = ' --threads {}'.format(threads)
        cmd1 = 'mmseqs createdb {} {}'
        cmd2 = 'mmseqs createdb {} {}'
        cmd3 = 'mmseqs search {} {} {} {} --search-type 3 -k {}'
        cmd4 = 'mmseqs convertalis {} {} {} {} --format-output target,query,tseq,tstart,tend --search-type 3'
        out1, err1 = cls.runCommand(cmd1.format(inputFasta, searchDB), verbose=True)
        out2, err2 = cls.runCommand(cmd2.format(strGenomeFasta, strainDB), verbose=True)
        out3, err3 = cls.runCommand(cmd3.format(searchDB, strainDB, aln, tempDir, kmer, threads) + cmd0, verbose=True)
        out4, err4 = cls.runCommand(cmd4.format(searchDB, strainDB, aln, result, threads) + cmd0, verbose=True)
        df = pd.read_csv(result, sep='\t', header=None)
        df.columns = ['substr', 'snp', 'strseq', 'start', 'end']
        df['aln'] = df.apply(lambda x: x[2][int(x[3]-1):int(x[4])], axis=1)
        df['len'] = df.aln.apply(lambda x: len(x)-1)
        df = df[['substr', 'snp', 'len', 'aln']]
        df.to_csv(result, header=False, index=False, sep='\t')
        print(err1 + err2 + err3 + err4)
        ProbeitUtils.delDir(searchDir)
        return out1 + out2 + out3 + out4

    # PARSE GENMAP KMER PATTERNS 
    @staticmethod
    def parseKmers(line):
        return re.findall(r'[0-9]+,[0-9]+', line)

    @staticmethod
    def parseGenmapPattern(header):
        p1 = re.compile('[0-9]+,')
        p2 = re.compile(',[0-9]+')
        return [int(i[:-1]) for i in p1.findall(header)], [int(i[1:]) for i in p2.findall(header)]

    # ARGUMENTS PARSING
    @staticmethod
    def getArgList(value, isInt=False):
        if isInt:
            return [int(i.strip()) for i in value.split(',')]

        return value.split(',')

    # OTHERS
    @classmethod
    def setCover(cls, coverage, length, eStop, dist, reps, mapCSV, genome, lookup, setcoverResultBed, probeLen):
        w = open(setcoverResultBed, 'w')
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
        for line in stdOut.strip().split('\n'):
            matchedKmers = line.split(';')[1].strip()
            idx = int(line.split(';')[0].split(',')[0])
            pos = line.split(';')[0].split(',')[1]
            w.write('\t'.join([genomeAndIdx[idx], pos, str(int(pos) + probeLen), matchedKmers]) + '\n')
        w.close()
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
        self.query = ['p1_intrinsic_probs', f'p1_tm<{self.minProbeTm}', f'p1_GC_perc<{self.minGC}', f'p1_GC_perc>{self.maxGC}', f'p1_homodimer_tm>{self.maxHomoDimerTm}', f'p1_hairpin_tm>{self.maxHairpinTm}']
        self.sorter = ['p1_hairpin_tm', 'p1_homodimer_tm']
        self.ascending = [True, True, True, True] if self.isLigational else [True, True]
        if self.isLigational:   
            self.query += ['p2_intrinsic_probs', f'p2_tm<{self.minProbeTm}', f'p2_GC_perc<{self.minGC}', f'p2_GC_perc>{self.maxGC}', f'p2_homodimer_tm>{self.maxHomoDimerTm}', f'p2_hairpin_tm>{self.maxHairpinTm}']
            self.sorter += ['p2_hairpin_tm', 'p2_homodimer_tm']

    def makeLogMessage(self, numKmers, numImpKmers):
        msg = '[INFO]filter probes with thermodynamic features\n' 
        msg += f"\tMinimum Tm: {self.minProbeTm}\n"
        msg += f"\tMinimum GC percentage: {self.minGC }\n"
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
                p1.append(rcKmer[0 : mid_pos])
                p2.append(rcKmer[mid_pos : self.pLen])

        self.kmerList = p1 + p2  
        cols = ['id', 'chromStart', 'chromEnd', 'genome_segment', 'rc', 'p1']      
        outputDF = pd.DataFrame(list(zip(identity, posStartList, posEndList, kmers, rc, p1)), columns=cols)
        if self.isLigational:
            outputDF['p2'] = p2
        return outputDF

    def getThermoFeaturesDF(self):
        thermoFeatures = pd.DataFrame(list(set(self.kmerList)), columns=['p'])
        thermoFeatures = thermoFeatures.assign(ultimate_base=thermoFeatures['p'].str[-1])
        thermoFeatures = thermoFeatures.assign(penultimate_base=thermoFeatures['p'].str[-2])
        thermoFeatures = thermoFeatures.assign(tm=np.vectorize(primer3.calcTm)(thermoFeatures['p']))
        thermoFeatures['hairpin_tm'] = list(map(primer3.calcHairpinTm, thermoFeatures['p']))
        thermoFeatures = thermoFeatures.assign(homodimer_tm=np.vectorize(primer3.calcHomodimerTm)(thermoFeatures['p']))
        thermoFeatures = thermoFeatures.assign(intrinsic_probs=np.vectorize(self.hasLowComplexity)(thermoFeatures['p']))
        thermoFeatures = thermoFeatures.assign(GC_perc=np.vectorize(self.getContentGC)(thermoFeatures['p']))
        return thermoFeatures

    def getThermoFeaturedKmersDF(self, kmersDF, thermoFeaturesDF):
        joinedDF = pd.merge(kmersDF, thermoFeaturesDF.add_prefix('p1_'), how='left', left_on='p1', right_on='p1_p')
        if self.isLigational:   
            joinedDF = pd.merge(joinedDF, thermoFeaturesDF.add_prefix('p2_'), how='left', left_on='p2', right_on='p2_p')
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


class Probeit:
    args = []
    subWork = None

    def __init__(self, args):
        self.args = args

    def do(self):
        if self.args == []:
            self.printUsage()

        workflow = self.args[0]
        match workflow:
            case '-h' | '--help':
                self.printUsage()

            case 'posnegset':
                print('CURRENT: ', os.getcwd())
                subWork = PosNegSet(self.args[1:])
                subWork.run()
                return

            case 'snp':
                print('CURRENT: ', os.getcwd())
                subWork = SNP(self.args[1:])
                subWork.run()
                return

            case 'primer':
                print('CURRENT: ', os.getcwd())
                subWork = Primer(self.args[1:])
                subWork.run()
                return

            case _:
                self.printUsage()

    @staticmethod
    def printUsage():
        print("PROBEIT V2.2")
        print('probeit <workflow> [<args>]')
        print("WORKFLOWS")
        print("posnegset: make two-sets probes with positive and negative sets")
        print("snp: make two-sets probes with wild type genome, strain genome and SNPs")
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
        'not-make-probe2', 'remove-redundancy', 'not-thermo-filter', 'ligation-probe', 'probe1-len=', 'probe2-len=', 'window-size=', 'threads=',
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
    needProbe2 = True
    doThermoFilter1 = True # thermo-filter
    doThermoFilter2 = True
    isLigationProbe = False # thermo-filter
    windowSize = 200
    scEarlyStop1 = 0.9
    scScore1 = 11
    scRepeats1 = 1
    scEarlyStop2 = 0.99
    scScore2 = 20
    scRepeats2 = 10
    threads = 8
    impKmers1 = set()
    impKmers2 = set()
    probe1Index = []
    lenSeq = {}

    def __init__(self, args):
        args = getopt.getopt(args, self.shortParams, self.longParams)[0]
        print(f"Your arguments: probeit posnegset {ProbeitUtils.getUserArgs(args)}")
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
                self.doRemoveRedundancy = True if opt == '--remove-redundancy' else self.doRemoveRedundancy
                self.needProbe2 = False if opt == '--not-make-probe2' else self.needProbe2
                self.doThermoFilter1 = False if opt == '--not-thermo-filter' else self.doThermoFilter1
                self.doThermoFilter2 = False if opt == '--not-thermo-filter' else self.doThermoFilter2
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
                # identity in mmeqs search
                self.ridNegId = float(val) if opt == '--rid-neg-id' else self.ridNegId
                # setcover similarity
                self.scScore1 = int(val) if opt == '--probe1-dist' else self.scScore1
                self.scScore2 = int(val) if opt == '--probe2-dist' else self.scScore2
            except Exception as e:
                self.printUsage()

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
        self.genomeFASTA  = self.inputGenome 
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
            self.window2FASTA = ProbeitUtils.defineFile(self.inputDir2, Config.window)
            self.window2PosBED = ProbeitUtils.defineFile(self.inputDir1, 'window.bed') 
            self.lookup2 = ProbeitUtils.defineFile(self.workDir, 'probe1.lookup')
            self.kmers2FASTA = ProbeitUtils.defineFile(self.inputDir2, 'kmers.fa')
            self.thermoCalc2TSV = ProbeitUtils.defineFile(self.thermoFiltering2Dir, 'thermo_uncut.tsv')
            self.tempProbe2 = ProbeitUtils.defineFile(self.workDir, 'temp2.fa')
            self.probe2 = ProbeitUtils.defineFile(self.workDir, Config.probe2)
            self.rcProbe2 = ProbeitUtils.defineFile(self.workDir, Config.rcprobe2)
        self.logUpdate('[INFO]Your arguments: probeit posnegset ' + ' '.join(['{} {}'.format(i[0], i[1]).strip() for i in args]), False)
        
    def logUpdate(self, msg, doPrint=True):
        if doPrint:
            print(msg)

        with open(self.log, 'a') as w:
            w.write(msg + '\n')

    def negFilter(self):
        w = open(self.window1PosBED, 'w')
        negKmerResult = ProbeitUtils.defineFile(self.maskingDir, 'search.tsv')
        negKmerPosBED = ProbeitUtils.defineFile(self.maskingDir, 'search.bed') 
        # REMOVE NEGATIVE K-MERS
        self.logUpdate('[INFO]remove probes found in negative genome')
        self.logUpdate(ProbeitUtils.ridNegKmers(self.posKmers1FASTA, self.negGenome, negKmerResult, self.maskingDir, self.ridNegId, self.threads))
        
        # MAKE DEDUPLICATED POSITIVE GENOME COORDINATE FILE(BED)
        idx = 0
        for _, seq in SimpleFastaParser(open(self.window1FASTA)):
            w.write(f'{idx}\t0\t{str(len(seq.strip()))}\n')
            idx += 1
        
        # EXTRACT NEGATIVE REMOVED POSITIONS BED
        w = open(negKmerPosBED, 'w')
        for i in open(negKmerResult):
            kmers = ProbeitUtils.parseKmers(i)[1:]
            for kmer in kmers:
                g, s = [int(i) for i in kmer.split(',')]
                w.write(f'{g}\t{s}\t{s+40}\n')       
        w.close()
        self.logUpdate(ProbeitUtils.getSubtractedBed(self.window1PosBED, negKmerPosBED, self.negRemPosBED), False)

    def getNegativeKmers(self): 
        negativeKmers = set()
        if not self.negGenome:
            return negativeKmers

        prevGenome = ""
        prevEnd = 0
        for line in open(self.negRemPosBED):
            g, s, e = line.strip().split()
            s = int(s)
            e = int(e)
            if g == prevGenome:
                for i in range(prevEnd, s):
                    negativeKmers.add(f'{g},{i}')
                prevEnd = e
                continue

            for i in range(prevEnd, self.lenSeq[prevGenome]):
                negativeKmers.add(f'{prevGenome},{i}')
            prevGenome = g
            prevEnd = e
            for i in range(0,s):
                negativeKmers.add(f'{g},{i}')
        else:
            for i in range(prevEnd, self.lenSeq[prevGenome]):
                negativeKmers.add(f'{prevGenome},{i}')
        return negativeKmers
        
    def easyComMap(self, uniqComMap):
        w = open(uniqComMap, 'w')
        alphabets = {'A', 'C', 'G', 'T'}
        seqsWithLength = dict()
        idx = 0
        for h, s in SimpleFastaParser(open(self.window1FASTA)) :
            seqsWithLength[str(idx)] = len(s)
            idx += 1     
        # negativeKmers = set()
        # prevGenome = '0'
        # prevEnd = 0
        negativeKmers = self.getNegativeKmers()
        for h, s in SimpleFastaParser(open(self.posKmers1FASTA)):
            kmers = set(ProbeitUtils.parseKmers(h))
            isThermoImproper = kmers&self.impKmers1
            notOnlyATGC = set(s) != alphabets&set(s)
            # TODO
            # isNegative = kmers&negativeKmers
            isNegative = {ProbeitUtils.parseKmers(h)[0]} & negativeKmers
            if isNegative or isThermoImproper or notOnlyATGC:
                continue
            w.write(h+'\n')
        w.close()
    
    def makeWindow1(self):
        ProbeitUtils.sortFasta(self.genomeFASTA)
        ProbeitUtils.simplifyFastaHeaders(self.genomeFASTA, self.window1FASTA)
        ProbeitUtils.makeLookup(self.window1FASTA, self.lookup1, self.window1PosBED)

    def makeWindow2(self):
        def _getScKmer(scKmer, window):
            temp = scKmer.split(',')
            g = int(temp[0])
            s = int(temp[1])
            e = s + self.pLen1
            window[g] = window[g][0], window[g][1][:s]+ 'N'*self.pLen1 + window[g][1][e:]
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

        w = open(self.window2FASTA, 'w')
        window1Seqs =[(h, s) for h, s in SimpleFastaParser(open(self.window1FASTA))]
        window1Lens = [len(i[1]) for i in window1Seqs]
        probe1List = [[_getScKmer(kmer, window1Seqs) for kmer in ProbeitUtils.parseKmers(line)] for line in open(self.scPosBed1)]
        for i, kmers in enumerate(probe1List):
            for j, kmer in enumerate(kmers):
                w.write(_getWinSeq(kmer, window1Lens, f'probe_{i}_{j}'))
        w.close()
        ProbeitUtils.makeLookup(self.window2FASTA, self.lookup2, self.window2PosBED)


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

        probe1Writer = open(self.probe1, 'w')
        rcprobe1Writer = open(self.rcProbe1, 'w')
        probeIdx = 0
        for seq, kmers in {s:ProbeitUtils.parseKmers(h) for h, s in SimpleFastaParser(open(self.tempProbe1))}.items():
            rcSeq = Seq(seq).reverse_complement()
            parsedKmers = []
            for kmerIdx, kmer in enumerate(kmers):
                [genome,pos] = kmer.split(',')
                genome = genomeKeys[int(genome)]
                parsedKmers.append(f'probe_{probeIdx}_{kmerIdx}:{genome}:{pos}')
                kmerIdx += 1
            self.probe1Index += [kmer.split(':')[0] for kmer in parsedKmers]
            probe1Writer.write(  f">probe_{probeIdx}\t{'|'.join(parsedKmers)}\n{seq}\n"  )
            rcprobe1Writer.write(f">probe_{probeIdx}\t{'|'.join(parsedKmers)}\n{rcSeq}\n")
            probeIdx += 1
        probe1Writer.close()
        rcprobe1Writer.close()

    def makeProbe2(self):
        uniqComMap = ProbeitUtils.defineFile(self.cmDir2, 'uniq.genmap.csv')
        match self.cmError2:
            case 0:
                ProbeitUtils.simpleComputeMappability(self.window2FASTA, self.lookup2, self.window2PosBED, uniqComMap, self.pLen2, improperKmers=self.impKmers2)
            case _:
                ProbeitUtils.computeMappability(self.window2FASTA, self.idxDir2, self.cmError2, self.pLen2, self.cmDir2, uniqComMap, threads=self.threads, improperKmers=self.impKmers2)

        scBed = ProbeitUtils.defineFile(self.scDir2, 'result.bed')
        msg = ProbeitUtils.makeProbe(
            self.tempProbe2, scBed, self.window2FASTA, self.lookup2, self.pLen2,
            self.scCoverage2, self.scEarlyStop2, self.scScore2, self.scRepeats2, uniqComMap, 
            overlap=True
        )
        self.logUpdate(msg, False)
        probeIdx = 0
        probe2Writer = open(self.probe2, 'w')
        rcProbe2Writer = open(self.rcProbe2, 'w')
        for h, seq in sorted([(h,s) for h, s in SimpleFastaParser(open(self.tempProbe2))]):
            rcSeq = Seq(seq).reverse_complement()
            keys, _ = ProbeitUtils.parseGenmapPattern(h)
            keys = [self.probe1Index[k] for k in sorted(keys)]
            probe2Writer.write(  '>cap_{}\t|{}\n{}\n'.format(probeIdx, '|'.join(keys), seq) )
            rcProbe2Writer.write('>cap_{}\t|{}\n{}\n'.format(probeIdx, '|'.join(keys), rcSeq) )
            probeIdx += 1
        probe2Writer.close()
        rcProbe2Writer.close()
    
    def run(self):
        # MAKE PROBE1
        self.logUpdate("[INFO] make 1st probes")
        
        # REMOVE REDUNDANCY FROM INPUT GENOME
        if self.doRemoveRedundancy:
            self.logUpdate('[INFO]deduplicate positive fasta')
            clusteredGenome = ProbeitUtils.defineFile(self.inputDir1, 'dedup') 
            tempDir = ProbeitUtils.defineDirectory('temp', make=False, root=self.inputDir1)
            msg, self.genomeFASTA= ProbeitUtils.clusterGenome(self.inputGenome, clusteredGenome, tempDir, self.cluIdentity, threads=self.threads)
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
        self.makeProbe1()
        if not self.needProbe2:
            print('COMPLETE!!!')
            return
        
        # MAKE PROBE2
        self.logUpdate("[INFO] make 2nd probes")
            
        # MAKE WINDOW FOR PROBE2
        self.makeWindow2()
            
        # THERMO FILTER FOR PROBE2
        if self.doThermoFilter2: 
            ProbeitUtils.extractKmers(self.window2FASTA, self.kmers2FASTA, self.pLen2)
            thermoFilter2 = ThermoFilter(self.kmers2FASTA, self.pLen2, self.thermoCalc2TSV)
            msg, self.impKmers2 = thermoFilter2.run()
            self.logUpdate(msg)

        # COMPLETE MAKE PROBE2
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
        for opt, val in args:
            if opt in ('-h', '--help'):
                self.printUsage()
            try:
                # required
                self.refGenome = val if opt in ('-r', '--reference') else self.refGenome
                self.strGenome = val if opt in ('-s', '--strain') else self.strGenome
                self.workDir = val if opt in ('-o', '--output') else self.workDir
                self.posList = ProbeitUtils.getArgList(val, isInt=True) if opt in ('-p', '--positions') else self.posList
                self.snpList = ProbeitUtils.getArgList(val) if opt in ('-m', '--mutations') else self.snpList
                self.refGenomeAnnot = val if opt in ('-a', '--annotation') else self.refGenomeAnnot
                # optional
                self.needProbe2 = False if opt == '--not-make-probe2' else self.needProbe2
                self.doThermoFilter2 = False if opt == '--not-thermo-filter' else self.doThermoFilter2
                self.threads = int(val) if opt == '--threads' else self.threads
                self.isMaxWindow = True if opt == '--max-window' else self.isMaxWindow
                self.windowSize = int(val) if opt == '--window-size' else self.windowSize
                self.pLen1 = int(val) if opt == '--probe1-len' else self.pLen1
                self.pLen2 = int(val) if opt == '--probe2-len' else self.pLen2
                self.cmError2 = int(val) if opt == '--probe2-error' else self.cmError2
                self.scCoverage2 = int(val) if opt == '--probe2-cover' else self.scCoverage2
                self.scEarlyStop2 = float(val)/100 if opt == '--probe2-earlystop' else self.scEarlyStop2
                self.scRepeats2 = int(val) if opt == '--probe2-repeat' else self.scRepeats2
                # hidden args
                self.scScore2 = int(val) if opt == '--probe2-dist' else self.scScore2
                self.searchKmer = int(val) if opt == '--search-kmer' else self.searchKmer
            except Exception as e:
                print(e)
                # print("Your arguments: probeit snp {}".format(ProbeitUtils.getUserArgs(args)))
                self.printUsage() 
        #
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
            else:
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
            self.maskedGenomeFASTA= ProbeitUtils.defineFile(self.inputDir2, 'masked.fa')
            self.window2PosBED = ProbeitUtils.defineFile(self.inputDir2, 'window2.bed')
            self.window2FASTA = ProbeitUtils.defineFile(self.inputDir2, Config.window)
            self.kmers2FASTA = ProbeitUtils.defineFile(self.thermoFilterDir2, 'kmers.fa' )
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
    def _parseSearchResult(searchResult):
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
        annotationFile = open(annotation).readlines()
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
                if result.empty:
                    return result
                return result

            case _:
                return pd.DataFrame()

    def findSNPsAA(self, orf, mutation, snp, refSeq, maxPos, minPos):
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
        seqWithSNP = refSeq[codonStartPos - (maxPos-1): codonEndPos + (self.pLen1 - minPos)]
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
        window1 = [f'>probe{i}_{p.ntSnp}{ "=" + p.aaSnp if p.aaSnp else ""}\n{p.stSeq}\n' for i, p in enumerate(self.probesByPos[-1])]           
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
            ProbeitUtils.simpleComputeMappability(self.window2FASTA, self.lookup, self.window2PosBED, uniqComMap, self.pLen2, improperKmers=self.impKmers2)
        else:
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
    def printUsage():
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


class Primer:
    args = []
    shortParams = 'hp:n:o:'
    longParams = [
        # usage
        'help',
        # required
        'positive=', 'negative=', 'output=',
        # optional
        'remove-redundancy', 'not-thermo-filter', 'primer-len=', 'max-amp-len=', 'min-amp-len=', 'threads=', 'primer-cover=',
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
    ridNegId = 0.90   # for mmseqs easy-search
    cmError1 = 0
    scCoverage1 = 30
    doRemoveRedundancy = False
    doThermoFilter1 = True # thermo-filter
    isLigationProbe = False # thermo-filter
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
            if opt in ('-h', '--help'):
                self.printUsage()
            try:
                self.inputGenome = str(val) if opt in ('-p', '--positive') else self.inputGenome
                self.negGenome = str(val) if opt in ('-n', '--negative') else self.negGenome
                self.workDir = str(val) if opt in ('-o', '--output') else self.workDir
                # optional args
                self.threads = int(val) if opt == '--threads' else self.threads
                self.doRemoveRedundancy = True if opt == '--remove-redundancy' else self.doRemoveRedundancy
                self.doThermoFilter1 = False if opt == '--not-thermo-filter' else self.doThermoFilter1
                self.pLen1 = int(val) if opt == '--primer-len' else self.pLen1
                self.maxAmpLen = int(val) if opt == '--max-amp-len' else self.maxAmpLen
                self.minAmpLen = int(val) if opt == '--min-amp-len' else self.minAmpLen
                self.primerCoverage = int(val) if opt == '--primer-cover' else self.primerCoverage
                self.cmError1 = int(val) if opt == '--error' else self.cmError1
                # HIDDEN args
                self.scCoverage1 = int(val) if opt == '--probe1-cover' else self.scCoverage1
                self.scRepeats1 = int(val) if opt == '--probe1-repeat' else self.scRepeats1
                self.scEarlyStop1 = float(val)/100 if opt == '--probe1-earlystop' else self.scEarlyStop1
                # identity in mmseqs cluster
                self.cluIdentity = float(val) if opt == '--dedup-id' else self.cluIdentity
                # identity in mmeqs search
                self.ridNegId = float(val) if opt == '--rid-neg-id' else self.ridNegId
                # setcover similarity
                self.scScore1 = int(val) if opt == '--probe1-dist' else self.scScore1
            except Exception as e:
                self.printUsage()

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
        self.genomeFASTA  = self.inputGenome 
        self.window1FASTA = ProbeitUtils.defineFile(self.inputDir1, Config.window)
        self.window1PosBED = ProbeitUtils.defineFile(self.inputDir1, 'window.bed') 
        self.posKmers1FASTA = ProbeitUtils.defineFile(self.inputDir1, 'kmers.fa')
        self.negRemPosBED = ProbeitUtils.defineFile(self.maskingDir, 'negRemPos.bed')
        self.thermoCalc1TSV = ProbeitUtils.defineFile(self.thermoFiltering1Dir, "thermo_uncut.tsv")
        self.scPosBed1 = ProbeitUtils.defineFile(self.scDir1, 'result.bed')       
        self.tempProbe1 = ProbeitUtils.defineFile(self.workDir, 'temp1.fa')
        self.primerFASTA = ProbeitUtils.defineFile(self.workDir, Config.primer)
        # self.rcProbe1 = ProbeitUtils.defineFile(self.workDir, Config.rcProbe1)
        self.logUpdate('[INFO]Your arguments: probeit primer ' + ' '.join(['{} {}'.format(i[0], i[1]).strip() for i in args]), False)

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
        self.logUpdate(ProbeitUtils.ridNegKmers(self.posKmers1FASTA, self.negGenome, negKmerResult, self.maskingDir, self.ridNegId, self.threads))
        
        # MAKE DEDUPLICATED POSITIVE GENOME COORDINATE FILE(BED)
        with open(self.window1FASTA) as f:
            with open(self.window1PosBED, 'w') as w:
                idx = 0
                for _, seq in SimpleFastaParser(f):
                    w.write(f'{idx}\t0\t{str(len(seq.strip()))}\n')
                    idx += 1
        
        # EXTRACT NEGATIVE REMOVED POSITIONS BED
        w = open(negKmerPosBED, 'w')
        for i in open(negKmerResult):
            kmers = ProbeitUtils.parseKmers(i)[1:]
            for kmer in kmers:
                g, s = [int(i) for i in kmer.split(',')]
                w.write(f'{g}\t{s}\t{s+40}\n')       
        w.close()
        self.logUpdate(ProbeitUtils.getSubtractedBed(self.window1PosBED, negKmerPosBED, self.negRemPosBED), False)

    def easyComMap(self, uniqComMap):
        w = open(uniqComMap, 'w')
        alphabets = {'A', 'C', 'G', 'T'}
        self.lenSeq.clear()
        idx = 0
        for h, s in SimpleFastaParser(open(self.window1FASTA)) :
            self.lenSeq[str(idx)] = len(s)
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
                for i in range(prevEnd, self.lenSeq[prevGenome]):
                    negativeKmers.add(f'{prevGenome},{i}')
                prevGenome = g
                prevEnd = e
                for i in range(0,s):
                    negativeKmers.add(f'{g},{i}')
        else:
            for i in range(prevEnd, self.lenSeq[prevGenome]):
                negativeKmers.add(f'{prevGenome},{i}')

        for h, s in SimpleFastaParser(open(self.posKmers1FASTA)):
            kmers = set(ProbeitUtils.parseKmers(h))
            isThermoImproper = kmers&self.impKmers1
            notOnlyATGC = set(s) != alphabets&set(s)
            # TODO
            # isNegative = kmers&negativeKmers
            isNegative = {ProbeitUtils.parseKmers(h)[0]} & negativeKmers
            if isNegative or isThermoImproper or notOnlyATGC:
                continue
            w.write(h+'\n')
        w.close()
    
    def makeWindow1(self):
        ProbeitUtils.sortFasta(self.genomeFASTA)
        ProbeitUtils.simplifyFastaHeaders(self.genomeFASTA, self.window1FASTA)
        ProbeitUtils.makeLookup(self.window1FASTA, self.lookup1, self.window1PosBED)

    def makePrimers(self):
        def kmerParser(rawKmer):
            [g,p] = rawKmer.split(',')
            return int(g),int(p)

        def counter(l, key, cov):
            cnt = 0
            while cnt<cov:
                # p = -1
                try:
                    p = l.index(key)
                except ValueError:
                    return False
                l = l[p+1:]
                cnt +=1
            else:
                return True

        uniqComMap = ProbeitUtils.defineFile(self.cmDir1, 'uniq.genmap.csv')
        self.easyComMap(uniqComMap)
        msg = ProbeitUtils.makeProbe(
            self.tempProbe1, self.scPosBed1, self.window1FASTA, self.lookup1, self.pLen1, 
            self.scCoverage1, self.scEarlyStop1, self.scScore1, self.scRepeats1, uniqComMap,  
        )
        self.logUpdate(msg, False)        

        # PRIMER
        numGenomes = int([i.split()[0] for i in open(self.lookup1)][-1]) + 1
        genomeAndPairs = {i:[] for i in range(numGenomes)}
        kmers = dict()
        seqs = dict()
        idx = 0
        for h, s in SimpleFastaParser(open(self.tempProbe1)):
            for kmer in h.split()[1].split('|'):
                kmers[kmerParser(kmer)] = idx
                seqs[idx] = s
            idx += 1
        
        for i in range(numGenomes):
            primerCandidates = sorted([(g,p) for g,p in kmers.keys() if g==i])
            for j, cand1 in enumerate(primerCandidates):
                pos1 = cand1[1]
                for cand2 in primerCandidates[j+1:]:
                    pos2 = cand2[1]
                    ampLen = pos2 - pos1 + self.pLen1
                    if self.minAmpLen <= ampLen <= self.maxAmpLen:
                        genomeAndPairs[i].append((kmers[cand1],kmers[cand2]))
                  
        # primers
        with open(self.lookup1)as f:
            genomeKeys = {int(i.split()[0]):i.split()[1] for i in f}
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
        for size in sorted(pairsCoverage)[::-1]:
            for pair in pairsCoverage[size]:
                genomes = pairsAndCoveredGenomes[pair]
                for g in genomes:
                    if not counter(coveredGenomes, g, self.primerCoverage):
                        coveredGenomes.append(g)
                        primerPairsSet.add(pair)
                    if len(coveredGenomes)==numGenomes:
                        break
        idx = 0
        with open(self.primerFASTA, 'w') as w:
            for i, j in primerPairsSet:
                # primer1 = [k for k, v in keys.items() if v==i]
                # primer2 = [k for k, v in keys.items() if v==j] 
                genomes = [genome for genome, primerList in genomeAndPairs.items() if (i,j) in primerList]
                # print(genomes)
                genomes = [genomeKeys[i] for i in genomes]
                genomes = ';'.join(genomes)
                # w.write(f'>primer_pair{idx}\n{seqs[i]}\n{reverse_complement(seqs[j])}\n>primer_pair{idx}\n{seqs[j]}\n{reverse_complement(seqs[i])}\n') 
                w.write(f'>primer_pair{idx}\t{genomes}\n{seqs[i]}\n{reverse_complement(seqs[j])}\n') 
                idx += 1

    def run(self):
        # MAKE PROBE1
        self.logUpdate("[INFO] make 1st probes")
        
        # REMOVE REDUNDANCY FROM INPUT GENOME
        if self.doRemoveRedundancy:
            self.logUpdate('[INFO]deduplicate positive fasta')
            clusteredGenome = ProbeitUtils.defineFile(self.inputDir1, 'dedup') 
            tempDir = ProbeitUtils.defineDirectory('temp', make=False, root=self.inputDir1)
            msg, self.genomeFASTA= ProbeitUtils.clusterGenome(self.inputGenome, clusteredGenome, tempDir, self.cluIdentity, threads=self.threads)
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
    def printUsage():
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

    

