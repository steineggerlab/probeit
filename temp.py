from Bio.SeqIO.FastaIO import SimpleFastaParser
import re

genomeKeys = {}
ligKeys = {}
ligNewNames = {}
with open('posnegset01/genome.lookup')as f:
    for i in f:
        genomeKeys[int(i.split()[1])] = i.split()[0]
with open('posnegset01/lig_probe.lookup')as f:
    for i in f:
        ligKeys[int(i.split()[1])] = i.split()[0]

cnt = 0
for h, s in SimpleFastaParser(open('posnegset01/probe1.fa')):
    print(h)
    p1 = re.compile('[0-9]+,')
    p2 = re.compile(',[0-9]+')

    genomes = [genomeKeys[int(i[:-1])] for i in p1.findall(h)]

    pos = [int(i[1:]) for i in p2.findall(h)]
    probes = ['{}:{}:{}'.format(genomes[i], pos[i], pos[i] + 40) for i in range(len(genomes))]
    for i in probes:
        ligNewNames[i] = 'LP{}'.format(cnt)
    print(
        'LP{}\t'.format(cnt) + ';'.join(probes)
    )

    cnt += 1