import os
import sys


def propN(fasta):
    infile = open(fasta,'r')
    lines = infile.readlines()
    infile.close()
    i = 0
    nDict = {}
    while i < len(lines):
        sN = lines[i]
        sN = sN[1:]
        while sN [-1] == '\n' or sN[-1] == '\t' or sN[-1] == '\r':
            sN = sN[0:-1]
        seq = lines[i+1]
        while seq[-1] == '\t' or seq[-1] == '\r' or seq[-1] == '\n':
            seq = seq[0:-1]
        nDict[sN] = seq.count('N')/float(len(seq))
        i += 2
    return nDict


def numStops(AAfasta):
    infile = open(AAfasta,'r')
    lines = infile.readlines()
    infile.close()
    i = 0
    stopDict = {}
    while i < len(lines):
        sN = lines[i]
        sN = sN[1:]
        while sN [-1] == '\n' or sN[-1] == '\t' or sN[-1] == '\r':
            sN = sN[0:-1]
        seq = lines[i+1]
        while seq[-1] == '\t' or seq[-1] == '\r' or seq[-1] == '\n':
            seq = seq[0:-1]
        stopDict[sN] = seq.count('*')
        i += 2
    return stopDict

fofn = open(sys.argv[1],'r')
fofnLines = fofn.readlines()
fofn.close()
seqList = []
fileName = fofnLines[0]
infile = open(fileName[0:-1],'r')
for line in infile:
    if line[0] == '>':
        sN = line[1:]
        while sN [-1] == '\n' or sN[-1] == '\t' or sN[-1] == '\r':
            sN = sN[0:-1]
        seqList.append(sN)

infile.close()
propNFile = open('Agenome.proportionNs.txt','w')
propNFile.write('Gene\tmeanPropN')
stopFile = open('Agenome.numStops.txt','w')
stopFile.write('Gene\tNumSeqs_with_Stop')
for seq in seqList:
    propNFile.write('\t' + seq)
    stopFile.write('\t' + seq)
propNFile.write('\n')
stopFile.write('\n')
for f in fofnLines:
    fileSplit = f.split('.')
    currNDict = propN(f[0:-1])
    currStopDict = numStops(fileSplit[0] + '.' + fileSplit[1] + '.cds_aaSeqs.fasta')
    geneSplit = fileSplit[0].split('/')
    meanPropN = 0
    numSeqsWStop = 0
    for seq in seqList:
        meanPropN += currNDict[seq]
        if currStopDict[seq] > 0:
            numSeqsWStop += 1
    meanPropN = round(meanPropN/len(seqList),4)
    propNFile.write(geneSplit[1] + '.' + fileSplit[1] + '\t' + str(meanPropN))
    stopFile.write(geneSplit[1] + '.' + fileSplit[1] + '\t' + str(numSeqsWStop))
    for seq in seqList:
        propNFile.write('\t' + str(round(currNDict[seq],4)))
        stopFile.write('\t' + str(currStopDict[seq]))
    propNFile.write('\n')
    stopFile.write('\n')
    
    



        

