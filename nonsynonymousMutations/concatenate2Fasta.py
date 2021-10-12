import sys
def concatenate2Fasta(fofn):
    infile = open(fofn,'r')
    seqDict = {}
    for line in infile:
        newLine = line
        while newLine[-1] == '\t' or newLine[-1] == '\n' or newLine[-1] == '\r':
            newLine = newLine[0:-1]
        currSeqDict,currSeqList = buildSeqDict(newLine)
        for seq in currSeqList:
            if seq in seqDict:
                currSeq = seqDict[seq] + currSeqDict[seq]
                seqDict[seq] = currSeq
            else:
                seqDict[seq] = currSeqDict[seq]
    infile.close()
    for seq in seqDict:
        sys.stdout.write('>' + seq + '\n' + seqDict[seq] + '\n')
                
        

def buildSeqDict(fasta):
    infile = open(fasta,'r')
    scaffoldDict = {}
    scaffoldList = []
    seqName = ''
    currSeq = ''
    for line in infile:
        if line[0] == '>':
            if seqName != '':
                scaffoldDict[seqName] = currSeq
            seqName = line
            while seqName[-1] == '\n' or seqName[-1] == '\t' or seqName[-1] == '\r':
                seqName = seqName[0:-1]
            while seqName[0] == '>':
                    seqName = seqName[1:]
            scaffoldList.append(seqName)
            currSeq = ''
        else:
            currSeq += line
            while currSeq[-1] == '\n' or currSeq[-1] == '\t' or currSeq[-1] == '\r':
                currSeq = currSeq[0:-1]
    scaffoldDict[seqName] = currSeq 
    return scaffoldDict, scaffoldList



concatenate2Fasta(sys.argv[1])
