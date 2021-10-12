def reverseComplement(seq):
    seq_revc = ''
    for nuc in seq:
        if nuc == 'A':
            seq_revc = 'T' + seq_revc
        elif nuc == 'T':
            seq_revc = 'A' + seq_revc
        elif nuc == 'C':
            seq_revc = 'G' + seq_revc
        elif nuc == 'G':
            seq_revc = 'C' + seq_revc
        else:
            seq_revc = nuc + seq_revc
    return seq_revc

def fileInputRC(Infile):
    #File must be single-line formatted fasta 
    infile = open(Infile,'r')
    outfile = open(Infile[0:-6] + '_RC.fasta','w')
    for line in infile:
        if line[0] == '>':
            outfile.write(line)
        else:
            currLine = line
            while currLine[-1] == '\n' or currLine[-1] == '\t' or currLine[-1] == '\r':
                currLine = currLine[0:-1]
            seq = reverseComplement(currLine)
            outfile.write(seq + '\n')
    infile.close()
    outfile.close()