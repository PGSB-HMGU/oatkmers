import os
import sys
import csv

#-------------------------------------------------------------------
def kmer2fasta(filename,outfilename):
    """transforms sequences of kmc dump to fasta formatted sequences"""
    # a few path checks
    if not os.path.isfile(filename):
        sys.stderr.write("Kmer input file does not exist\n")
        raise OSError
    if not os.path.isdir(os.path.dirname(outfilename)):
        sys.stderr.write("Output directory does not exist\n")
        raise OSError

    with open(filename) as fhd:
        reader = csv.reader(fhd,delimiter='\t')
        out = open(outfilename,'w')
        for seqid,datum in enumerate(reader,1):
            seq,cnt = datum
            if int(cnt) > 1: continue  # take only single copy kmers
            out.write(">km%i\n%s\n" %(seqid,seq))
        out.close()


#-------------------------------------------------------------------
if __name__ == '__main__':

    fn = sys.argv[1]
    outfile = sys.argv[2]
    kmer2fasta(fn,outfile)