#################################################################################
#                                                                               #
#   Script takes vmatch kmer mapping output files for subgenome specific kmers  #
#   and combines them into simplified format.                                   #
#                                                                               #
#################################################################################

import re
import os
import csv
import sys
import itertools
from operator import attrgetter
from dataclasses import dataclass

@dataclass
class KmerPosition:
    __slots__ = ["contigid","position","ancestor"]
    contigid : str
    position : int
    ancestor : int


def combine_simplify(cnffile,outdir,outbase):

    data = []
    with open(cnffile) as fhd:
        reader = csv.reader(fhd,delimiter='\t')
        for infile,label in reader:
            data.append((infile,label))

    mapping = []
    # now per vmatch output for each of the subgenome specific kmer mappings
    for fn,label in data:
        subgenome_specifity = int(label)
        with open(fn) as fhd:
            for line in fhd:
                if line.startswith("#") : continue  # this is a vmatch comment
                s,ctg,pos,*_ = re.split("\s+",line.strip())
                if ctg == "chrUn": continue
                mapping.append(KmerPosition(contigid=ctg,position=int(pos),ancestor=subgenome_specifity))

    mapping.sort(key=attrgetter('contigid','position'))
    os.chdir(outdir)
    for ctgid,g in itertools.groupby(mapping,key=attrgetter('contigid')):
        outfile = "%s.%s.simplemap.txt" %(outbase,ctgid)
        with open(outfile,'w') as out:
            for m in g:
                # we skip the contigid info, not necessary as file
                # contains only kmer positions from exactly one contigid
                out.write("%i\t%i\n" %(m.position,m.ancestor))



if __name__ == '__main__':

    argcount = len(sys.argv)
    cnffile = sys.argv[1]   # contains information for kmer mappings, see README
    outdir  = sys.argv[2]   # output directory where to save result files per contigID
    outbase = sys.argv[3]   # contains combined simplified kmer position per contigID


