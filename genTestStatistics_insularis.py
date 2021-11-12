import os
import csv
import sys
from scipy.stats import binom_test


#-------------------------------------------------------------------------------
def main(mapfilebase,datadir,dtotal,ctotal):

    total  = dtotal + ctotal
    pd = dtotal/total
    pc = ctotal/total
    os.chdir(datadir)
    for ctgno in range(1,8):
        for sg in ['D','C']:
            trgctg = "chr%i%s" %(ctgno,sg)
            fn = "%s.%s.simplemap.txt" % (mapfilebase, trgctg)
            fhd = open(fn)
            reader = csv.reader(fhd, delimiter='\t')
            outfile = "%s.%s.mapstats.txt" % (mapfilebase, trgctg)
            with open(outfile,'w') as out:
                positions = []
                counter = 0
                for p, cl in reader:
                    positions.append((int(p),int(cl)))
                    # for insularis we require 2000 consecutive kmers
                    if len(positions) >= 2000:
                        beg = positions[0][0]
                        end = positions[-1][0]
                        c = sum(1 for x in positions if x[1] == 0)  # label 0 for C subgenome
                        d = sum(1 for x in positions if x[1] == 1)  # label 1 for D subgenome
                        N = d + c
                        obsd = binom_test(d, n=N, p=pd, alternative='greater')
                        obsc = binom_test(c, n=N, p=pc, alternative='greater')
                        out.write("%i\t%i\t%i\t%i\t%.6e\t%.6e\n" %(beg,end,c,d,obsc,obsd))
                        positions = positions[1000:]
                        counter += 1


#-------------------------------------------------------------------------------
if __name__ == '__main__':

    mapfile_base = sys.argv[1]      # basenames for simplified comined vmatch result file
    mapfile_dir  = sys.argv[2]      # base directory for mapfiles and output stat files
    dtotal = int(sys.argv[3])       # total number of mapped D-specific unique kmers
    ctotal = int(sys.argv[4])       # total number of mapped C-specific unique kmers
    main(mapfile_base,mapfile_dir,dtotal,ctotal)