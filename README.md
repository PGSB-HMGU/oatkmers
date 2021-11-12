# 0. Introduction

Directory contains workflow and scripts applied for the detection of sub-genome specific kmers in oat
as described in our oat reference genome manuscript (submitted). Workflow is provided in this README.

Required software tools (third party, please install them according to instructions provided at respective
links; e.g. by a conda environment):

1. python 3.7 (or higher, tested with 3.7)
2.  KMC and KMC tools\
http://sun.aei.polsl.pl/REFRESH/index.php?page=projects&project=kmc&subpage=about\
https://github.com/refresh-bio/KMC
3. vmatch\
http://www.vmatch.de/

You will also need the following python site packages:
- scipy
- numpy
- matplotlib
- biopython

You will also need access to a high end workstation, in our case we used either a slurm batch queue, or a
server with 80 cores and 1 TB of RAM (20-40 cores and 100 Gb should be sufficient)


# 1. Generation of subgenome-specific kmer databases by KMC/KMC tools

Here, only the workflow for the D/C-subgenome specific kmers of A.insularis are described. For hexaploid
oat, steps are analogous to this workflow, a full description of the applied set operations can be found
in the manuscript, see Supplementary Methods.

`$genome` path variable to genome fasta sequence path

`$subgenome` path variable to all fasta sequences of one subgenome (eg,. D,C for insularis, A,C,D for oat)

`$workdir` path to temporary working directory

`$outbase` path to basename for output databases

`$KMER` kmer size, for manuscript k=27


## a. Database generation for each genome/subgenome:

```
kmc -m30 -t10 -k${KMER} -fm -ci1 -cs1000000000 $genome $outbase.k$KMER $workdir
rm -rf $workdir
```

Given three Avena species, _A.longiglumis_ (A-genome), _A.eriantha_ (C-genome) and _A.insularis_ (D(A)- and C-genome),
three databases are generated, `longiglumis.k27`, `eriantha.k27` and `insularis.k27`

## b. Subgenome kmer lists

`$def_file` (example provided for D-specific kmers of _A.insularis_):
```
INPUT:
lon = longiglumis.k27
eri = eriantha.k27
ins = insularis.k27
OUTPUT:
insD.k27 = (ins - eri) * left lon
OUTPUT_PARAMS:
-cx1
```


```
kmc_tools complex $def_file
kmc dump insD.k27 insD.k27.txt
# dumps text file of kmers generated from def_file rules
python kmerlist2fasta.py $insD.k27.txt $insD.k27.fasta
```


# 2. Positioning of unique subgenome specific kmers

First generate an suffixarray index of the genome, then map subgenome specific kmer fasta list to it.
For mapping, we recommend to split the kmer fastas and run batches in a queue, then cat vmatch result
files.

`$kmerfasta` fasta file containing specific kmers from above

`$outfile` vmatch output of mapping `$kmerfasta` to index

## make genomic index
```
mkvtree -db $genome -dna -indexname ${genome}.vmidx -allout
```

## map kmer list to genome/index
```
vmatch -showdesc 20 -d -p -complete -noidentity -noevalue -nodist -q $kmerfasta ${genome}.vmidx > $outfile
```


next, we simplify and combine subgenome-specific kmer mappings into position ordered hit files
for this we generate a tab-delimited `config_file` that list for each subgenome specific mapping the following infos:
```
path_vmatch_output_C-specific_kmers     integer_label[eg. 0]
path_vmatch_output_D-specific_kmers     integer_label[eg. 1]
path_vmatch_output_A-specific_kmers     integer_label[eg. 2]  # optional for hexaploid oat
```
consider `$input_file` the config file above, and `$outdir`, `$simple_outfile_base` the path definitions you store
the combined mapping results per contig ID (for convenience) in, then run

```
python orderedKmermappings.py $input_file $outdir $simple_outfile_base
```


# 3. Window-based kmer statistics

To get the statistics for a potential subgenome specific kmer overrepresentation in a particular genomic region,
binomial testing within a user-defined range of consecutive kmers (2000 in the case of _A.insularis_, and 5000 for
_A.sativa_). Two scripts are provided, one for _A.insularis_, and one for hexaploid oat, each gets the respective path
definitions as in `orderedKmermappings.py` (see above), a/c/d-totals are provided in the Supplement of manuscript:

```
python genTestStatistics_insularis.py $mapfiles_base $mapfiles_dir $dtotal $ctotal
python genTestStatistics_sativa.py $mapfiles_base $mapfiles_dir $atotal $dtotal $ctotal
```





