#!/bin/bash
#
#  Same Species Lift Over procedure, the first step: run.blat procedure
#  query and target sequences are chunked up into smaller pieces
#  each query chunk is run with blat against each target chunk
#  the result of this is a number of psl output files which are
#  chained together in the run.chain second procedure.

#  ASSUMPTION HERE ! - this is for simple small genome sequences,
#        small as in several 10's of million bases.  Larger genome sequences
#        require a cluster computer to run 10's of thousands of blat jobs.
#        and some different parameters on the chunking operation, not
#        listed here.

#  Buyer Beware - this was written without being tested, it may have errors.
#        The general theme and procedure is correct, some specifics may be off.

#       exit on any error
set -beEu -o pipefail

# kent source programs required:
#   twoBitInfo twoBitToFa partitionSequence.pl gensub2 blat calc

# use 2bit files for your target and query sequences:
#  e.g.: ${targetDb}.2bit and ${queryDb}.2bit
# 2bit files are the most convenient format for this type of work.

if [ $# -eq 0 ]
  then
    echo "Usage: liftover_chain-gen1.sh work_dir 2bit_dir targetdb querydb"
    exit
fi

#######################################################################
# First part of the procedure is a blat run between target and query chunks
# decide on a work-directory where this will take place:
export workDirectory=$1
# two directories will be created here where work-steps take place:

mkdir -p ${workDirectory}/run.blat ${workDirectory}/run.chain

# the resulting liftOver chain file will end up in ${workDirectory}/

# location where your two .2bit files are located:
export twoBitDirectory=$2

#######################################################################
# targetDb is the sequence you have good coordinates you want to convert
export targetDb=$3
# typical chunk size for a target sequence is 10,000,000 bases
export targetChunkSize="10000000"
# no need to overlap chunks, either target or query:
export chunkOverlapSize="0"
# allow only a single sequence to exist in a single chunk group
export chunkCountLimit="1"

# queryDb is the sequence you want to convert coordinates from the target
#    to the query
export queryDb=$4
# typical chunk size for a query sequence is also 10,000,000 bases
export queryChunkSize="10000000"

#######################################################################
# Compute partitions (coordinate ranges) for cluster jobs.
# working in the run.blat directory

cd ${workDirectory}/run.blat

# directory tParts will be created, make sure it is clean
rm -fr tParts
twoBitInfo ${twoBitDirectory}/${targetDb}.2bit stdout | sort -k2nr > ${targetDb}.chrom.sizes
# check the usage message for this command in the kent source tree
#   to see what all the arguments mean:
partitionSequence.pl ${targetChunkSize} \
   ${chunkOverlapSize} ${twoBitDirectory}/${targetDb}.2bit ${targetDb}.chrom.sizes \
   ${chunkCountLimit} \
   -lstDir=tParts > t.lst
# output t.lst is the listing of the target chunks to work with

# same procedure for the query bits:
rm -fr qParts
twoBitInfo ${twoBitDirectory}/${queryDb}.2bit stdout | sort -k2nr > ${queryDb}.chrom.sizes
# check the usage message for this command in the kent source tree
#   to see what all the arguments mean:
partitionSequence.pl ${queryChunkSize} \
   ${chunkOverlapSize} ${twoBitDirectory}/${queryDb}.2bit ${queryDb}.chrom.sizes \
   ${chunkCountLimit} \
   -lstDir=qParts > q.lst
# output q.lst is the listing of the query chunks to work with

#######################################################################
# Construct the 11.ooc file for the target sequence, assuming

cd ${workDirectory}/run.blat

#  the printout of the faSize is a line something like:
# 23332831 bases (10 N's 23332821 real 18400480 upper 4932341 lower) in 16 sequences in 1 files
# use the fifth field there, the 23332821 "real" number of bases:
export realSize=`twoBitToFa ${twoBitDirectory}/${targetDb}.2bit stdout | faSize stdin | grep real | awk '{print $5}'`
# calculate repMatch based on a proportion of the hg19 genoe size
export tmp=`calc $realSize / 2897316137 | sed 's/\ /\t/g' | cut -f6`
export repMatch=`calc $tmp \\\* 1024`
# That gives a really small number for small genomes, ignore that result
# and simply use a repMatch of 100:
blat ${twoBitDirectory}/${targetDb}.2bit /dev/null \
        /dev/null -tileSize=11 -makeOoc=11.ooc -repMatch=100

#######################################################################
# Construct the result output PSL directories

cd ${workDirectory}/run.blat
mkdir ${workDirectory}/run.blat/psl
for F in `cat t.lst`
do
  B=`basename ${F}`
  mkdir ${workDirectory}/run.blat/psl/${B}
done

#######################################################################
# construct gensub2 template file
#  ( allow only ${workDirectory} to shell expand here, the other $ variables
#	are for gensub2 use, hence \$ to keep them literal )
# the "blatJob.csh" script can be found at:

cat << EOF > template
#LOOP
bsub -o blat_jobs.log -q hour blatJob.csh ${workDirectory}/run.blat/\$(path1) ${workDirectory}/run.blat/\$(path2) ${workDirectory}/run.blat/psl/\$(file1)/\$(file2).psl
#ENDLOOP
EOF

## construct jobList for each query chunk blat to each target chunk:

gensub2 t.lst q.lst template jobList

# the resulting jobList is a listing of commands to be run which will perform
# the blat on each specified target/query chunk pair of sequences
# With the blatJob.csh in place in this run.blat directory, you can simply
# run each job in the jobList:
chmod +x jobList
./jobList

#######################################################################