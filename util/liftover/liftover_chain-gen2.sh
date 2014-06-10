#!/bin/bash
#
#  Same Species Lift Over procedure, this second step: run.chain procedure
#  using the PSL file results from the run.blat first step, chain
#  those results into raw chain files in preparation for the net step

#  Buyer Beware - this was written without being tested, it may have errors.
#        The general theme and procedure is correct, some specifics may be off.

#       exit on any error
set -beEu -o pipefail

# kent source programs required:
#   gensub2 axtChain

# use 2bit files for your target and query sequences:
#  e.g.: ${targetDb}.2bit and ${queryDb}.2bit
# 2bit files are the most convenient format for this type of work.
if [ $# -eq 0 ]
  then
    echo "Usage: liftover_chain-gen1.sh work_dir 2bit_dir targetdb querydb"
    exit
fi
#######################################################################

# same workDirectory as from the run.blat first step
export workDirectory=$1
# the run.chain directory ${workDirectory}/run.chain
# was already created by the run.blat first step script
# the resulting liftOver chain file will end up in ${workDirectory}/

# location where your two .2bit files are located:
export twoBitDirectory=$2

#######################################################################
# same target and query sequences used in the run.blat first step
# targetDb is the sequence you have good coordinates you want to convert
export targetDb=$3

# queryDb is the sequence you want to convert coordinates from the target
#    to the query
export queryDb=$4

#######################################################################
# simple construction of the parts list, is merely a listing of
#   the files from the run.blat/psl/ directory hierarchy
#   (this only works because the parts are all separate chromosome bits,
#   there is a more complicated script in the source tree that can manage
#   chunks for large genomes where chromosomes are split in chunks)
cd ${workDirectory}/run.chain
mkdir -p chainRaw
find ../run.blat/psl -type f | sed -e 's#../run.blat/psl/##' > pslParts.lst
# construct the individual result directories
for D in ../run.blat/psl/*
do
   B=`basename $D`
   mkdir chainRaw/${B}
done

#######################################################################
# construct gensub2 template file and the chainJob.csh script:

cat << EOF > template
#LOOP
bsub -o chain_jobs.log -q hour ${workDirectory}/run.chain/chainJob.csh \$(path1) ${workDirectory}/run.chain/chainRaw/\$(path1).chain
#ENDLOOP
EOF

cat << EOF > chainJob.csh
#!/bin/csh -ef

set inPattern = \$1
set outChain = \$2
set targetDb = "${targetDb}"
set queryDb = "${queryDb}"
set tmpOut = \`mktemp -p /dev/shm doSame.chain.XXXXXX\`

cat ../run.blat/psl/\$inPattern* \
  | axtChain -verbose=0 -linearGap=medium -psl stdin \
    ${twoBitDirectory}/${targetDb}.2bit \
    ${twoBitDirectory}/${queryDb}.2bit \
    \$tmpOut
mv \$tmpOut \$outChain
chmod 664 \$outChain
EOF

chmod +x chainJob.csh

gensub2 pslParts.lst single template jobList

# the resulting jobList is a listing of commands to be run which will perform
# the chain on each specified pslParts list
# With the chainJob.csh in place in this run.chain directory, you can simply
# run each job in the jobList:
chmod +x ./jobList
./jobList

#######################################################################
#######################################################################
#######################################################################
# setup net script

cd ${workDirectory}/run.chain

cat << EOF > doNet.csh
#!/bin/csh -efx
set targetDb = "${targetDb}"
set queryDb = "${queryDb}"
# Use local scratch disk... this can be quite I/O intensive:
set tmpDir = \`mktemp -d -p /dev/shm doSame.net.XXXXXX\`
# Merge up the hierarchy and assign unique chain IDs:
mkdir \$tmpDir/chainMerged
foreach d (chainRaw/*)
  set tChunk = \$d:t
  chainMergeSort \$d/*.chain > \$tmpDir/chainMerged/\$tChunk.chain
end
chainMergeSort \$tmpDir/chainMerged/*.chain \
| chainSplit -lump=100 \$tmpDir/chainSplit stdin
endsInLf \$tmpDir/chainSplit/*.chain
rm -rf \$tmpDir/chainMerged/
mkdir \$tmpDir/netSplit \$tmpDir/overSplit
foreach f (\$tmpDir/chainSplit/*.chain)
  set split = \$f:t:r
  chainNet \$f \
  ../run.blat/${targetDb}.chrom.sizes \
  ../run.blat/${queryDb}.chrom.sizes \
    \$tmpDir/netSplit/\$split.net /dev/null
  netChainSubset \$tmpDir/netSplit/\$split.net \$f stdout \
  | chainStitchId stdin \$tmpDir/overSplit/\$split.chain
end
endsInLf \$tmpDir/netSplit/*.net
endsInLf \$tmpDir/overSplit/*.chain
cat \$tmpDir/chainSplit/*.chain | gzip -c > ${targetDb}.${queryDb}.all.chain.gz
cat \$tmpDir/netSplit/*.net     | gzip -c > ${targetDb}.${queryDb}.noClass.net.gz

cat \$tmpDir/overSplit/*.chain | gzip -c > ../${targetDb}To${queryDb}.over.chain.gz

rm -rf \$tmpDir/
EOF

chmod +x doNet.csh
./doNet.csh > net.log 2>&1

### Note the result file is contructed by the doNet.csh script:
###   ../${targetDb}To${queryDb}.over.chain.gz

#######################################################################