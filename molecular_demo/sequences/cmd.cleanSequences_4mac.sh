export seqFile=$1

gsed -i "/>/d; /structureX/d" $seqFile
tr -d '\n' < $seqFile  > a
tr "*" "\n" < a > b
mv b $seqFile
# Get rid of the chain separating line;
gsed -i "s/\///g" $seqFile

