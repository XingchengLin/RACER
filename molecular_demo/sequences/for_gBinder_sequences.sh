cat native.seq > gBinder_sequences.txt 
sed -i "/>/d; /structureX/d" gBinder_sequences.txt
tr -d '\n' < gBinder_sequences.txt  > a
tr "*" "\n" < a > b
mv b gBinder_sequences.txt
sed -i "s/\///g" gBinder_sequences.txt
