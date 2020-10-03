#########################################################################
# Author: Xingcheng Lin
# Created Time: Wed Jun 10 22:48:33 2020
# File Name: cmd.preprocessing.sh
# Description: 
#########################################################################
#!/bin/bash

export PDBid=$1
export alphaChain=$2
export betaChain=$3

# Add fake CB atoms, AWSEM format requirement
cp data/native.pdb native_structures_pdbs_with_virtual_cbs/
python add_fakeCB.py



# Generate randomized sequence for the decoys;
cp proteins_list.txt sequences/
cp native_structures_pdbs_with_virtual_cbs/native.pdb sequences/
cd sequences/

# Build the sequence for native.pdb
python buildseq.py native
bash cmd.cleanSequences_4mac.sh native.seq 

bash for_gBinder_sequences_4mac.sh

# Create the label for the peptide sequence indices;
#echo "Type the startind and ending indices of peptide residues, separated by space:"
#read -rsp $'Press any key to continue...\n' -n1 key

#read -r pResid_startID pResid_endID

pResid_startID=$4
pResid_endID=$5

rm randomize_position_file.txt
for ((i=$pResid_startID; i<$pResid_endID; i++))
do 
    printf "%d " $i >> randomize_position_file.txt
done
printf "%d" $pResid_endID >> randomize_position_file.txt


mkdir -p TCR_randomization 
mv randomize_position_file.txt native.seq gBinder_sequences.txt TCR_randomization/

cd ../

python generate_decoy_seq.py

# Create the tms file, where the peptide is labeled as '2', while TCRs are labeled as '1';
grep 'CA' native_structures_pdbs_with_virtual_cbs/native.pdb > tmp.txt
# Get the total number of residues;

tot_resnum=`cat tmp.txt | awk 'END{print $6}'`
python create_tms.py sequences/TCR_randomization/randomize_position_file.txt $tot_resnum



gsed "s/TCR_NAME/$PDBid/g; s/TCR_ALPHACHAIN/$alphaChain/g; s/TCR_BETACHAIN/$betaChain/g" template_evaluate_phi.py > evaluate_phi.py
python evaluate_phi.py

