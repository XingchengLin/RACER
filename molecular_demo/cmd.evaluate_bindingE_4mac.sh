#########################################################################
# Author: Xingcheng Lin
# Created Time: Thu Jun 18 22:32:00 2020
# File Name: cmd.for_evaluation_E.sh
# Description: 
#########################################################################
#!/bin/bash


PDBid=$1
export alphaChain=$2
export betaChain=$3



mkdir -p test_structures_pdbs_with_virtual_cbs/
cp data/testBinder.*.pdb test_structures_pdbs_with_virtual_cbs/

rm testSetFiles.txt
for ((i=1; i<=4; i++))
do 
    echo testBinder.$i >> testSetFiles.txt
done
cp testSetFiles.txt test_structures_pdbs_with_virtual_cbs/

cd tms/
for ((i=1; i<=4; i++))
do 
    cp native.tm testBinder.$i.tm
done
cd ../

mkdir -p evaluated_binding_E

gsed "s/TCR_NAME/$PDBid/g; s/TCR_ALPHACHAIN/$alphaChain/g; s/TCR_BETACHAIN/$betaChain/g" template_evaluate_binding_E.py > evaluate_binding_E.py
python evaluate_binding_E.py
