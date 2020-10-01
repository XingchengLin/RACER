#########################################################################
# Author: Xingcheng Lin
# Created Time: Mon Jun  8 22:14:51 2020
# File Name: cmd.createtm.sh
# Description: 
#########################################################################
#!/bin/bash


rm native.tm

for ((i=1; i<=588; i++))
do
    echo "1" >>native.tm
done

for ((i=589; i<=604; i++))
do
    echo "2" >> native.tm
done
