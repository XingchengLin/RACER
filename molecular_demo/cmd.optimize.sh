#########################################################################
# Author: Xingcheng Lin
# Created Time: Wed Jun 17 22:29:58 2020
# File Name: cmd.process.sh
# Description: 
#########################################################################
#!/bin/bash

# The cut off for noise filtering
abs_cutoff=74

python optimize_gamma.py $abs_cutoff

