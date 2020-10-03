#########################################################################
# Author: Xingcheng Lin
# Created Time: Sun Jun 14 16:58:44 2020
# File Name: cmd.sh
# Description: 
#########################################################################
#!/bin/bash

bash cmd.preprocessing_4mac.sh 3qib C D 782 794

bash cmd.optimize.sh

echo "Check the TCR chain IDs of the testBinder file (because Modeller changed it when rebuilding...) !!!"
read -rsp $'Press any key to continue...\n' -n1 key

bash cmd.evaluate_bindingE_4mac.sh 3qib C D
