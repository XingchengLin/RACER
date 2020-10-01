####################################################################################
# This script will help add fake CB atoms for Glycines
#
# Written by Xingcheng Lin, 06/10/2020
####################################################################################

import math
import subprocess
import os
import time
import sys

import numpy as np


from common_function import *

################################################


def my_lt_range(start, end, step):
    while start < end:
        yield start
        start += step


def my_le_range(start, end, step):
    while start <= end:
        yield start
        start += step


def RepresentsFloat(s):
    try:
        float(s)
        return True
    except ValueError:
        return False
###########################################



os.chdir('./test_structures_pdbs_with_virtual_cbs')


add_virtual_glycines_list("testSetFiles.txt")
