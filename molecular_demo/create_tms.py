###########################################################################
# This script will segment PDB according to user provided start and end IDs
#
# Written by Xingcheng Lin, 12/12/2016;
###########################################################################

import time
import subprocess
import os
import math
import sys
import numpy as np

################################################


def my_lt_range(start, end, step):
    while start < end:
        yield start
        start += step


def my_le_range(start, end, step):
    while start <= end:
        yield start
        start += step
#############################################


def create_tms(random_position_file, tot_resnum):

    # directory for the tms file
    tms_directory = "./tms/"
    # Get current working directory
    pwd = os.getcwd()

    infile = open(random_position_file, 'r')
    outfile = open(tms_directory + 'native.tm', 'w')

    # Read in the only one line from the random_position_file;

    lines = [line.strip() for line in infile]

    infile.close()

    random_position = lines[0].split(" ")
    random_position_int = [int(integer) for integer in random_position]
    print(random_position_int)

    for i in my_le_range(1, tot_resnum, 1):

        if (i in random_position_int):
            # There is no newline after the last residue
            if (i == tot_resnum):
                outfile.write("2")
            else:
                outfile.write("2" + "\n")
        else:
            # There is no newline after the last residue
            if (i == tot_resnum):
                outfile.write("1")
            else:
                outfile.write("1" + "\n")

    outfile.close()
    return

############################################################################

if __name__ == "__main__":
    random_position_file = sys.argv[1]
    tot_resnum = int(sys.argv[2])

    create_tms(random_position_file, tot_resnum)
