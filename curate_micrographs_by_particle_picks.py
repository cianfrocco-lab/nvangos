#!/usr/bin/env python
from collections import Counter
from statistics import mode 
import sys 
import shutil
import re

### This script is for curating the micrographs_ctf.star file based on micrographs that have (autopick) particle picks on them, identified via the Autopick/summary.star file.
### Run this script in the same directory as the summary.star (symlinked) and micrographs_ctf.star (original) files, ie. the CtfFind job output directory.


# Set input filenames - must be in same directory as script

it00_file = 'summary.star'            # Output from autopick job - symlink this file to the directory containing this script
it01_file = 'micrographs_ctf.star'    # Output from CTF estimation job 


# Set empty lists to fill with variables of interest later on

it00_value_list = []
it01_value_list = []


# Set output filename

it01_reset_file = '%s_CURATED.star' %(it01_file[:-5])


# Copy contents of micrographs_ctf.star file to variable for later (this is a string)

fileopen_it01=open(it01_file,'r')
it01_contents = fileopen_it01.read()
fileopen_it01.close()


# Open input files for reading

fileopen_it00=open(it00_file,'r') 
fileopen_it01=open(it01_file,'r') 


#Get _rlnMicrograph_name column # 

for line in fileopen_it00:
    if '_rlnMicrographName' in line: 
       mic_name_col  = line.split('#')[-1]


# Save names of micrographs containing picks (from summary.star) to a list for later repurposing, using previous index values for future proofing

    if '.mrc' in line:
        mic_name = line.split()[int(mic_name_col)-1]

        it00_value_list.append(mic_name)


it00_value_list.append('padding_not_in_list')

#print(it00_value_list)

# Close input files

fileopen_it00.close()
fileopen_it01.close()


# Set index with which to check line number in subsequent for loop

line_index = -1


# Set counter to use to pull sequential values from value lists

counter = 0


# Split micrographs_ctf.star string variable into list to iterate through by line

it01_contents = it01_contents.split('\n')

it01_contents_new = []

# Check each line of run_it001_data.star contents

for line in it01_contents:
    line_index += 1


# Set variable to work with specific line info based on iterative index number
    if '.mrc' not in line:
        line_info = it01_contents[line_index]
        it01_contents_new.append(line_info)

# For every line that contains a micrograph with picks (as determined by summary.star file earlier), copy that line info. This excludes lines that don't have any picks.

    if str(it00_value_list[counter]) in line:
        line_info = it01_contents[line_index]
        it01_contents_new.append(line_info) 
        counter += 1


# Join previously split string back together using return characters, creating an identically formatted string to the original it01_contents

it01_modified_contents = '\n'.join(it01_contents_new)


# Write updated string, which is basically the modified run_it001_data.star file, into a new file

fileopen_it01_reset = open(it01_reset_file,'w')

fileopen_it01_reset.write(it01_modified_contents)


# Close modified file

fileopen_it01_reset.close() 

