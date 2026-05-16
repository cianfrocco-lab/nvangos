#!/usr/bin/env python
from collections import Counter
from statistics import mode 
import sys 
import shutil

### This script is for creating different optics groups for each patch collected via SerialEM in the run_data.star file from a 3D Refinement job prior to the final CTF Refinement.


# Set input filenames - must be in same directory as script

data_file = 'run_data.star'


# Set empty lists to fill with variables of interest later on

image_shift_list = []


# Set output filename

data_modified_file = '%s_optics_enumerated.star' %(data_file[:-5])


# Copy contents of run_data.star file to variable for later (this is a string)

fileopen_data=open(data_file,'r')
data_contents = fileopen_data.read()
fileopen_data.close()


# Open input files for reading

fileopen_data=open(data_file,'r') 


#Get _rlnMicrographName for immediate use and _rlnOpticsGroup _rlnOpticsGroupName index values for later use

executed_once = False

for line in fileopen_data:
    if '_rlnMicrographName' in line: 
        micrograph_name_col = int(line.split('#')[-1]) - 1
        print("micrograph_name_col = " + str(micrograph_name_col))
    if '_rlnOpticsGroup ' in line and executed_once:
        optics_group_col = int(line.split('#')[-1]) - 1
        print("optics_group_col = " + str(optics_group_col))
    if '_rlnOpticsGroup ' in line and not executed_once:
        optics_group_header_col = int(line.split('#')[-1]) - 1
        executed_once = True
        print("optics_group_header_col = " + str(optics_group_header_col))
    if '_rlnOpticsGroupName' in line:
        optics_group_name_header_col = int(line.split('#')[-1]) - 1
        print("optics_group_name_header_col = " + str(optics_group_name_header_col))
    
# Close input files

fileopen_data.close()


# Split previously saved string of run_data.star contents by return characters, making a list of text lines

data_contents = data_contents.split('\n')


#print(data_contents)

# Set index with which to check line number in subsequent for loop

line_index = -1


# Check each line of run_data.star contents

for line in data_contents:
    line_index += 1


# Set variable to work with specific line info based on iterative index number

    line_info = data_contents[line_index]


# For every non-header line (ie. every line with real info), check if there is an unique image shift. If not, continue. If so, create a new optics group and optics group header for that new image shift.
    
    if ".mrc" in line:
        line_info = line_info.split()
        image_shift = line_info[micrograph_name_col][-12:-4]
        if image_shift not in image_shift_list:
            image_shift_list.append(image_shift)
            print(image_shift)
        image_shift_index = image_shift_list.index(image_shift) + 1
        line_info[optics_group_col] = str(image_shift_index)
        line_info = '   '.join(line_info)      



# Change content in actual list

    data_contents[line_index] = line_info

# Change header of updated file to reflect updated list of opticsGroups

line_index = -1

new_line = []

for line in data_contents:
    line_index += 1

    line_info = data_contents[line_index]
    
    extended_line = line_info

    if "opticsGroup" in line:
        for n in range(1,len(image_shift_list)+1):
            line_info = line_info.split()
            line_info[optics_group_name_header_col] = line_info[optics_group_name_header_col][:11] + str(n)
            line_info[optics_group_header_col] = str(n)
            line_info = '   '.join(line_info)
            new_line.append(line_info)
        line_info = '\n'.join(new_line)

    data_contents[line_index] = line_info

# Join previously split string back together using return characters, creating an identically formatted string to the original data_contents

data_modified_contents = '\n'.join(data_contents)


# Write updated string, which is basically the modified run_data.star file, into a new file

fileopen_data_modified = open(data_modified_file,'w')

fileopen_data_modified.write(data_modified_contents)


# Close modified file

fileopen_data_modified.close() 

