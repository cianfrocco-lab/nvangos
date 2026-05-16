from collections import Counter
from statistics import mode 
import sys 
import shutil

### This script is for inner joining two protofilament stacks from either side of a 2pf seam stack with alternating registers.


# Set input filenames - must be in same directory as script

it00_file = 'particles.star'
it01_file = 'particles_alternate_pf.star'


# Set empty lists to fill with variables of interest later on

it00_value_list = []
it01_value_list = []


# Set output filename

it01_modified = '%s_seam_combined.star' %(it00_file[:-5])


# Copy contents of run_it001_data.star file to variable for later (this is a string)

fileopen_it01=open(it01_file,'r')
it01_contents = fileopen_it01.read()
fileopen_it01.close()


# Open input files for reading

fileopen_it00=open(it00_file,'r') 
fileopen_it01=open(it01_file,'r') 


#Get rlnAngleTiltPrior, rlnAnglePsiPrior, rlnAngleRot, rlnAngleTilt, rlnAnglePsi, rlnOriginXAngst, and rlnOriginYAngst index values

for line in fileopen_it00:
    if '_rlnImageName' in line: 
        particle_name_col = line.split('#')[-1]


# Save important baseline variables from run_it000_data.star to a list of lists for later repurposing, using previous index values for future proofing

    if '.mrc' in line:
        particle_name = line.split()[int(particle_name_col)-1]

        it00_value_list.append(particle_name)


# Close input files

fileopen_it00.close()
fileopen_it01.close()


# Split previously saved string of run_it001_data.star contents by return characters, making a list of text lines

it01_contents = it01_contents.split('\n')

new_contents = []

# Set index with which to check line number in subsequent for loop

line_index = -1


# Set counter to use to pull sequential values from value lists

counter = 0


it00_value_set = set(it00_value_list)


# Check each line of run_it001_data.star contents

for line in it01_contents:
    line_index += 1


# Set variable to work with specific line info based on iterative index number

    line_info = it01_contents[line_index]


# For every non-header line (ie. every line with real info), populate new file with  

    if '.mrc' in line:
        particle_name = line.split()[int(particle_name_col)-1]
        if particle_name in it00_value_set:
            new_contents.append(line_info)
    else:
        new_contents.append(line_info)


# Join previously split string back together using return characters, creating an identically formatted string to the original it01_contents

it01_modified_contents = '\n'.join(new_contents)


# Write updated string, which is basically the modified run_it001_data.star file, into a new file

fileopen_it01_modified = open(it01_modified,'w')

fileopen_it01_modified.write(it01_modified_contents)


# Close modified file

fileopen_it01_modified.close() 

