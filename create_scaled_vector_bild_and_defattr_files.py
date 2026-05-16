#!/usr/bin/env python

#Import necessary modules

import math
import numpy as np
import Bio.Align
from Bio.Align import Alignment

one_letter_codes = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
     'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N', 
     'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W', 
     'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}

#Inputs:

file_1='input_file_1.pdb'                               ### CHANGE INPUT FILE NAME 1 HERE ###
file_1_number='1'
file_2='input_file_2.pdb'              ### CHANGE INPUT FILE NAME 2 HERE ###
file_2_number='2'

color = "green"                                         ### CHANGE PSEUDOBOND COLOR HERE (can also edit file header later) ###
radius = "0.3"                                          ### CHANGE PSEUDOBOND RADIUS HERE (can also edit file header later) ###
dashes = "1"

chains=['B','C','H','F','P','N','J','K']                ### CHANGE CHAINS TO ANALYZE HERE ###


#output
outfile1='ca_distances.defattr'
outfile2='colored_vectors_XYZ.bild'
outfile3='colored_vectors_XYZ.pb'
outfile4='ca_distances_XY_only.defattr'
outfile5='colored_vectors_XY_only.bild'
outfile6='ca_distances_Z_only.defattr'
outfile7='colored_vectors_Z_only.bild'


'''This script will generate .bild and .pb files with per-residue CA vectors between two models colored by orientation angle. ChimeraX can display these files.
 '''



o1 = open(outfile1,'w')
o1.write('''#  Per-residue alpha carbon rmsd values in Angstroms. These are calculated from 
#  the .pdb files themselves, and it is assumed these were aligned properly beforehand.
#
#  Use this file to color by attribute in ChimeraX.
#
attribute: ca_distances
recipient: residues
''')

o2 = open(outfile2,'w')
o2.write('''#  .bild file to display displacement vectors between residue alpha carbons 
# colored by cartesian coordinate displacement mapped to RGB.
#
#
''')

o3 = open(outfile3,'w')
o3.write("; color = "+color+"\n; radius = "+radius+"\n; dashes = "+dashes+"\n")

o4 = open(outfile4,'w')
o4.write('''#  Per-residue alpha carbon rmsd values of lateral (X & Y) displacement in 
# Angstroms. These are calculated from the .pdb files themselves, and it is assumed 
# these were aligned properly beforehand.
#
#  Use this file to color by attribute in ChimeraX.
#
attribute: ca_XY_distances
recipient: residues
''')

o5 = open(outfile5,'w')
o5.write('''#  .bild file to display displacement vectors between residue alpha carbons 
colored by lateral (X & Y) displacement.
#
#
''')

o6 = open(outfile6,'w')
o6.write('''#  Per-residue alpha carbon rmsd values of longitudinal (Z) displacement in 
# Angstroms. These are calculated from the .pdb files themselves, and it is assumed 
# these were aligned properly beforehand.
#
#  Use this file to color by attribute in ChimeraX.
#
attribute: ca_Z_distances
recipient: residues
''')

o7 = open(outfile7,'w')
o7.write('''#  .bild file to display displacement vectors between residue alpha carbons
colored by longitudinal (Z) displacement (this is directional).
#
#
''')

#Harvest maximum distance and vector components for scaling/coloring arrows.

max_distance = 0
max_distance_XY = 0
max_vector_componentX = 0
max_vector_componentY = 0
max_vector_componentZ = 0

for chain in chains:

    chain_info_file1 = {}
    file1_index = []
    file1_sequence = ''
    chain_info_file2 = {}
    file2_index = []
    file2_sequence = ''

    pdbopen=open(file_1,'r')
    
    for line in pdbopen:
        if 'ATOM' in line:
            if line.split()[4] == chain:
                atom=line.split()[2]
                if atom=='CA':
                    residue_num = line.split()[5]   
                    residue = line.split()[3]
                    positionX = float(line.split()[6])
                    positionY = float(line.split()[7])
                    positionZ = float(line.split()[8])

#                    print({residue_num:(residue,positionX,positionY,positionZ)}) # FOR TESTING

                    chain_info_file1.update({residue_num:(residue,positionX,positionY,positionZ)})

                    file1_index.append(residue_num)

                    file1_sequence = file1_sequence + one_letter_codes.get(residue)
                    
    pdbopen.close()

    pdbopen=open(file_2,'r')

    for line in pdbopen:
        if 'ATOM' in line:
            if line.split()[4] == chain:
                atom=line.split()[2]
                if atom=='CA':
                    residue_num = line.split()[5]   
                    residue = line.split()[3]
                    positionX = float(line.split()[6])
                    positionY = float(line.split()[7])
                    positionZ = float(line.split()[8])

    #                print({residue_num:(residue,positionX,positionY,positionZ)}) # FOR TESTING

                    chain_info_file2.update({residue_num:(residue,positionX,positionY,positionZ)})

                    file2_index.append(residue_num)

                    file2_sequence = file2_sequence + one_letter_codes.get(residue)
                    
    pdbopen.close()

    '''
    if int(list(chain_info_file1.keys())[-1]) > int(list(chain_info_file2.keys())[-1]):
        shortest_chain = chain_info_file2
        print("file 2 shorter")
    else:
        shortest_chain = chain_info_file1
        print("file 1 shorter")
    '''
#    print(file1_sequence)  # FOR TESTING
#    print(file2_sequence)  # FOR TESTING
    
    aligner = Bio.Align.PairwiseAligner()
    aligner.mode = 'global'
    aligner.open_gap_score = -10
    aligner.extend_gap_score = -0.5
    aligner.match_score = 1
    aligner.mismatch_score = -1 
    
    pairwise_alignments = aligner.align(file1_sequence,file2_sequence)

    pairwise_alignment = pairwise_alignments[0]  # Get the best alignment

#    pairwise_alignment = Alignment(pairwise_alignment)  # Convert to Biopython Alignment object

#    print(pairwise_alignment) # FOR TESTING
#    print(pairwise_alignment.indices) # FOR TESTING

    indices = np.array(pairwise_alignment.indices)
#    print(common_residues)

    index1_counter = 0
    index2_counter = 0

    '''print('file_index1 = '+str(len(file1_index)))
    print(file1_index)
    print('file_index2 = '+str(len(file2_index)))
    print(file2_index)

    print(indices)'''

    for n in range(indices.shape[1]):

        index1 = indices[0,n]
#        print("index1:", index1)  # FOR TESTING
        index2 = indices[1,n]
#        print("index2:", index2)  # FOR TESTING

        index1_str = file1_index[index1 + index1_counter] 
        index2_str = file2_index[index2 + index2_counter]

#        print(index1) # FOR TESTING
#        print(index1_str) # FOR TESTING

        '''if index1 == -1 and index2 != -1:
            index2_counter = index2_counter + 1
            print("index2 counter =" +str(index2_counter))  # FOR TESTING
            continue
        elif index1 != -1 and index2 == -1:
            index1_counter = index1_counter + 1
            print("index1 counter =" +str(index1_counter))  # FOR TESTING
            continue'''
        
        if index1 != -1 and index2 != -1:
            GDP_ca_positionX = chain_info_file1.get(index1_str)[1]
            GDP_ca_positionY = chain_info_file1.get(index1_str)[2]
            GDP_ca_positionZ = chain_info_file1.get(index1_str)[3]

            TAX_ca_positionX = chain_info_file2.get(index2_str)[1]
            TAX_ca_positionY = chain_info_file2.get(index2_str)[2]
            TAX_ca_positionZ = chain_info_file2.get(index2_str)[3]

            distance = math.sqrt((TAX_ca_positionX - GDP_ca_positionX)**2 + (TAX_ca_positionY - GDP_ca_positionY)**2 + (TAX_ca_positionZ - GDP_ca_positionZ)**2)

            if distance > max_distance:
                max_distance = distance
            
            distance_XY = math.sqrt((TAX_ca_positionX - GDP_ca_positionX)**2 + (TAX_ca_positionY - GDP_ca_positionY)**2)

            if distance_XY > max_distance_XY:
                max_distance_XY = distance_XY

            vector_componentX = (TAX_ca_positionX - GDP_ca_positionX)
            vector_componentY = (TAX_ca_positionY - GDP_ca_positionY)
            vector_componentZ = (TAX_ca_positionZ - GDP_ca_positionZ)

            if abs(vector_componentX) > max_vector_componentX:
                max_vector_componentX = abs(vector_componentX)
            if abs(vector_componentY) > max_vector_componentY:
                max_vector_componentY = abs(vector_componentY)
            if abs(vector_componentZ) > max_vector_componentZ:
                max_vector_componentZ = abs(vector_componentZ)

# print(max_distance)           
# print(max_distance_XY)

#Generate lists of distances between atoms.

for chain in chains:

    chain_info_file1 = {}
    file1_index = []
    file1_sequence = ''
    chain_info_file2 = {}
    file2_index = []
    file2_sequence = ''
    ca_distances = {}
    vector_components = {}
    colored_vectors = {} 

    pdbopen=open(file_1,'r')
    
    for line in pdbopen:
        if 'ATOM' in line:
            if line.split()[4] == chain:
                atom=line.split()[2]
                if atom=='CA':
                    residue_num = line.split()[5]   
                    residue = line.split()[3]
                    positionX = float(line.split()[6])
                    positionY = float(line.split()[7])
                    positionZ = float(line.split()[8])

                    chain_info_file1.update({residue_num:(residue,positionX,positionY,positionZ)})

                    file1_index.append(residue_num)

                    file1_sequence = file1_sequence + one_letter_codes.get(residue)
                    
#                    print(chain_info_file1)  # FOR TESTING
    pdbopen.close()

    pdbopen=open(file_2,'r')

    for line in pdbopen:
        if 'ATOM' in line:
            if line.split()[4] == chain:
                atom=line.split()[2]
                if atom=='CA':
                    residue_num = line.split()[5]   
                    residue = line.split()[3]
                    positionX = float(line.split()[6])
                    positionY = float(line.split()[7])
                    positionZ = float(line.split()[8])

                    chain_info_file2.update({residue_num:(residue,positionX,positionY,positionZ)})

                    file2_index.append(residue_num)

                    file2_sequence = file2_sequence + one_letter_codes.get(residue)                    
#                    print(chain_info_file1)  # FOR TESTING
    pdbopen.close()

#    print(file1_sequence)  # FOR TESTING
#    print(file2_sequence)  # FOR TESTING
    
    aligner = Bio.Align.PairwiseAligner()
    aligner.mode = 'global'
    aligner.open_gap_score = -10
    aligner.extend_gap_score = -0.5
    aligner.match_score = 1
    aligner.mismatch_score = -1 
    
    pairwise_alignments = aligner.align(file1_sequence,file2_sequence)

    pairwise_alignment = pairwise_alignments[0]  # Get the best alignment

#    print(pairwise_alignment) # FOR TESTING
#    print(pairwise_alignment.indices) # FOR TESTING

    indices = np.array(pairwise_alignment.indices)
#    print(common_residues)

    index1_counter = 0
    index2_counter = 0

    for n in range(indices.shape[1]):

        index1 = indices[0,n]
        index2 = indices[1,n]

        index1_str = file1_index[index1 + index1_counter] 
        index2_str = file2_index[index2 + index2_counter]

#        print(index1) # FOR TESTING
#        print(index1_str) # FOR TESTING

        if index1 != -1 and index2 != -1:
            GDP_ca_positionX = chain_info_file1.get(index1_str)[1]
            GDP_ca_positionY = chain_info_file1.get(index1_str)[2]
            GDP_ca_positionZ = chain_info_file1.get(index1_str)[3]

            TAX_ca_positionX = chain_info_file2.get(index2_str)[1]
            TAX_ca_positionY = chain_info_file2.get(index2_str)[2]
            TAX_ca_positionZ = chain_info_file2.get(index2_str)[3]

            distance = math.sqrt((TAX_ca_positionX - GDP_ca_positionX)**2 + (TAX_ca_positionY - GDP_ca_positionY)**2 + (TAX_ca_positionZ - GDP_ca_positionZ)**2)
            ca_distances.update({index1_str:distance})

            vector_componentX = (TAX_ca_positionX - GDP_ca_positionX)
            vector_componentY = (TAX_ca_positionY - GDP_ca_positionY)
            vector_componentZ = (TAX_ca_positionZ - GDP_ca_positionZ)

            vector_components.update({index1_str:(vector_componentX,vector_componentY,vector_componentZ)})

            R_value = abs(vector_componentX / max_vector_componentX)
            B_value = abs(vector_componentY / max_vector_componentY)
            G_value = abs(vector_componentZ / max_vector_componentZ)

            colored_vectors.update({index1_str:(R_value,B_value,G_value)})


#        print(R_value,B_value,G_value)  # FOR TESTING



    for n in ca_distances.keys():
        o1.write('\t#%s/%s:%s\t%s\n' %(file_1_number,chain,n,ca_distances.get(n)))
        o4.write('\t#%s/%s:%s\t%s\n' %(file_1_number,chain,n,math.sqrt((vector_components.get(n)[0])**2 + (vector_components.get(n)[1])**2)))
        o6.write('\t#%s/%s:%s\t%s\n' %(file_1_number,chain,n,vector_components.get(n)[2]))


#    print(ca_distances.keys())  # FOR TESTING

    index1_counter = 0
    index2_counter = 0

    for n in range(indices.shape[1]):

        index1 = indices[0,n]
        index2 = indices[1,n]

        index1_str = file1_index[index1 + index1_counter] 
        index2_str = file2_index[index2 + index2_counter]

#        print(index1) # FOR TESTING
#        print(index1_str) # FOR TESTING

        if index1 != -1 and index2 != -1:

            c1 = chain_info_file1.get(index1_str)
            c2 = chain_info_file2.get(index2_str)

            radius = ca_distances.get(index1_str)/(max_distance*2)  #scale arrow size to max distance

    #       o2.write('.color %s %s %s\n.arrow %s %s %s %s %s %s %s %s %s\n' %(cv[0],cv[1],cv[2],c1[1],c1[2],c1[3],c2[1],c2[2],c2[3],radius,radius*2,0.5))
            
            x = vector_components.get(index1_str)[0]
            y = vector_components.get(index1_str)[1]
            z = vector_components.get(index1_str)[2]

            x_ratio = abs(x)/(abs(x)+abs(y)+abs(z))
            y_ratio = abs(y)/(abs(x)+abs(y)+abs(z))
            z_ratio = abs(z)/(abs(x)+abs(y)+abs(z))

            x_ratio_scaled = x_ratio / max(x_ratio,y_ratio,z_ratio)
            y_ratio_scaled = y_ratio / max(x_ratio,y_ratio,z_ratio)
            z_ratio_scaled = z_ratio / max(x_ratio,y_ratio,z_ratio)

            o2.write('.color %s %s %s\n.arrow %s %s %s %s %s %s %s %s %s\n' %(x_ratio_scaled,y_ratio_scaled,z_ratio_scaled,c1[1],c1[2],c1[3],c2[1],c2[2],c2[3],radius,radius*2,0.5))

            o3.write("#"+file_1_number+"/"+chain+":"+index1_str+"@ca"+"\t"+"#"+file_2_number+"/"+chain+":"+index2_str+"@ca"+"\n")

            x_ratio_2D = abs(x)/(abs(y)+abs(x))
            y_ratio_2D = abs(y)/(abs(y)+abs(x))
            
            x_ratio_2D_scaled = x_ratio_2D / max(x_ratio_2D,y_ratio_2D)
            y_ratio_2D_scaled = y_ratio_2D / max(x_ratio_2D,y_ratio_2D)

    #        z_ratio_1D = abs(z)/(max_vector_componentZ)

            radius_Z = abs(z)/(max_distance*2)
            radius_XY = math.sqrt(x**2 + y**2)/(max_distance*2)

            if z < -1:
                o7.write('.color %s %s %s\n.arrow %s %s %s %s %s %s %s %s %s\n' %(0.75,0.39,0.31,c1[1],c1[2],c1[3],c1[1],c1[2],c2[3],radius_Z,radius_Z*2,0.5))

            elif z >= -1 and z <= 1:
                o7.write('.color %s %s %s\n.arrow %s %s %s %s %s %s %s %s %s\n' %(1,1,1,c1[1],c1[2],c1[3],c1[1],c1[2],c2[3],radius_Z,radius_Z*2,0.5))

            elif z > 1:
                o7.write('.color %s %s %s\n.arrow %s %s %s %s %s %s %s %s %s\n' %(0.126,0.66,0.83,c1[1],c1[2],c1[3],c1[1],c1[2],c2[3],radius_Z,radius_Z*2,0.5))
            
            o5.write('.color %s %s %s\n.arrow %s %s %s %s %s %s %s %s %s\n' %(x_ratio_2D_scaled,y_ratio_2D_scaled,0,c1[1],c1[2],c1[3],c2[1],c2[2],c1[3],radius_XY,radius_XY*2,0.5))
            

o1.close()
o2.close()
o3.close()
o4.close()
o5.close()
o6.close()
o7.close()
