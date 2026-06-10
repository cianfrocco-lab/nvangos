#!/usr/bin/env python

#Import necessary modules

import argparse
import math
import os
import numpy as np
import Bio.Align
from Bio.Align import Alignment

one_letter_codes = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
     'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N', 
     'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W', 
     'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}


def parse_ca_atom_line(line, target_chain):
    """Return CA atom info for standard PDB or packed labels like METB1."""
    if not line.startswith("ATOM"):
        return None

    fields = line.split()
    if len(fields) < 8 or fields[2] != "CA":
        return None

    residue_field = fields[3]

    if len(residue_field) > 3 and residue_field[:3] in one_letter_codes:
        residue = residue_field[:3]
        chain = residue_field[3:]
        residue_num = fields[4]
        coord_start = 5
    else:
        if len(fields) < 9:
            return None
        residue = residue_field
        chain = fields[4]
        residue_num = fields[5]
        coord_start = 6

    if chain != target_chain:
        return None

    try:
        return (
            residue_num,
            residue,
            float(fields[coord_start]),
            float(fields[coord_start + 1]),
            float(fields[coord_start + 2]),
        )
    except (IndexError, ValueError):
        return None


def load_chain_ca_atoms(pdb_file, chain):
    chain_info = {}
    residue_index = []
    sequence = ''

    with open(pdb_file, 'r') as pdbopen:
        for line in pdbopen:
            parsed = parse_ca_atom_line(line, chain)
            if parsed is None:
                continue

            residue_num, residue, positionX, positionY, positionZ = parsed
            one_letter = one_letter_codes.get(residue)
            if one_letter is None:
                continue

            chain_info.update({residue_num: (residue, positionX, positionY, positionZ)})
            residue_index.append(residue_num)
            sequence = sequence + one_letter

    return chain_info, residue_index, sequence


def get_ca_chain_ids(pdb_file):
    chains = []
    seen = set()

    with open(pdb_file, 'r') as pdbopen:
        for line in pdbopen:
            if not line.startswith("ATOM"):
                continue

            fields = line.split()
            if len(fields) < 8 or fields[2] != "CA":
                continue

            residue_field = fields[3]
            if len(residue_field) > 3 and residue_field[:3] in one_letter_codes:
                chain = residue_field[3:]
            elif len(fields) >= 9:
                chain = fields[4]
            else:
                continue

            if chain not in seen:
                seen.add(chain)
                chains.append(chain)

    return chains


def parse_args():
    parser = argparse.ArgumentParser(
        description="Generate ChimeraX .defattr, .bild, and .pb files from per-residue CA vectors between two aligned PDB models."
    )
    parser.add_argument(
        "file_1",
        help="First aligned PDB model.",
    )
    parser.add_argument(
        "file_2",
        help="Second aligned PDB model.",
    )
    parser.add_argument(
        "--chains",
        nargs="+",
        default=None,
        help="Chain IDs to analyze, in order. If omitted, all matching CA chain IDs present in both models are used.",
    )
    parser.add_argument("--file-1-number", default="1", help="ChimeraX model number for file 1. Default: %(default)s")
    parser.add_argument("--file-2-number", default="2", help="ChimeraX model number for file 2. Default: %(default)s")
    parser.add_argument("--color", default="blue", help="Pseudobond color written to the .pb header. Default: %(default)s")
    parser.add_argument("--radius", default="0.3", help="Pseudobond radius written to the .pb header. Default: %(default)s")
    parser.add_argument("--dashes", default="1", help="Pseudobond dashes value written to the .pb header. Default: %(default)s")
    parser.add_argument("--out-ca", default=None, help="Output CA distance .defattr file. Default: ca_distances.defattr beside file_1")
    parser.add_argument("--out-bild-xyz", default=None, help="Output XYZ .bild file. Default: colored_vectors_XYZ.bild beside file_1")
    parser.add_argument("--out-pb-xyz", default=None, help="Output XYZ .pb file. Default: colored_vectors_XYZ.pb beside file_1")
    parser.add_argument("--out-ca-xy", default=None, help="Output XY distance .defattr file. Default: ca_distances_XY_only.defattr beside file_1")
    parser.add_argument("--out-bild-xy", default=None, help="Output XY .bild file. Default: colored_vectors_XY_only.bild beside file_1")
    parser.add_argument("--out-ca-z", default=None, help="Output Z distance .defattr file. Default: ca_distances_Z_only.defattr beside file_1")
    parser.add_argument("--out-bild-z", default=None, help="Output Z .bild file. Default: colored_vectors_Z_only.bild beside file_1")
    return parser.parse_args()


def resolve_chains(file_1, file_2, requested_chains=None):
    file_1_chains = get_ca_chain_ids(file_1)
    file_2_chains = get_ca_chain_ids(file_2)
    file_1_set = set(file_1_chains)
    file_2_set = set(file_2_chains)

    if file_1_set != file_2_set:
        only_file_1 = sorted(file_1_set - file_2_set)
        only_file_2 = sorted(file_2_set - file_1_set)
        raise ValueError(
            "Input models do not have the same CA chain IDs.\n"
            "%s chains: %s\n"
            "%s chains: %s\n"
            "Only in %s: %s\n"
            "Only in %s: %s"
            % (file_1, file_1_chains, file_2, file_2_chains, file_1, only_file_1, file_2, only_file_2)
        )

    if requested_chains is None:
        return sorted(file_1_chains)

    missing_file_1 = [chain for chain in requested_chains if chain not in file_1_set]
    missing_file_2 = [chain for chain in requested_chains if chain not in file_2_set]

    if missing_file_1 or missing_file_2:
        raise ValueError(
            "Requested chains are missing from one or both models.\n"
            "Missing from %s: %s\n"
            "Missing from %s: %s"
            % (file_1, missing_file_1, file_2, missing_file_2)
        )

    return requested_chains


def output_path_for_file_1(file_1, output_path, default_name):
    file_1_dir = os.path.dirname(os.path.abspath(file_1))
    output_name = output_path if output_path is not None else default_name

    if os.path.isabs(output_name):
        return output_name

    return os.path.join(file_1_dir, output_name)

args = parse_args()

file_1 = args.file_1
file_1_number = args.file_1_number
file_2 = args.file_2
file_2_number = args.file_2_number

color = args.color
radius = args.radius
dashes = args.dashes

chains = resolve_chains(file_1, file_2, args.chains)
print("Analyzing chains: %s" % ", ".join(chains))

outfile1 = output_path_for_file_1(file_1, args.out_ca, "ca_distances.defattr")
outfile2 = output_path_for_file_1(file_1, args.out_bild_xyz, "colored_vectors_XYZ.bild")
outfile3 = output_path_for_file_1(file_1, args.out_pb_xyz, "colored_vectors_XYZ.pb")
outfile4 = output_path_for_file_1(file_1, args.out_ca_xy, "ca_distances_XY_only.defattr")
outfile5 = output_path_for_file_1(file_1, args.out_bild_xy, "colored_vectors_XY_only.bild")
outfile6 = output_path_for_file_1(file_1, args.out_ca_z, "ca_distances_Z_only.defattr")
outfile7 = output_path_for_file_1(file_1, args.out_bild_z, "colored_vectors_Z_only.bild")


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

    chain_info_file1, file1_index, file1_sequence = load_chain_ca_atoms(file_1, chain)
    chain_info_file2, file2_index, file2_sequence = load_chain_ca_atoms(file_2, chain)

    if len(file1_sequence) == 0 or len(file2_sequence) == 0:
        raise ValueError(
            "No CA sequence found for chain %s. %s has %s CA residues; %s has %s CA residues."
            % (chain, file_1, len(file1_sequence), file_2, len(file2_sequence))
        )

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

    chain_info_file1, file1_index, file1_sequence = load_chain_ca_atoms(file_1, chain)
    chain_info_file2, file2_index, file2_sequence = load_chain_ca_atoms(file_2, chain)

    if len(file1_sequence) == 0 or len(file2_sequence) == 0:
        raise ValueError(
            "No CA sequence found for chain %s. %s has %s CA residues; %s has %s CA residues."
            % (chain, file_1, len(file1_sequence), file_2, len(file2_sequence))
        )

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
