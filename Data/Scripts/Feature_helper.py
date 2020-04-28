###################################################################################
# Contain functions and classes to help to extract secondary structure information
###################################################################################

import sys
import numpy as np
from Bio.PDB.DSSP import dssp_dict_from_pdb_file
from group import smoothForSS

def remove(string,char):
    '''
    Remove the space in a string.
    '''
    string_char = [i for i in string.split(char) if i != '']
    return string_char[0]

def read_pdb(pdb_file):
    '''
    Extract the residue information inside a pdb file.
    '''
    protein_dict = {}
    with open(pdb_file,'r') as p_file:
        lines = p_file.readlines()
        for line in lines:
            if line[0:4] == 'ATOM' and line[26] == ' ':
                atom = remove(line[13:17],' ')
                resi = line[17:20]
                chain = line[21]
                index = int(remove(line[22:26],' '))
                x = float(remove(line[30:38],' '))
                y = float(remove(line[38:46],' '))
                z = float(remove(line[46:54],' '))
        ############ Judge whether a new chain begins. ########################      
                if not chain in protein_dict.keys():
                    protein_dict[chain] = {}
        ############ Save the sequence infomation. ######################## 
                if not index in protein_dict[chain].keys():
                    protein_dict[chain][index] = {'resi':resi}
                elif resi != protein_dict[chain][index]['resi']:
                    print('PDB read error! The residue kind of resi %d is not consistent!')%index
                    return 0
                protein_dict[chain][index][atom] = [x,y,z]  
    return protein_dict


class SS_extraction:
   
    def __init__(self,pdb_file): 
        self.AA_dict = {'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C','GLN':'Q','GLU':'E','GLY':'G','HIS':'H','ILE':'I','LEU':'L','LYS':'K','MET':'M','PHE':'F','PRO':'P','SER':'S','THR':'T','TRP':'W','TYR':'Y','VAL':'V'}
        self.ss_dict_8_3 = {'H':'H','G':'H','I':'H','E':'E','B':'E','S':'C','T':'C','-':'C','C':'C'} # 8-classes to 3-classes and 3 to 3
        self.ss_edge_label_dict = {('H','H'):'000',('H','E'):'001',('H','C'):'010',
                              ('E','E'):'011',('E','H'):'001',('E','C'):'100',
                              ('C','C'):'101',('C','H'):'010',('C','E'):'100'}

        self.pdb_file = pdb_file
        self.protein_dict = read_pdb(pdb_file)
        self.Seq_dict = {}
        self.SS_dict_8 = {}
        self.SS_dict_3 = {}
        dssp_dict = dssp_dict_from_pdb_file(pdb_file)[0]

        for chain in self.protein_dict.keys():
            Complete_Seq = ''
            Complete_SS_8 = ''
            Complete_SS_3 = ''
            index_min = min(self.protein_dict[chain].keys())
            index_max = max(self.protein_dict[chain].keys())
            for index in range(index_min,index_max + 1):
                if index in self.protein_dict[chain].keys():
                    self.protein_dict[chain][index]['SeconStru'] = dssp_dict[chain, (' ', index, ' ')][1]
                    if self.AA_dict[self.protein_dict[chain][index]['resi']] != dssp_dict[chain, (' ', index, ' ')][0]:
                        print('Residue Error! %s and %s do not match!'%(self.protein_dict[chain][index]['resi'],dssp_dict[chain, (' ', index, ' ')][0]))
                    else:
                        self.protein_dict[chain][index]['AminoAci'] = dssp_dict[chain, (' ', index, ' ')][0]
                    Complete_Seq += dssp_dict[chain, (' ', index, ' ')][0]
                    Complete_SS_8 += dssp_dict[chain, (' ', index, ' ')][1]
                    Complete_SS_3 += self.ss_dict_8_3[dssp_dict[chain, (' ', index, ' ')][1]]
                else:
                    Complete_Seq += 'x'
                    Complete_SS_8 += 'x'
                    Complete_SS_3 += 'x'
            self.Seq_dict[chain] = Complete_Seq
            self.SS_dict_8[chain] = Complete_SS_8
            self.SS_dict_3[chain] = Complete_SS_3

