import os
import sys
import re
import numpy as np
from Detect_CloseContact import Atomic_distance ,Get_Chain


def Get_Mol2Files (dir):
    directory = [files for folders,subfolders,files in os.walk(dir)] ## navegar pelo diretorio
    cleaned_directory = [file for file in directory[0] if file.endswith('.mol2')]
    
    return cleaned_directory

def Get_ProteinCoords (Complex_file, Aminoacids= ['ALA', 'ASX', 'CYS', 'ASP', 'GLU', 'PHE', 'GLY', 
                                                'LYS', 'LEU', 'MET', 'ASN', 'PRO', 'GLN', 'ARG', 
                                                'SER', 'THR', 'VAL', 'TRP', 'TYR', 'GLX','HIS','ILE']): 

    Complex = open('.\\'+Complex_file,'r').readlines()
    line = 0
    while Complex[line].startswith('@<TRIPOS>ATOM') == False:
        line+=1
    line+=1
    Protein_coords = []
    while Complex[line].startswith('@<TRIPOS>BOND') == False:
        Complex_line = list(filter(None, Complex[line].split(' ')))

        if Complex_line[-2][:3] in Aminoacids: ## protein
            #chain = (Get_Chain(Complex_file ,Complex_line[-3]+' '+Complex_line[-2]))
            #Complex_line.append(chain)    ## use quando o numero de ligantes x  30 seja maior do que o numero de atomos da proteina
            Protein_coords.append(Complex_line)
                               
        line+=1
        
        
    return Protein_coords

def Get_LigWaterCoords(Complex_file, Aminoacids= ['ALA', 'ASX', 'CYS', 'ASP', 'GLU', 'PHE', 'GLY', 
                                                'LYS', 'LEU', 'MET', 'ASN', 'PRO', 'GLN', 'ARG', 
                                                'SER', 'THR', 'VAL', 'TRP', 'TYR', 'GLX','HIS','ILE']): 
    Ligand_coords = []
    Water_coords = []
    Complex = open('.\\'+Complex_file,'r').readlines()
    line = 0
    while Complex[line].startswith('@<TRIPOS>ATOM') == False:
        line+=1
    line+=1
    while Complex[line].startswith('@<TRIPOS>BOND') == False:
        Complex_line = list(filter(None, Complex[line].split(' ')))

        if Complex_line[-2].startswith('HOH'): ## Water
            Water_coords.append(Complex_line)
        elif Complex_line[-2][:3] not in Aminoacids: ## Ligand
            Ligand_coords.append(Complex_line)
        line+=1    
    return Ligand_coords, Water_coords    


def Get_IntermolecularContacts (Complex_directory,Distance, Output_file):
    Complex_files = Get_Mol2Files(Complex_directory)
    print(Complex_files)
    Protein = Get_ProteinCoords(Complex_files[0])
    Contacts = []
    for Complex_file in Complex_files:
        Ligand_water = Get_LigWaterCoords(Complex_file)
        for Lig_atom in Ligand_water[0]:
            for prot_atom in Protein:
                distance = Atomic_distance(prot_atom[-7:-4], Lig_atom[-7:-4])
                if distance <= Distance:
                    
                    chain = (Get_Chain(Complex_file ,prot_atom[-3]+' '+prot_atom[-2])) ## use quando o numero de ligantes x  30 seja menor do que o numero de atomos da proteina
                    print(chain+':'+prot_atom[-2]+':'+prot_atom[1]+' - '+ Complex_file.split('.mol2')[0]+':'+Lig_atom[1]+'\t'+str(distance)+'\n')
                    Contacts.append(chain+':'+prot_atom[-2]+':'+prot_atom[1]+' - '+ Complex_file.split('.mol2')[0]+':'+Lig_atom[1]+'\t'+str(distance)+'\n')

            for water_atom in Ligand_water[1]:
                distance = Atomic_distance(water_atom[-7:-4], Lig_atom[-7:-4]) 
                if distance <= Distance:   
                    Contacts.append('Water:'+water_atom[-2]+':'+water_atom[1]+' - '+ Complex_file.split('.mol2')[0]+':'+Lig_atom[1]+'\t'+str(distance)+'\n')    



    contacts_file = open('.\\'+Output_file+'Contacts_file_'+str(Distance)+'_.txt','w')
    for line in Contacts:
        contacts_file.writelines(line)



if __name__ == '__main__':
 
    try:
        Complex_Directory = sys.argv[1]
        Distance = float(sys.argv[2])
        Output_file = sys.argv[3]
        
        Get_IntermolecularContacts (Complex_Directory, Distance, Output_file)
    except:
        print('USAGE ERROR: Use: python .\Detect_CloseContact_Water.py [.\Complexes_directory] [Contacts distance (Angstroms)] [Output name]')
        print('\nComplexes directory is the directory where you have your docked complexes from the molecular docking - GOLD\n\n')
        print('''\nThe script only works with files complexes in ".mol2" format\n''')
    
    
    
