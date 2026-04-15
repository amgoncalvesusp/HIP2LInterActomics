import sys
import re
import numpy as np


def Get_ProtAtoms (Prot_file):

    Prot_data = open (Prot_file, 'r').readlines ()

    line = 0    
    while Prot_data[line].startswith('@<TRIPOS>ATOM') == False:
        line+=1
    line+=1   
    Prot_Atoms = []
    while Prot_data[line].startswith('@<TRIPOS>BOND') == False:
        Inf_ProtAtom = re.findall('[0-9][0-9]*.[A-Z][A-Za-z]*[0-9]*', Prot_data[line])
        Amcd_ProtAtom = Inf_ProtAtom[-1]
        ID_ProtAtom = Inf_ProtAtom[0].split(' ')[1]
        Coords_ProtAtom = re.findall('\-*[0-9][0-9]*\.\d{4}', Prot_data[line])[:-1]
        Prot_Atoms.append ([ID_ProtAtom, Amcd_ProtAtom, Coords_ProtAtom])
        line+=1

    return Prot_Atoms

def Get_Chain (Prot_file , Amcd_ProtAtom):
    Prot_data = open (Prot_file, 'r').readlines ()

    line = 0
    while Prot_data[line].startswith('@<TRIPOS>SUBSTRUCTURE') == False:
        line+=1
    line+=1    
      
    while len(re.findall(Amcd_ProtAtom, Prot_data[line])) == 0:
        line+=1
    
    return re.findall('[0-9][0-9]*.[A-Z].....',Prot_data[line])[-1].split(' ')[1]     

def Get_LigsAtom (Ligs_file):
    Ligs_data = open (Ligs_file, 'r').readlines ()

    No_molecules = np.unique(Ligs_data, return_counts = True)[1][np.where(np.unique(Ligs_data, return_counts = True)[0] =='@<TRIPOS>MOLECULE\n')][0]
    line = 0
    Mols_name = []
    Mols_Score = []
    Mols_ID = []
    Mols_Coords = []
    
    for molecule in range(No_molecules):
        
        while Ligs_data[line].startswith('#       Name:') == False:
            line+=1
        Name = re.findall('[0-9][0-9]*', Ligs_data[line])[0]
        line+=1
        
        if Name in Mols_name:
            idx =Mols_name.index(Name)
            
            line_score_comparativo = line
            while Ligs_data[line_score_comparativo].startswith('> <Gold.Score>') == False:
                line_score_comparativo+=1
              
            if Ligs_data[line_score_comparativo+2].split(' ')[5] >= Mols_Score[idx]:
                Mols_Score[idx] = Ligs_data[line_score_comparativo+2].split(' ')[5]
                
                line_replace = line
                while Ligs_data[line_replace].startswith('@<TRIPOS>MOLECULE') == False:
                    line_replace+=1
                Mols_ID[idx] = Ligs_data[line_replace+1].split('|')[0]
                
                while Ligs_data[line_replace].startswith('@<TRIPOS>ATOM') == False:
                    line_replace+=1
                Coord_list= []
                line_replace+=1    
                while Ligs_data[line_replace].startswith('@<TRIPOS>BOND') == False:
                    Atom_ID = re.findall('[A-Z][A-Za-z]*[0-9]*', Ligs_data[line_replace])[0]
                    if Atom_ID != 'LP':  
                        coords = re.findall('\-*[0-9][0-9]*\.\d{4}', Ligs_data[line_replace])[:-1]
                        Coord_list.append([Atom_ID, coords[0], coords[1], coords[2]])
                    line_replace+=1    
                
                Mols_Coords[idx] = Coord_list
                

        else:    
            Mols_name.append(Name)  
            line_ID = line

            while Ligs_data[line_ID].startswith('@<TRIPOS>MOLECULE') == False:
                line_ID+=1
            Mols_ID.append(Ligs_data[line_ID+1].split('|')[0])

            line_coord = line_ID+1
            while Ligs_data[line_coord].startswith('@<TRIPOS>ATOM') == False:
                line_coord+=1
            line_coord+=1
            Coord_list= []    
            while Ligs_data[line_coord].startswith('@<TRIPOS>BOND') == False:
                Atom_ID = re.findall('[A-Z][A-Za-z]*[0-9]*', Ligs_data[line_coord])[0]
                
                if Atom_ID != 'LP':  
                    coords = re.findall('\-*[0-9][0-9]*\.\d{4}', Ligs_data[line_coord])[:-1]
                    Coord_list.append([Atom_ID, coords[0],coords[1],coords[2]])
                line_coord+=1    
            Mols_Coords.append(Coord_list)
            line_score = line_coord
            while Ligs_data[line_score].startswith('> <Gold.Score>') == False:
                line_score+=1
            
            Mols_Score.append(Ligs_data[line_score+2].split(' ')[5])
    
    return Mols_ID, Mols_Score, Mols_Coords

    
def Atomic_distance(Coords_ProtAtom,Coords_LigAtom):

    """calculate the interatomic distance in Angstroms""" 
    Coords_ProtAtom = np.array(Coords_ProtAtom).astype(float)
    Coords_LigAtom = np.array(Coords_LigAtom).astype(float)
    
    return ((Coords_ProtAtom[0]-Coords_LigAtom[0])**2 + (Coords_ProtAtom[1]-Coords_LigAtom[1])**2  + (Coords_ProtAtom[2]-Coords_LigAtom[2])**2)**0.5

def Get_IntermolecularContacts (Prot_file, Ligs_file, Distance, Output_file):
        

    Protein = Get_ProtAtoms(Prot_file)
    
    Molecule_ = Get_LigsAtom (Ligs_file)
    Molecule_names = Molecule_[0]
    Molecule_scores = Molecule_[1]
    Molecule_coords = Molecule_[2]
    Contacts = []
    for prot_atom in Protein:
        for indx, molecule in enumerate(Molecule_coords):
            for molecule_atom in molecule:
                distance = Atomic_distance(prot_atom[2], molecule_atom[1:])
                if distance <= Distance:
                    Contacts.append(Get_Chain(Prot_file , prot_atom[1])+':'+prot_atom[1].split(' ')[1]+':'+prot_atom[0]+' - '+Molecule_names[indx]+':'+molecule_atom[0]+'\t'+str(distance)+'\n')


    contacts_file = open('.\\'+Output_file+'Contacts_file_'+str(Distance)+'_.txt','w')
    for line in Contacts:
        contacts_file.writelines(line)




if __name__=='__main__':
    
    try:
        Prot_file = sys.argv[1]
        Ligs_file = sys.argv[2]
        Distance = float(sys.argv[3])
        Output_file = sys.argv[4]
        Get_IntermolecularContacts (Prot_file, Ligs_file, Distance, Output_file)
    except:
        print('USAGE ERROR: Use: python .\Detect_CloseContact.py [.\Protein_file.mol2] [.\Molecule_solution.mol2] [Contacts distance (A)] [Output name]')
        print('\nProtein file is the used for the molecular docking - GOLD in ".mol2" format\n\n')
        print('''\nMolecule_solution file is the result file with the poses from the molecular docking - GOLD,\n
        the Script works with files with Multiple-molecule in ".mol2" format\n''')
       
    #Prot_file = sys.argv[1]
    #Ligs_file = sys.argv[2]
    #Distance = float(sys.argv[3])
    #Output_file = sys.argv[4]
    #Get_IntermolecularContacts (Prot_file, Ligs_file, Distance, Output_file)
    
    