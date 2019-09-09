from __future__ import print_function
import codecs
import scipy.io as sio
from molml.features import CoulombMatrix 
from molml.features import BagOfBonds 

from rdkit import Chem
import numpy as np
import scipy.io as sio
from rdkit.Chem import Descriptors
from rdkit.Chem import AllChem
import codecs
from molml.features import CoulombMatrix
import time
from progressbar import *


filetype=input('1:xyz / 2:sdf')

if filetype=='1':
    molsum=input("what's the max molecule number?")
    molsum=int(molsum)
    filepath=input("what's xyz file name?")
    fd=codecs.open(filepath+'.xyz')   
    text=fd.readlines()
    index=0
    xyz_dict={}
    mol_no=molsum
    CM_dict=[]
    for i in range(0,molsum):    
        atom_nums=int(text[index])
        a=text[index+2:index+atom_nums+2]
        name='mol'+ str(mol_no)
        xyz_dict[name]=a
        index=index+atom_nums+2
        mol_no=mol_no-1
        
        atom_list=[]
        axis=[]
        for j in range(0,len(a)):
            atom=a[j].split()        
            atom_list.append(atom[0])
            axis.append([float(atom[1]),float(atom[2]),float(atom[3])])
            
        feat=CoulombMatrix()
        mole=(atom_list,axis)
        feat.fit([mole])
        CM=feat.transform([mole])[0]
        t=CM.reshape((atom_nums,atom_nums)).tolist()
        CM_dict.append(t)
    fd.close()
    sio.savemat('CM.mat',{'CM':CM_dict})
    print("successfully transferred!")
else:
    #get CM
    def smi2cm(m):
        m1 = Chem.MolFromSmiles(m)
        m = Chem.AddHs(m1)
        AllChem.EmbedMolecule(m,AllChem.ETKDG())
        n_atoms = m.GetNumAtoms()
        m1=Chem.MolToMolBlock(m)
        m1=m1.split()
        axis=[]
        atom_list=[]
        for i in range(0,n_atoms):
            axis.append([float(m1[13+16*i]),float(m1[14+i*16]),float(m1[15+16*i])])
            atom_list.append(m1[16+16*i])
        feat=CoulombMatrix()
        mole=(atom_list,axis)
        feat.fit([mole])
        t=feat.transform([mole])[0]
        CM=t.reshape((n_atoms,n_atoms)).tolist()
        
        return CM
    #print(len(smi2cm('c1ccccc1')))
    
    #supply = Chem.SDMolSupplier('F:\#project_ML\dataset\saturated_hydrocarbons.sdf')
    filepath=input("what's sdf file name?")+'.sdf'
   # filepath='saturated_hydrocarbons.sdf'
    supply = Chem.SDMolSupplier(filepath)
    #supply = Chem.SDMolSupplier('tran_data/sdf_db.sdf')
    molsum=len(supply)
    CM_dict=[];
    n=1
    total =int(input("how many molecules you want to transfer?(maximun:%s)" %(molsum)))
    def dosomework():
        time.sleep(0.01)
    pbar = ProgressBar().start()
    
    #while 10000<n<=20000:
    for mol in supply:
        if total>n:
#         print(mol.GetNumAtoms())
            Chem.Kekulize(mol)
            s=Chem.MolToSmiles(mol,kekuleSmiles=True)
            c=smi2cm(s)
            CM_dict.append(c)
            n=n+1
            pbar.update(int((n / (total +1)) * 100)) 
            dosomework()
        else:
            break
    pbar.finish()

#CMall=list(c)
sio.savemat('sdf2cm.mat',{'CM':CM_dict})
print("successfully transferred!")