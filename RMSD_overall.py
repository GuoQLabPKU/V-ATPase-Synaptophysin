import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import itertools

from biopandas.pdb import PandasPdb
from Bio import pairwise2
from Bio.Seq import Seq
###parameter setting###
data_dir = 'Z:/4wenhong/v-ATPase-MS/RMSD_RatWithOtherSpecies/RMSD_state3_intact/'
##make the combination of file 
files = ['6wm3-mutated_delete-change-chain-coot-0.pdb', '7u4t_mutated_deleteH-change-chain-coot-0.pdb', 
          '7unf_mutated_delete-change-chain-coot-0.pdb',
          'state2_rat_deleteH.pdb', 'state2_mus_delH-change-chain.pdb']

# data_dir = 'Z:/4wenhong/v-ATPase-MS/RMSD_RatWithOtherSpecies/RMSD_state3_intact/'
# ##make the combination of file 
# files = ['KO_State3Ca_new.pdb', 'state3_autoSSpredictbyPymolCa.pdb' ]

##the main part of this script
three_letter = {'VAL': 'V',
 'ILE': 'I',
 'LEU': 'L',
 'GLU': 'E',
 'GLN': 'Q',
 'ASP': 'D',
 'ASN': 'N',
 'HIS': 'H',
 'TRP': 'W',
 'PHE': 'F',
 'TYR': 'Y',
 'ARG': 'R',
 'LYS': 'K',
 'SER': 'S',
 'THR': 'T',
 'MET': 'M',
 'ALA': 'A',
 'GLY': 'G',
 'PRO': 'P',
 'CYS': 'C'}
step1 = 1
if step1:
    for filePair in itertools.combinations(files, 2):
        file1, file2 = filePair
        #read the PDB file using      
        file1_PDB = PandasPdb().read_pdb(data_dir+file1)
        file2_PDB = PandasPdb().read_pdb(data_dir+file2)
        #extract the dataframe 
        pdb1_df = file1_PDB.df['ATOM']
        pdb2_df = file2_PDB.df['ATOM']
        ##get the chain name
        chainIds = np.unique(pdb1_df['chain_id'].values)
        #chainIds_v0 = ['d', 'a', 'p', 'b', 'g','h','i','j','k','l','m','n','o']
        ##should we focuso only on the CA chain?
        ifCA = 1
        if ifCA:
            file1_data = pdb1_df[pdb1_df['atom_name'] == 'CA']
            file2_data = pdb2_df[pdb2_df['atom_name'] == 'CA']
        ##record the RMSD
        distArray = [ ]
        ##cycle each chain :)
        for sgChain in chainIds:
            sgChain1 = file1_data[file1_data['chain_id'] == sgChain]
            sgChain2 = file2_data[file2_data['chain_id'] == sgChain]
            ##align two sequence 
            res_chain1 = [three_letter[i] for i in sgChain1['residue_name'].values]
            res_chain2 = [three_letter[i] for i in sgChain2['residue_name'].values]
            res_seq1 = Seq(''.join(res_chain1))
            res_seq2 = Seq(''.join(res_chain2))
            alignments = pairwise2.align.globalxx(res_seq1, res_seq2)
            ali = alignments[0]
            ali_seq1 = ali[0]
            ali_seq2 = ali[1]
            ##store the coordinates of each resiude
            coord1 = np.zeros((len(ali_seq1), 3))
            coord2 = np.zeros((len(ali_seq2), 3))
            
            ##record the coord1
            label_seq1 = np.zeros(sgChain1.shape[0])
            for m, sgRes in enumerate(ali_seq1):
                if sgRes == '-':
                    continue
                else:
                    for n, sgRes_ori in enumerate(res_chain1):
                        if (label_seq1[n] == 0) & (sgRes_ori == sgRes):
                            sgPosXYZ1 = sgChain1[['x_coord','y_coord','z_coord']].values[n]
                            sgPosXYZ1 = np.array([float(j) for j in sgPosXYZ1])               
                            coord1[m,:] = sgPosXYZ1
                            label_seq1[n] = 1
                            break
                        else:
                            continue
                                   
            ##record the coord2                
            label_seq2 = np.zeros(sgChain2.shape[0])
            for p, sgRes in enumerate(ali_seq2):
                if sgRes == '-':
                    continue
                else:
                    for q, sgRes_ori in enumerate(res_chain2):
                        if (label_seq2[q] == 0) & (sgRes_ori == sgRes):
                            sgPosXYZ2 = sgChain2[['x_coord','y_coord','z_coord']].values[q]
                            sgPosXYZ2 = np.array([float(j) for j in sgPosXYZ2])               
                            coord2[p,:] = sgPosXYZ2
                            label_seq2[q] = 1
                            break
                        else:
                            continue                        
         
            ##double check if there exist gap, else skip the coordinate
            if ('-' in ali_seq1) | ('-' in ali_seq2):
                gapCoord1 = [i for i,j in enumerate(ali_seq1) if j=='-']
                gapCoord2 = [i for i,j in enumerate(ali_seq2) if j=='-']
                #combine the coordinates 
                gapCoord1.extend(gapCoord2)
                for sgGap in gapCoord1:
                    coord1[sgGap,:] = np.zeros(3)
                    coord2[sgGap,:] = np.zeros(3)
            ##calculate the distance of this chain:
            distPairs = np.linalg.norm(coord2-coord1, axis=1)
            distArray.extend(distPairs)
                
        ##calculate the RMSD
        RMSD = np.mean(distArray)
        print(filePair, RMSD)
                
                
##in step2, we just make the plot like heatmap to show that 
##create a dataframe to vis the heatmap
names = ['this work','7unf','6vqg','6wm3','7u4t']       
values_df = pd.DataFrame({names[0]:[12.498, 11.39, 14.41, 4.072, 0],
                          names[1]:[12.484,11.65,14.37,0,1.69],
                          names[2]:[3.67,4.35,0,4.85,4.72],
                          names[3]:[2.3,0,1.45,4.36,4.30],
                          names[4]:[0,1.24,1.60,4.89,4.99]})
names.reverse()
values_df['names'] =  names
##reset the index 
values_df.set_index('names', drop=True, inplace = True)
step2 = 0
if step2:
    ##create mask 
    mask = np.zeros((5,5))
    mask[np.triu_indices_from(mask)] = True
    row,col = np.diag_indices_from(mask)
    mask[row,col] = 0
    sns.heatmap(values_df,vmin = 0, vmax = 14.651191125571325,cmap="RdBu_r",  linewidths=0.3, 
                square=True, cbar_kws = {'shrink':0.8})
    plt.yticks(rotation = 0, fontsize = 15)
    plt.xticks(fontsize = 15)
    plt.xlabel(None)



#in step3, calculate the RMSD for  each subunit of eahc PDB pair



    
            
            
            
