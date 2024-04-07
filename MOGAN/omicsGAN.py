import os
import sys
import pandas as pd 
import numpy as np 
import argparse
from omics1 import omics1
from omics2 import omics2
from omics3 import omics3

total_update = int(sys.argv[1])
omics1_name = sys.argv[2]
omics2_name = sys.argv[3]
omics3_name = sys.argv[4]
adj_file1 = sys.argv[5]
adj_file2 = sys.argv[6]
label = sys.argv[7]

omics1_result = []
omics2_result = []
new_omics1_result = []
omics3_result = []
flag = 0
for i in range(1, total_update+1):
	omics1_result.append(omics1(i,omics1_name,omics2_name,adj_file1,label,flag))
	omics2_result.append(omics2(i,omics1_name,omics2_name,adj_file1,label,flag))

omics1_result = np.array(omics1_result)
omics2_result = np.array(omics2_result)

keep_mRNA = (np.argsort(np.mean(omics1_result,axis=1))[::-1][0])
keep_protein = (np.argsort(np.mean(omics2_result,axis=1))[::-1][0])

print('Best prediction for omics1: ',omics1_result[keep_mRNA])
print('Best update for omics1: ',keep_mRNA+1)
print('Best prediction for omics2: ',omics2_result[keep_protein])
print('Best update for omics2: ',keep_protein+1)

for i in range(1, total_update+1):
	if i!=keep_mRNA+1:
		os.remove("omics1_"+str(i)+".csv") 
	os.remove("omics2_"+str(i)+".csv")

flag = 1
new_omics1_file = 'omics1_' + str((keep_mRNA+1)) + '.csv'
data = pd.read_csv(new_omics1_file,index_col=0)
data = data.transpose()
mRNA_protein_fusion = 'mRNA_Protein.csv'
data.to_csv(mRNA_protein_fusion)
os.remove('omics1_' + str((keep_mRNA+1)) + '.csv')


for i in range(1, total_update+1):
	new_omics1_result.append(omics1(i,mRNA_protein_fusion,omics3_name,adj_file2,label,flag))
	omics3_result.append(omics3(i,mRNA_protein_fusion,omics3_name,adj_file2,label))

new_omics1_result = np.array(new_omics1_result)
omics3_result = np.array(omics3_result)

keep_mRNA = (np.argsort(np.mean(new_omics1_result,axis=1))[::-1][0])
keep_meth = (np.argsort(np.mean(omics3_result,axis=1))[::-1][0])

print('Best prediction for generate_omics1: ',new_omics1_result[keep_mRNA])
print('Best update for generate_omics1: ',keep_mRNA+1)
print('Best prediction for omics3: ',omics3_result[keep_meth])
print('Best update for omics3: ',keep_meth+1)

for i in range(1, total_update+1):
	if i!=keep_mRNA+1:
		os.remove("omics1_"+str(i)+".csv") 

	if i!=keep_meth+1:
		os.remove("omics3_"+str(i)+".csv")