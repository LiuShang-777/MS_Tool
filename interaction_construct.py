print("Start Interaction Network Construction")
print("Import Packages")
import pandas as pd
from scipy.stats import pearsonr
import sys
#Fetch input file 
target_gene=sys.argv[1]
gene_list_f=sys.argv[2]
dat_f=sys.argv[3]
outputf=sys.argv[4]

def get_network(gene,gene_list,express):
    source,target,weight=[],[],[]
    gene_array=express.loc[express.index==gene].iloc[0,:]
    for i in gene_list:
        tmp_array=express.loc[express.index==i].iloc[0,:]       
        tmp=pearsonr(gene_array,tmp_array)
        if tmp[1]<=0.05:
            source.append(gene)
            target.append(i)
            weight.append(tmp[0])
    result=pd.DataFrame()
    result['source']=source
    result['target']=target
    result['weight']=weight    
    return (result)

gene_list=[]
with open(gene_list_f,'r') as file:
    for line in file:
        line=line.strip()
        gene_list.append(line)
dat=pd.read_csv(dat_f,sep=',',index_col='ID')
result=get_network(target_gene,gene_list,dat)
result.to_csv(outputf,sep=',',index=False)