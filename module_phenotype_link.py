print('Module-Phenotype relation Analysis Starting!')
print("Importing packages")
#import packages
import scipy.stats as stats
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import sys
from sklearn.decomposition import  PCA

#fetch input files
expression_f=sys.argv[1]
wgcna_f=sys.argv[2]
pheno_f=sys.argv[3]
output_f=sys.argv[4]

#define functions
def cor_dataframe(dat1,dat2):
    columns1=dat1.columns.tolist()
    columns2=dat2.columns.tolist()
    column_,index_,array,parray=[],[],np.zeros([len(columns1),len(columns2)]),np.zeros([len(columns1),len(columns2)])
    for i in range(len(columns1)):
        index_.append(columns1[i])
        for j in range(len(columns2)):
            if columns2[j] not in column_:
                column_.append(columns2[j])
            array[i,j]=stats.pearsonr(dat1[columns1[i]],dat2[columns2[j]])[0]
            parray[i,j]=stats.pearsonr(dat1[columns1[i]],dat2[columns2[j]])[1]
    result=pd.DataFrame(array)
    resultp=pd.DataFrame(parray)
    result.columns=column_
    resultp.columns=column_
    result.index=index_
    resultp.index=index_
    for c in result.columns.tolist():
        result[c]=result[c].fillna(0)
        resultp[c]=resultp[c].fillna(0)
    return (result,resultp)       
def split_dat(dataframe,kme):
    dat=pd.read_csv(dataframe,sep=',',index_col='ID')
    kme_dat=pd.read_csv(kme,sep=' ',header=None)
    modules=list(set([i for i in kme_dat[2]]))
    dic={}
    for module in modules:
        tmp=kme_dat.loc[kme_dat[2]==module]
        dic[module]=np.log2(dat.loc[dat.index.isin(tmp[0])]+1)
    return dic    
def pca_dataframe(dic_dataframe):
    dic={}
    for i in dic_dataframe.keys():
        
        pca=PCA(n_components=1)
        pca.fit(dic_dataframe[i].T)
        tmp=pca.transform(dic_dataframe[i].T)
        dic[i]=[i for i in tmp[:,0]]
    result=pd.DataFrame(dic)
    result.index=dic_dataframe[i].columns
    return result

#generate correlation matrix and p value matrix 
dic_dataframe=split_dat(expression_f,wgcna_f)
test2=pca_dataframe(dic_dataframe)
pheno=pd.read_csv(pheno_f,index_col='Sample',sep=',')
test3,test4=cor_dataframe(test2,pheno)
test3.to_csv(output_f+'.csv',sep='\t')
test4.to_csv(output_f+'_pvalue.csv',sep='\t')
#plot heatmap for correlation
plt.figure(figsize=(6.68,8),dpi=600)
sns.heatmap(test3,annot=True,cmap='OrRd',alpha=1,linewidths=0.1,fmt='.3g')
plt.savefig(output_f+'png',bbox_inches='tight')
plt.clf()