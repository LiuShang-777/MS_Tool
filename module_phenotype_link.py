print('Module-Phenotype relation Analysis Starting!')
print("Importing packages")
import scipy.stats as stats
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import sys
print("Reading input parameters")

expression_f=sys.argv[1]
wgcna_f=sys.argv[2]
pheno_f=sys.argv[3]
output_f=sys.argv[4]

print("Loading input data")
expression_dat=pd.read_csv(expression_f,sep=',',index_col='ID')
kme_dat=pd.read_csv(wgcna_f,sep=' ',header=None)
pheno_dat=pd.read_csv(pheno_f,sep=',',index_col='Sample')

print("Profiling transcription feature for modules")
dic_module={}
modules=list(set(kme_dat[2].tolist()))
for module in modules:
    module_genes=list(set(kme_dat.loc[kme_dat[2]==module][0].tolist()))
    tmp_dat=expression_dat.loc[expression_dat.index.isin(module_genes)]
    dic_module[module]=tmp_dat.median(axis=0)
    print("Generating heatmap for %s"%module)
    plt.figure(dpi=600)
    sns.heatmap(np.log2(tmp_dat+1),cmap='RdYlBu_r',yticklabels=False)
    plt.tight_layout()
    plt.savefig(output_f+'_%s.png'%module)
    plt.clf()        

print("Calculating Pearson correlationship between modules and phenotype")
pheno=pheno_dat.columns.tolist()
p_dat,r_dat=pd.DataFrame(),pd.DataFrame()

for module in modules:
    p_list,r_list=[],[]
    for phe in pheno:
        tmp_r,tmp_p=stats.pearsonr(dic_module[module],pheno_dat[phe])
        r_list.append(tmp_r)
        p_list.append(tmp_p)
    p_dat[module]=p_list
    r_dat[module]=r_list
p_dat.index,r_dat.index=pheno,pheno
p_dat=p_dat.T
r_dat=r_dat.T
r_dat.to_csv(output_f+'_correlation.csv',sep=',')
p_dat.to_csv(output_f+'_pvalue.csv',sep=',')
p_dat=-np.log10(p_dat)

print ("Plot heatmap for pearson correlationship ")
plt.figure(dpi=600)
sns.heatmap(p_dat,cmap='Blues',annot=True)
plt.tight_layout()
plt.savefig(output_f+'_pvalue.png')
plt.clf()        
plt.figure(dpi=600)
sns.heatmap(r_dat,cmap='Reds',annot=True)
plt.tight_layout()
plt.savefig(output_f+'_cvalue.png')
plt.clf()    

print("Analysis has been finished")
