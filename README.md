# MS_Tool
Module-related subsequent analysis tool for WGCNA

This tool is totally developed by Python script which means it could run directly by Python in command line mode.

Prerequisite packages:
python (>=3.7)

numpy (>=1.17.4)

pandas (>=0.23.4)

seaborn (>=0.9.0)

matplotlib (>=2.2.3)

sklearn (>=0.19.2)

scipy (>=1.1.0)

The function of modulate-phenotype-link could be used like this:

#modulate-phenotype-link.py

#1 python modulate-phenotype-link.py expression.csv module.txt phenotype.csv output_prefix

Input file:

expression.csv (seperated by tab, the first column is "ID")

ID  SAMPLE1 SAMPLE2...

gene1 12  20

gene2 30  40

module.txt (seperated by space)

gene1 0.9 yellow

gene2 0.8 blue

phenotype.txt (seperated by comma, the index of the column is "Sample")

Sample,tissue1,tissue2...

sample1,pheno1,pheno2

Output file:

output_prefix.csv (Pearson correlationship matrix)

output_prefix_pvalue.csv (P value of pearson correlation)

output_prefix.png

#interaction-construct.py

#2 python interaction-construct.py gene gene_list_file expression.csv

Input file:

gene: gene name

gene_list_file:

gene1
gene2

expression.csv (The same format as the expression.csv in modulate-phenotype-link.py function)
