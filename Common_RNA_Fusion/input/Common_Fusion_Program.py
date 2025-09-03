import pandas as pd
import glob
import os
from functools import reduce
#import re
#import numpy as np

###################################### Start PART I ################################
#set the working directory
path = os.getcwd()
files = glob.glob(os.path.join(path,"*.preliminary"))
print(files)

#Modifying each gene fusion preliminary file  #Modified Edited by Vyomesh 16-Feb-22
def concat_arrange(fp):
    df = pd.read_csv(fp, sep='\t')
    df["Score"] = df["Score"].astype(float).round(3)
    df["Score"] = df["Score"].astype(str)
    df['Read_Support_Count']=df['ReadNames'].str.split(';').apply((len)) 
    #if next value is blank
    df['Read_Support_Count']=df['Read_Support_Count']-1
    df['Combined'] = df['Score'].map(str) + '|(' + df['LeftBreakpoint'].map(str) + ' ' + df['RightBreakpoint'].map(str) +')|'+ df['Read_Support_Count'].map(str)
    #row corresponding to same gene merged on first row as coloumns
    #df1_intrem=df.groupby('#FusionGene',as_index=False).agg({'Combined':list,'#FusionGene':'first'})
    df1_intrem = df.groupby('#FusionGene')['Combined'].apply(list).reset_index()
    df1_intrem=df1_intrem.join(pd.DataFrame(df1_intrem.pop('Combined').tolist()).rename(columns=lambda x:f"Combined{x+1}"))
    df_col_list=list(df1_intrem.columns)
    df_col_list.pop(0)
    #Combining the info of multiple records and add it into single coloumn
    df1_intrem['combined'] = df1_intrem[df_col_list].apply(lambda row:'|#|'.join(row.values.astype(str)), axis=1)
    #Renaming 2nd coloum to Sample ID 
    sample_col_name=os.path.basename(fp).split('.')[0]
    a='combined'
    df1_intrem.rename(columns = {a: sample_col_name}, inplace = True) 
    df_final = df1_intrem[['#FusionGene',sample_col_name]]
    #replace unwanted substring from values in 2nd coloumn
    df_final=df_final.replace(to_replace ='\|\#\|None', value = '', regex = True)
    df_final.to_csv(r'./rna_final_files/'+sample_col_name+'_merged_FUS.csv', index = False)
  
#Calling the above function on all the files
for f in files:
    concat_arrange(f)
###################################### END PART I ################################

###################################### Start PART II ################################                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     os.chdir("C:\\Users\\LENOVO\\Downloads\\Gene_Fusion\\Gene_Fusion\\final_files")
#taking modified files to merge 
path = os.getcwd()+"/rna_final_files"
path
file = glob.glob(os.path.join(path,"*.csv"))
#print(file)
#print(path)
df_2=[]
for fp_1 in file:
  df_1 = pd.read_csv(fp_1, sep=',')
  df_2.append(df_1)
sample_all_1 = reduce(lambda left,right: pd.merge(left,right,on='#FusionGene',how='outer'), df_2)
#droping unwanted cloumns
sample_all_1 = sample_all_1[sample_all_1.columns.drop(list(sample_all_1.filter(regex='Unnamed')))]
#print(sample_all_1.columns)
sample_all_1.to_csv(path+r'/final_2/'+'gene_fusion_batch_final.csv', index = False)
###################################### End PART II ################################

###################################### Start PART III ######################################
path = os.getcwd()+"/rna_final_files/final_2"
os.chdir(path)

#df_fus = pd. read_csv("Fusion_list.csv",header='infer')
df_fus = pd. read_csv("Fusion_list.csv")
df_fus.columns = ['#FusionGene']
df_fus = df_fus[df_fus.columns.drop(list(df_fus.filter(regex='Unnamed')))]
df_fus = df_fus['#FusionGene'].astype(str)

df_scan = pd.read_csv('gene_fusion_batch_final.csv', sep=',')
filtered=pd.DataFrame(columns=df_scan.columns)
#taking gene list and scan it against the batch file created in PART II section 
for i in range(len(df_fus)):
    for j in df_fus[[i]]:
        gene_name_str = str(j)
    filtered = pd.concat([filtered, df_scan[df_scan['#FusionGene'].str.contains(gene_name_str) == True]])
filtered_1=filtered.reset_index(drop = True)

###################################### End PART III ######################################

###################################### Start PART IV ######################################
#gene_list = pd.read_csv("Genes_with_Test_name_Format_2.csv",header='infer')
gene_list = pd.read_csv("Genes_with_Test_name_Format_2.csv")
gene_list.columns = ["TarGT_CORE","TarGT_Abs","TarGT_INDIEGENE","TarGT_FIRST"]
columns_add = ['TC', 'TA', 'TI','TF']
for newcol in columns_add:
    filtered_1[newcol]= ''
for i in range(len(gene_list)):
    gene_name= gene_list.loc[i, "TarGT_CORE"]
    gene_name_str = str(gene_name)
    filtered_1['TC'][filtered_1['#FusionGene'].str.contains(gene_name_str)]=" TarGT_CORE " 
    gene_name= gene_list.loc[i, "TarGT_Abs"]
    gene_name_str = str(gene_name)
    filtered_1['TA'][filtered_1['#FusionGene'].str.contains(gene_name_str)]=" TarGT_Abs "
    gene_name= gene_list.loc[i, "TarGT_INDIEGENE"]
    gene_name_str = str(gene_name)
    filtered_1['TI'][filtered_1['#FusionGene'].str.contains(gene_name_str)]=" TarGT_INDIEGENE "
    gene_name= gene_list.loc[i, "TarGT_FIRST"]
    gene_name_str = str(gene_name)
    filtered_1['TF'][filtered_1['#FusionGene'].str.contains(gene_name_str)]=" TarGT_FIRST "
    
filtered_1['GENE_ASSOCIATED_TO_TEST']=filtered_1['TC'].map(str)+filtered_1['TA'].map(str)+filtered_1['TI'].map(str)+filtered_1['TF'].map(str)
filtered_1=filtered_1.drop(['TC', 'TA','TI','TF'], axis = 1)
filtered_1.to_csv('final_filtered_FUS.csv', index = False)
###################################### End PART IV ######################################
