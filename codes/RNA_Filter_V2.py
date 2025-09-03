# -*- coding: utf-8 -*-
"""
Created on Wed Jun 12 15:13:42 2024

@author: Naavarasi
"""

import pandas as pd
import os

# Define the input path
#input_path = r"C:\4BC\RNA_Filter\Example_Fus\CT"  # Replace with your actual directory path
input_path = os.getcwd()
output_folder = os.path.join(input_path, 'Filtered_Fus_Updated')
os.makedirs(output_folder, exist_ok=True)
DB_path = "/home/ubuntu/Programs/NGS3Pipeline/RNA/RNA_Filter/"

# Read list.txt to get the filenames
list_path = os.path.join(input_path, 'list.txt')
with open(list_path, 'r') as list_file:
    filenames = [line.strip() for line in list_file]

# Read known fusions and TI gene list
"""
TI_db_list_path = os.path.join(DB_path, 'TI_DB_list.txt')
with open(TI_db_list_path, 'r', encoding='utf-8') as file:
    known_fusions_CT = set(file.read().splitlines())

TA_db_list_path = os.path.join(DB_path, 'DB_list.txt')
with open(TA_db_list_path, 'r', encoding='utf-8') as file:
    known_fusions_ST8 = set(file.read().splitlines())
"""
ti_gene_list_path = os.path.join(DB_path, 'TI_gene_list.txt')
with open(ti_gene_list_path, 'r', encoding='utf-8') as ti_file:
    ti_gene_list = set(ti_file.read().splitlines())

# Load Fusion_DB.xlsx
fusion_db_path = os.path.join(DB_path, 'Fusion_DB.xlsx')
df_fusion_db = pd.read_excel(fusion_db_path, sheet_name=0)

# Split modified_fusion_pair into Gene1 and Gene2
df_fusion_db[['Gene1', 'Gene2']] = df_fusion_db['modified_fusion_pair'].str.split('--', expand=True)

# Read TI_gene_list.txt into a set
with open(ti_gene_list_path, 'r', encoding='utf-8') as ti_file:
    ti_gene_list = set(ti_file.read().splitlines())

# Filter modified fusion pairs where Gene1 or Gene2 is in TI gene list
filtered_fusion_pairs = df_fusion_db[
    (df_fusion_db['Gene1'].isin(ti_gene_list)) | (df_fusion_db['Gene2'].isin(ti_gene_list))
]['modified_fusion_pair'].unique()

# Set known_fusions_CT (filtered fusions)
known_fusions_CT = set(filtered_fusion_pairs)

# Set known_fusions_ST8 (all fusions from modified_fusion_pair)
known_fusions_ST8 = set(df_fusion_db['modified_fusion_pair'].unique())

# Function to check if a fusion is known
def check_fusion_status_TI(fusion):
    genes = fusion.split('--')
    reversed_fusion = f"{genes[1]}--{genes[0]}"
    if fusion in known_fusions_CT or reversed_fusion in known_fusions_CT:
        return 'KnownFusions'
    else:
        return 'Unknown'

def check_fusion_status_TA(fusion):
    genes = fusion.split('--')
    reversed_fusion = f"{genes[1]}--{genes[0]}"
    if fusion in known_fusions_ST8 or reversed_fusion in known_fusions_ST8:
        return 'KnownFusions'
    else:
        return 'Unknown'
    
# Function to check if fusion genes are in TI_gene_list
def check_ti_gene(fusion):
    genes = fusion.split('--')
    if genes[0] in ti_gene_list or genes[1] in ti_gene_list:
        return True
    else:
        return False

# Function to count ReadSupport
def count_read_support(read_names):
    if pd.isna(read_names) or read_names.strip() == '':
        return 0
    return len(read_names.split(';')) - 1

# Function to aggregate ReadSupport
def aggregate_read_support(df):
    df['ReadSupport'] = df['ReadNames'].apply(count_read_support)
    agg_df = df.groupby(['#FusionGene', 'Score', 'LeftBreakpoint', 'RightBreakpoint', 'FusionStatus'], as_index=False).agg({
        'ReadNames': ';'.join,
        'ReadSupport': 'sum'
    })
    return agg_df
"""      
def fetch_cancer_type_and_pmids(fusion):
    direct_match = df_fusion_db[df_fusion_db['modified_fusion_pair'] == fusion]
    reverse_match = df_fusion_db[df_fusion_db['modified_fusion_pair'] == '--'.join(reversed(fusion.split('--')))]
    match = pd.concat([direct_match, reverse_match], ignore_index=True)
    
    if not match.empty:
        cancer_type = match.iloc[0]['cancer_type']
        pmids = match.iloc[0]['pmids']
    else:
        cancer_type = None
        pmids = None
    
    return cancer_type, pmids
"""
def fetch_cancer_type_and_pmids(fusion):
    try:
        # Look for direct and reversed matches
        direct_match = df_fusion_db[df_fusion_db['modified_fusion_pair'] == fusion]
        reverse_match = df_fusion_db[df_fusion_db['modified_fusion_pair'] == '--'.join(reversed(fusion.split('--')))]
        match = pd.concat([direct_match, reverse_match], ignore_index=True)

        # Extract cancer type and PubMed IDs
        if not match.empty:
            cancer_type = match.iloc[0].get('cancer_type', pd.NA)
            #pmids = match.iloc[0].get('pmids', pd.NA)
        else:
            cancer_type = pd.NA
    except Exception as e:
        print(f"Error processing fusion: {fusion}. Error: {e}")
        cancer_type = pd.NA

    # Return as a tuple
    return (cancer_type)


for filename in filenames:
    base_filename = filename.split('.')[0]
    fusion_file = os.path.join(input_path, f"{base_filename}.fusion_candidates.preliminary")

    if not os.path.isfile(fusion_file):
        print(f"File not found: {fusion_file}")
        continue
    
    if '-CT' in filename:
        df_fusions_CT = pd.read_csv(fusion_file, sep='\t')
        df_all_fusions_CT = df_fusions_CT.copy()
        df_all_fusions_CT['FusionStatus'] = df_all_fusions_CT['#FusionGene'].apply(check_fusion_status_TI)

        # Step 4: Create df_sheet1
        df_sheet1a = df_fusions_CT.copy()

        # Step 5: Create df_sheet2 by checking against TI_gene_list.txt
        df_sheet2a = df_all_fusions_CT[df_all_fusions_CT['#FusionGene'].apply(check_ti_gene)].copy()

        # Add ReadSupport to df_sheet2
        df_sheet2a = aggregate_read_support(df_sheet2a)

        # Sort df_sheet2 by 'Score' from highest to lowest
        df_sheet2a = df_sheet2a.sort_values(by='Score', ascending=False)

        # Step 6: Create df_sheet3 with KnownFusions and filtering out the column 'Score' which is greater than and equal to 0.7
        df_sheet3a = df_sheet2a[(df_sheet2a['FusionStatus'] == 'KnownFusions') & (df_sheet2a['Score'] >= 0.7)].copy()

        # Add entries with the same fusion but a score less than 0.7 if they are repeated
        low_score_fusions_kf = df_sheet2a[(df_sheet2a['FusionStatus'] == 'KnownFusions') & (df_sheet2a['Score'] < 0.7)]
        for fusion in df_sheet3a['#FusionGene'].unique():
            repeated_entries = low_score_fusions_kf[low_score_fusions_kf['#FusionGene'] == fusion]
            df_sheet3a = pd.concat([df_sheet3a, repeated_entries], ignore_index=True)

        # Sort df_sheet3 by 'Score' from highest to lowest
        df_sheet3a = df_sheet3a.sort_values(by='Score', ascending=False)
        
        # Step 7: Create df_sheet3 with UnknownFusions and filtering out the column 'Score' which is greater than and equal to 0.7
        df_sheet4a = df_sheet2a[(df_sheet2a['FusionStatus'] == 'Unknown') & (df_sheet2a['Score'] >= 0.7)].copy()

        # Add entries with the same fusion but a score less than 0.7 if they are repeated
        low_score_fusions_uf = df_sheet2a[(df_sheet2a['FusionStatus'] == 'Unknown') & (df_sheet2a['Score'] < 0.7)]
        for fusion in df_sheet4a['#FusionGene'].unique():
            repeated_entries = low_score_fusions_uf[low_score_fusions_uf['#FusionGene'] == fusion]
            df_sheet4a = pd.concat([df_sheet4a, repeated_entries], ignore_index=True)

        # Sort df_sheet3 by 'Score' from highest to lowest
        df_sheet4a = df_sheet4a.sort_values(by='Score', ascending=False)
        
        # Add Cancer Type and PubMed ID columns to df_sheet3b
        """
        df_sheet3a[['Cancer Type', 'PubMed Id']] = df_sheet3a['#FusionGene'].apply(
            lambda fusion: pd.Series(fetch_cancer_type_and_pmids(fusion))
        )
        """
        df_sheet3a['Cancer Type'] = (df_sheet3a['#FusionGene'].apply(fetch_cancer_type_and_pmids).apply(pd.Series))

        # Save to an Excel file with multiple sheets
        output_file = os.path.join(output_folder, f"{base_filename}.fusions_output.xlsx")
        with pd.ExcelWriter(output_file, engine='openpyxl') as writer:
            df_sheet1a.to_excel(writer, sheet_name='All_fusions', index=False)
            df_sheet2a.to_excel(writer, sheet_name='TI_gene_filter', index=False)
            df_sheet3a.to_excel(writer, sheet_name='Known_Fusions', index=False)
            df_sheet4a.to_excel(writer, sheet_name='UnKnown_Fusions', index=False)
            
        print(f"Excel file '{output_file}' has been created successfully.")
    
    if '-ST8' in filename:
        df_fusions_ST8 = pd.read_csv(fusion_file, sep='\t')

        # Step 2: Access the '#FusionGenes' column and search against DB_list.txt
        df_all_fusions_ST8 = df_fusions_ST8.copy()
        df_all_fusions_ST8['FusionStatus'] = df_all_fusions_ST8['#FusionGene'].apply(check_fusion_status_TA)

        # Step 4: Create df_sheet1
        df_sheet1b = df_fusions_ST8.copy()

        # Add ReadSupport to df_sheet2
        df_sheet2b = df_all_fusions_ST8.copy()
        df_sheet2b = aggregate_read_support(df_sheet2b)

        # Sort df_sheet2 by 'Score' from highest to lowest
        df_sheet2b = df_sheet2b.sort_values(by='Score', ascending=False)

        # Step 6: Create df_sheet3 with KnownFusions and filtering out the column 'Score' which is greater than and equal to 0.7
        df_sheet3b = df_sheet2b[(df_sheet2b['FusionStatus'] == 'KnownFusions') & (df_sheet2b['Score'] >= 0.7)].copy()

        # Add entries with the same fusion but a score less than 0.7 if they are repeated
        low_score_fusions_k = df_sheet2b[(df_sheet2b['FusionStatus'] == 'KnownFusions') & (df_sheet2b['Score'] < 0.7)]
        for fusion in df_sheet3b['#FusionGene'].unique():
            repeated_entries = low_score_fusions_k[low_score_fusions_k['#FusionGene'] == fusion]
            df_sheet3b = pd.concat([df_sheet3b, repeated_entries], ignore_index=True)

        # Sort df_sheet3 by 'Score' from highest to lowest
        df_sheet3b = df_sheet3b.sort_values(by='Score', ascending=False)
        
        # Step 7: Create df_sheet3 with UnknownFusions and filtering out the column 'Score' which is greater than and equal to 0.7
        df_sheet4b = df_sheet2b[(df_sheet2b['FusionStatus'] == 'Unknown') & (df_sheet2b['Score'] >= 0.7)].copy()

        # Add entries with the same fusion but a score less than 0.7 if they are repeated
        low_score_fusions_u = df_sheet2b[(df_sheet2b['FusionStatus'] == 'Unknown') & (df_sheet2b['Score'] < 0.7)]
        for fusion in df_sheet4b['#FusionGene'].unique():
            repeated_entries = low_score_fusions_u[low_score_fusions_u['#FusionGene'] == fusion]
            df_sheet4b = pd.concat([df_sheet4b, repeated_entries], ignore_index=True)

        # Sort df_sheet3 by 'Score' from highest to lowest
        df_sheet4b = df_sheet4b.sort_values(by='Score', ascending=False)
        
        # Add Cancer Type and PubMed ID columns to df_sheet3b
        """
        df_sheet3b[['Cancer Type', 'PubMed Id']] = df_sheet3b['#FusionGene'].apply(
            lambda fusion: pd.Series(fetch_cancer_type_and_pmids(fusion))
        )
        """
        df_sheet3b['Cancer Type'] = (df_sheet3b['#FusionGene'].apply(fetch_cancer_type_and_pmids).apply(pd.Series))
        # Save to an Excel file with multiple sheets
        output_file = os.path.join(output_folder, f"{base_filename}.fusions_output.xlsx")
        with pd.ExcelWriter(output_file, engine='openpyxl') as writer:
            df_sheet1b.to_excel(writer, sheet_name='All_fusions', index=False)
            df_sheet2b.to_excel(writer, sheet_name='Fusion_Status', index=False)
            df_sheet3b.to_excel(writer, sheet_name='Known_Fusions', index=False)
            df_sheet4b.to_excel(writer, sheet_name='UnKnown_Fusions', index=False)
            
        print(f"Excel file '{output_file}' has been created successfully.")
        
    else:
           None

