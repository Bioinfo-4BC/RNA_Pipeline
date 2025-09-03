# -*- coding: utf-8 -*-
"""
Created on Mon Sep  9 14:49:05 2024

@author: Yukti, Naavarasi
@Modification: Naavarasi

Modified on Thur Jan 2 14:22:57 2025

Versions: Python 3.10, Pandas 1.5.3, Pysam 0.22.1 
"""
import pandas as pd
import os
import pysam  # for handling BAM files
import glob  # for expanding wildcards
import subprocess  # for executing shell commands
import shutil
import tempfile  # for handling file operations
import argparse  # for handling command-line arguments
import logging
import sys

# Setup logging configuration
logging.basicConfig(
    filename="rna_vis.log",  # Log file name
    level=logging.INFO,  # Log level
    format="%(asctime)s - %(levelname)s - %(message)s",  # Log message format
)

# Parse command-line arguments
parser = argparse.ArgumentParser(description="Process RNA fusion data.")
parser.add_argument("bs_path", help="Path of Basespace")
parser.add_argument("project_name", help="Project name")
args = parser.parse_args()

# Define paths
bs_path = args.bs_path
print(bs_path)
list_file = 'list.txt'
project_name = args.project_name

# Load sample IDs from list.txt
with open(list_file, 'r') as f:
    sample_ids = [line.strip() for line in f.readlines() if line.strip()]

# Create a main folder "RNA_html"
input_path = os.getcwd()
rna_html_dir = os.path.join(input_path, 'RNA_html')
os.makedirs(rna_html_dir, exist_ok=True)

logging.info(f"Created 'RNA_html' directory at {rna_html_dir}.")

# Create intermediate RNA_input.xlsx
data = []

for sample_id in sample_ids:
    # fusion_output_file = f"{sample_id}.fusions_output.xlsx"
    fusion_output_file = os.path.join(input_path, "Filtered_Fus_Updated", f"{sample_id}.fusions_output.xlsx")
    logging.info(f"Processing fusion output file: {fusion_output_file}")
    print(f"Processing fusion output file: {fusion_output_file}")

    if not os.path.exists(fusion_output_file):
        logging.warning(f"Fusion output file not found for {sample_id}. Skipping.")
        continue

    # Load the relevant subsheet from the Excel file
    try:
        xls = pd.ExcelFile(fusion_output_file)

        if "TI_gene_filter" in xls.sheet_names:
            df = pd.read_excel(xls, sheet_name="TI_gene_filter")
        elif "Fusion_Status" in xls.sheet_names:
            df = pd.read_excel(xls, sheet_name="Fusion_Status")
        else:
            logging.warning(
                f"Neither 'TI_gene_filter' nor 'Fusion_Status' subsheets found in {fusion_output_file}. Skipping.")
            continue

        # Filter rows with Score >= 0.85
        filtered_df = df[df['Score'] >= 0.85]

        # Check required columns exist
        required_columns = ['#FusionGene', 'Score', 'LeftBreakpoint', 'RightBreakpoint', 'ReadNames', 'ReadSupport']
        if not all(col in filtered_df.columns for col in required_columns):
            logging.warning(f"Missing required columns in {fusion_output_file}. Skipping.")
            continue

        # Append data to the intermediate list
        for _, row in filtered_df.iterrows():
            data.append({
                'Sample_ID': sample_id,
                '#FusionGene': row['#FusionGene'],
                'Score': row['Score'],
                'ReadSupport': row['ReadSupport']
            })

    except Exception as e:
        logging.error(f"Error processing {fusion_output_file}: {e}")
        print(f"Error processing {fusion_output_file}: {e}")

    # Create the RNA_input.xlsx file
    rna_input_path = os.path.join(rna_html_dir, 'RNA_input.xlsx')
    rna_input_df = pd.DataFrame(data, columns=['Sample_ID', '#FusionGene', 'Score', 'ReadSupport'])

    if not rna_input_df.empty:
        rna_input_df.to_excel(rna_input_path, index=False)
        print(f"RNA input file created: {rna_input_path}")
    else:
        logging.info("No data to write to RNA input file. Process complete.")

# Load the Excel file into a DataFrame
rna_ip = os.path.join(input_path, 'RNA_html', 'RNA_input.xlsx')
df = pd.read_excel(rna_ip)

# Display the first few rows of the DataFrame
print("Loaded RNA input data:")
print(f"\n{df.head()}")

# Iterate through each row of the DataFrame and process
for index, row in df.iterrows():
    sample_id = str(row[0])
    fusion = str(row[1])
    score = str(row[2])
    rs = str(row[3])

    logging.info(f"Processing Fusion {index + 1} for Sample ID: {sample_id}")
    logging.info(f"Fusion: {fusion}, Score: {score}, Read Support: {rs}")
    print(f"Processing Fusion {index + 1} for Sample ID: {sample_id} with {fusion}")

    # Create a directory for the current fusion
    fusion_folder = os.path.join(input_path, 'RNA_html', f'{sample_id}_{fusion}_{score}')
    os.makedirs(fusion_folder, exist_ok=True)
    logging.info(f"Created folder for fusion: {fusion_folder}")
    print(f"Created folder for fusion: {fusion_folder}")

    # Initialize matching_entries DataFrame for this fusion
    matching_entries = pd.DataFrame()

    # Use glob to find the BAM and prelim files
    bam_file_pattern = bs_path + f'/Projects/{project_name}/AppResults/{sample_id}*/Files/{sample_id}.bam'
    prelim_file_pattern = bs_path + f'/Projects/{project_name}/AppResults/{sample_id}*/Files/{sample_id}.fusion_candidates.preliminary'

    bam_files = glob.glob(bam_file_pattern)
    prelim_files = glob.glob(prelim_file_pattern)

    if not bam_files:
        logging.warning(f"BAM file not found for {sample_id}. Skipping.")
        print(f"BAM file not found for {sample_id}. Skipping.")
        continue

    if not prelim_files:
        logging.warning(f"Prelim file not found for {sample_id}. Skipping.")
        continue

    bam_file_path = bam_files[0]
    prelim_file_path = prelim_files[0]

    logging.info(f"Path of BAM: {bam_file_path}")
    logging.info(f"Path of prelim: {prelim_file_path}")

    # Accessing prelim file
    try:
        prelim_df = pd.read_csv(prelim_file_path, sep="\t")
        logging.info(f"Prelim file loaded into DataFrame")
    except Exception as e:
        logging.error(f"Error loading prelim file: {e}")
        continue

    # Search for the matching Fusion and Score in prelim DataFrame
    try:
        match = prelim_df[(prelim_df['#FusionGene'] == fusion) & (prelim_df['Score'] == float(score))]

        if not match.empty:
            for _, row in match.iterrows():
                left_bp = row['LeftBreakpoint']
                right_bp = row['RightBreakpoint']
                read_names = row['ReadNames']

                try:
                    chr1, pos1, dir1 = left_bp.split(':')
                    chr2, pos2, dir2 = right_bp.split(':')

                    pos1 = int(pos1)
                    pos2 = int(pos2)

                    start1 = pos1 - 250
                    end1 = pos1 + 250
                    start2 = pos2 - 250
                    end2 = pos2 + 250

                    new_entry_df = pd.DataFrame([{
                        'chr1': str(chr1),
                        'start1': start1,
                        'end1': end1,
                        'chr2': str(chr2),
                        'start2': start2,
                        'end2': end2,
                        'Fusion': f"{fusion} | {score} | {rs}"
                    }])

                    matching_entries = pd.concat([matching_entries, new_entry_df], ignore_index=True)

                    read_names_list = [name for name in read_names.split(';') if name]
                    if read_names_list:
                        with tempfile.NamedTemporaryFile(delete=False, mode='w') as temp_readnames_file:
                            for read_name in read_names_list:
                                temp_readnames_file.write(f"{read_name}\n")
                            temp_file_path = temp_readnames_file.name

                        intermediate_bam_path = os.path.join(fusion_folder,
                                                             f'{sample_id}_{fusion}_{score}_filtered.bam')
                        samtools_cmd = f"samtools view -N {temp_file_path} -b {bam_file_path} > {intermediate_bam_path}"

                        try:
                            subprocess.run(samtools_cmd, shell=True, check=True)
                            logging.info(f"Intermediate BAM file created for Fusion: {fusion}")

                            samtools_index_cmd = f"samtools index {intermediate_bam_path}"
                            subprocess.run(samtools_index_cmd, shell=True, check=True)
                            logging.info(f"Index file created for {intermediate_bam_path}")
                        except subprocess.CalledProcessError as e:
                            logging.error(f"Error running samtools command: {e}")
                    else:
                        logging.warning(f"No valid read names found for Fusion: {fusion}. Skipping BAM processing.")

                except ValueError as e:
                    logging.error(f"ValueError: {e}. Check the format of LeftBreakpoint or RightBreakpoint.")
                    continue

        else:
            logging.warning(f"No matching entry found for Fusion: {fusion} with Score: {score}")
    except KeyError as e:
        logging.error(f"KeyError: {e}. Check if column names are correct.")
    except ValueError as e:
        logging.error(f"ValueError: {e}. Check if the data types are correct.")

    bedpe_path = os.path.join(fusion_folder, f'{sample_id}_{fusion}_{score}_fusions.bedpe')
    matching_entries.to_csv(bedpe_path, sep='\t', index=False, header=False)
    logging.info(f"BEDPE file created at: {bedpe_path}")

    try:
        fasta_file = '/home/ubuntu/Freebayes_Test_Runs/Common_files/hg19.fasta'
        ideogram_file = '/home/ubuntu/Programs/IGV_Linux_2.9.2/igv/genomes/hg19/cytoBand.txt'
        gene_track = '/home/ubuntu/Programs/IGV_Linux_2.9.2/igv/genomes/hg19/ncbiRefGene.txt'
        report_file = os.path.join(fusion_folder, f'{sample_id}_{fusion}_{score}_report.html')

        create_report_cmd = (
            f" create_report"
            f" {bedpe_path} "
            f"--fasta {fasta_file} "
            f"--exclude-flag 512 "
            f"--ideogram {ideogram_file} "
            f"--flanking 300 "
            f"--tracks {intermediate_bam_path} {gene_track} "
            f"--output {report_file}"
        )

        subprocess.run(create_report_cmd, shell=True, check=True)
        logging.info(f"IGV report created at: {report_file}")
    except subprocess.CalledProcessError as e:
        logging.error(f"Error running IGV report command: {e}")
        print(f"Error running IGV report: {e}")

    # Move HTML report to a sample-wise folder and delete the fusion folder
    sample_folder = os.path.join(rna_html_dir, sample_id)
    os.makedirs(sample_folder, exist_ok=True)

    for file in os.listdir(fusion_folder):
        if file.endswith('.html'):
            shutil.move(os.path.join(fusion_folder, file), sample_folder)
            logging.info(f"Moved HTML report to: {sample_folder}")

    shutil.rmtree(fusion_folder)
    logging.info(f"Deleted fusion folder: {fusion_folder}")
    print(f"Process completed for {sample_id} with {fusion} \n")
logging.info("Processing complete.")
print("Processing complete.")

