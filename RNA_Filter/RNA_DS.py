# -*- coding: utf-8 -*-
"""
Main script to copy RNA_Filter_V2.py and RNA_Vis_V3.py to the current directory
and execute them sequentially.
"""

import subprocess
import os
import shutil
import sys
import logging

# Setup logging configuration
logging.basicConfig(
    filename="rna_ds.log",  # Log file name
    level=logging.INFO,              # Log level
    format="%(asctime)s - %(levelname)s - %(message)s",  # Log message format
)

def copy_scripts():
    
    script_sources = {
        "RNA_Filter_V2.py": "/home/bioinfo4/Patient_Samples/RNA_QC_25/RNA_Filter/RNA_Filter_V2.py",  # Replace with the actual path
        "RNA_Vis_V3.py": "/home/bioinfo4/Patient_Samples/RNA_QC_25/RNA_Filter/RNA_Vis_V3.py"                # Replace with the actual path
    }

    for script_name, script_path in script_sources.items():
        if os.path.isfile(script_path):
            shutil.copy(script_path, os.getcwd())
            logging.info(f"Copied {script_name} to the current directory.")
        else:
            logging.error(f"Error: {script_name} not found at {script_path}.")
            sys.exit(1)

def create_folders():
    """
    Create the required folders if they do not exist.
    """
    sj_folder = os.path.join(os.getcwd(), "SJ_Files")
    sf_folder = os.path.join(os.getcwd(), "SF_Files")

    os.makedirs(sj_folder, exist_ok=True)
    os.makedirs(sf_folder, exist_ok=True)

    logging.info(f"Created folders: {sj_folder}, {sf_folder}")
    return sj_folder, sf_folder

def move_files(sj_folder, sf_folder):
    """
    Move .out.tab files to SJ_Files and .sf files to SF_Files.
    """
    for filename in os.listdir(os.getcwd()):
        if filename.endswith(".SJ.out.tab"):
            shutil.move(filename, os.path.join(sj_folder, filename))
            logging.info(f"Moved {filename} to {sj_folder}")
        elif filename.endswith(".sf"):
            shutil.move(filename, os.path.join(sf_folder, filename))
            logging.info(f"Moved {filename} to {sf_folder}")

def main():
    if len(sys.argv) != 2:
        logging.error("Usage: python3 main.py <bs_path>")
        print("Usage: python3 main.py <bs_path>")
        sys.exit(1)

    bs_path = sys.argv[1]

    # Step 1: Copy the required scripts to the current directory
    logging.info("Copying scripts...")
    copy_scripts()

    # Step 2: Run RNA_Filter_V2.py
    logging.info("Executing RNA_Filter_V2.py...")
    try:
        subprocess.run(["python3", "RNA_Filter_V2.py"], check=True)
        logging.info("RNA_Filter_V2.py executed successfully.")
    except subprocess.CalledProcessError as e:
        logging.error(f"Error executing RNA_Filter_V2.py: {e}")
        sys.exit(1)

    # Step 3: Run RNA_Vis_V3.py with bs_path as argument
    logging.info("Executing RNA_Vis_V3.py...")
    try:
        subprocess.run(["python3", "RNA_Vis_V3.py", bs_path], check=True)
        logging.info("RNA_Vis_V3.py executed successfully.")
    except subprocess.CalledProcessError as e:
        logging.error(f"Error executing RNA_Vis_V3.py: {e}")
        sys.exit(1)

    # Step 4: Create folders for SJ_Files and SF_Files
    sj_folder, sf_folder = create_folders()

    # Step 5: Move the relevant files
    move_files(sj_folder, sf_folder)

if __name__ == "__main__":
    main()

