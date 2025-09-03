import pandas as pd
import os
import sys
import subprocess
import datetime
import time
import shutil
import logging

# Adjustments for each system
# bs install path and home folder location in this module
# home folder path in CRF_script.sh

bsinstall = ""
homefolderpath = "/home/ubuntu/Programs/NGS3Pipeline/RNA/"

# check the status of analysis
def stat_check(appsessionname, outfile):
    # check if the analysis is complete on basespace
    command_stat = bsinstall + 'bs appsession list | grep ' + appsessionname + ' | awk \'{print $6}\'' ' >' 'out.txt'
    output_stat = subprocess.check_output(command_stat, shell=True)
    print(output_stat.decode().strip())
    file1 = open(outfile + 'out.txt', 'r')
    stat = file1.read()
    file1.close()
    return stat

# copy the RNA files from basespace for QC
def copy_rna_files(bspath, project_name, sample):
    os.system("cp " + bspath + "/Projects/"+project_name+"/AppResults/"+sample+"_*/Files/"+sample+".fusion_candidates.preliminary .")
    os.system("cp " + bspath + "/Projects/"+project_name+"/AppResults/"+sample+"_*/Files/multiqc_data/multiqc_general_stats.txt .")
    os.system("mv multiqc_general_stats.txt "+sample+"_multiqc_general_stats.txt")
    os.system("cp " + bspath + "/Projects/"+project_name+"/AppResults/"+sample+"_*/Files/"+sample+".qc-coverage-region-1_overall_mean_cov.csv .")
    os.system("cp " + bspath + "/Projects/"+project_name+"/AppResults/"+sample+"_*/Files/"+sample+".quant.sf .")
    os.system("cp " + bspath + "/Projects/"+project_name+"/AppResults/"+sample+"_*/Files/"+sample+".quant.genes.sf .")
    os.system("cp " + bspath + "/Projects/"+project_name+"/AppResults/"+sample+"_*/Files/*.SJ.out.tab .")
    return

# coverconcat.sh
def coverconcat():
    os.system("tail -n +1 *csv | grep '==' > samples.txt")
    os.system("tail -n +1 *csv | grep 'Average alignment coverage over' > values.txt")
    os.system("paste samples.txt values.txt > coverage.txt")
    os.system("sed -i 's/.qc-coverage-region-1_overall_mean_cov.csv <==//g;s/,/\t/;s/==> //g' coverage.txt")
    os.system("cut -f 1,3 coverage.txt > Target_coverage.txt")
    os.system("sed -i '1i Sample\tTargetCoverage' Target_coverage.txt")
    return

# retrieve bsIDs
def get_bsids(sample_ID):
    # read the terminal output
    os.system(bsinstall + "bs get biosample -n "+ sample_ID +" –terse | grep \"Id\" | head -1 | grep -Eo '[0-9]{1,}'")


    os.system(bsinstall + "bs get biosample -n "+ sample_ID +" –terse") >> "bsids.txt"

    bsid = os.system(bsinstall + "bs get biosample -n "+ sample_ID +" –terse | grep \"Id\" | head -1 | grep -Eo '[0-9]{1,}'")
    bsid.strip()
    return bsid

def statcheckloop(appsessname, outvar):
    output_stat = stat_check(appsessname)
    while "Complete" not in output_stat:
        outvar = False
        # wait for 5 minutes
        time.sleep(300)
        output_stat = stat_check(appsessname)
        # define exit condition
        if "Complete" in output_stat:
            outvar = True
            break
    return outvar

def renaming(file):
    # only consider the filename after the last /
    file = os.path.basename(file)

    if '_R1' in file:
        new_name = file.replace('_R1', '-S1_R1')
        new_name = new_name.replace("_001", "")
        new_name = new_name.replace("_", "-")
        new_name = new_name.replace("-R1.", "_S1_L001_R1_001.")
        os.rename(file, new_name)
    if '_R2' in file:
        new_name = file.replace('_R2', '-S1_R2')
        new_name = new_name.replace("_001", "")
        new_name = new_name.replace("_", "-")
        new_name = new_name.replace("-R2.", "_S1_L001_R2_001.")
        os.rename(file, new_name)
    return

def awsdownload(foldername, filelist):
    for file in filelist:
        awscmd = f'aws s3 cp s3://4basecare-data/' + foldername + '/Reads/ ./ --recursive --exclude "*" --include "' + file + '" --profile=wiprodata'
        subprocess.run(awscmd, shell=True)
    return

def fixcovfile(file_path):
    with open(file_path, "r") as f:
        lines = [line.strip() for line in f if line.strip()]

    header = lines[0]
    data_lines = lines[1:]

    # First line (usually the .csv metadata) contains the first value
    try:
        first_coverage = data_lines[0].split()[-1]
        float(first_coverage)  # ensure it's numeric
    except (IndexError, ValueError):
        raise ValueError("Couldn't extract first coverage value from the .csv line.")

    # Extract sample names and coverage values from remaining lines
    sample_names = []
    coverage_values = [first_coverage]  # start with the orphaned first value

    for line in data_lines[1:]:
        parts = line.split()
        if len(parts) == 2:
            sample_names.append(parts[0])
            coverage_values.append(parts[1])
        elif len(parts) == 1:
            sample_names.append(parts[0])
        else:
            continue  # skip malformed lines

    # Check that number of samples matches number of coverage values
    if len(sample_names) != len(coverage_values):
        raise ValueError(f"Sample count ({len(sample_names)}) does not match coverage count ({len(coverage_values)}).")

    # Write back the corrected file
    with open(file_path, "w") as f:
        f.write(f"{header}\n")
        for sample, cov in zip(sample_names, coverage_values):
            f.write(f"{sample}\t{cov}\n")

def movesjsf(location, sj_folder, sf_folder):
    os.makedirs(sj_folder, exist_ok=True)
    os.makedirs(sf_folder, exist_ok=True)
    for filename in os.listdir(os.getcwd()):
        if filename.endswith(".SJ.out.tab"):
            shutil.move(filename, os.path.join(sj_folder, filename))
            logging.info(f"Moved {filename} to {sj_folder}")
        elif filename.endswith(".sf"):
            shutil.move(filename, os.path.join(sf_folder, filename))
            logging.info(f"Moved {filename} to {sf_folder}")