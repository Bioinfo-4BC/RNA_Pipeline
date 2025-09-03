import os
import pandas as pd
import glob
import sys
# load module
homefolderpath = "/home/ubuntu/Programs/NGS3Pipeline/RNA/"
sys.path.append(homefolderpath)
import codes.RNA_module as RNA_module

project_name = sys.argv[1]
sample_list = pd.read_csv("list.txt", sep="\t", header=None)
path = os.getcwd()
# Batch ID for QC
dirname = sys.argv[2]
bspath = sys.argv[3]



path_cov = os.path.join(path, "Target_coverage.txt")
pattern = "*stats.txt"
os.system("mkdir "+dirname)
# DRAGEN QC variables
tir = ('DRAGEN mapping_mqc-generalstats-dragen_mapping-Total_input_reads')
dups = ('DRAGEN mapping_mqc-generalstats-dragen_mapping-Number_of_duplicate_marked_reads_pct')
insertL = ('DRAGEN mapping_mqc-generalstats-dragen_mapping-Insert_length_median')

# FinalGUI2 script
for sample in sample_list[0]:
    RNA_module.copy_rna_files(bspath, project_name, sample)

os.system("mkdir FUS_P")
prelim_files = glob.glob(os.path.join(path,"*.preliminary"))
for fp in prelim_files:
    df = pd.read_table(fp, sep='\t')
    df["Score"] = df["Score"].astype(float)
    df["Score"] = df["Score"].round(5)
    df["Score"] = df["Score"].astype(str)
    f_name = os.path.basename(fp).split('.')[0]
    df.to_excel(path + r'/FUS_P/' + f_name + '_FUS_P' + '.xlsx', index=False)

gs_files = glob.glob(os.path.join(path,"*_general_stats.txt"))
for fg in gs_files:
    df = pd.read_csv(fg, sep="\t")
    df1 = df.transpose()
    f_name = os.path.basename(fg).split('.')[0]
    df1.to_excel(path + r'/' + f_name + '.xlsx', index=False)

# CRF script
crf_dir = os.path.join(path, "Common_RNA_Fusion")
os.system("mkdir Common_RNA_Fusion")
# rna files list
rna_files = glob.glob(os.path.join(path, "*.preliminary"))
os.system("cp " + " ".join(rna_files) + " Common_RNA_Fusion")
#os.system("find . -path './ignore_this_dir' -prune -o \\( -name '*preliminary' -o -name '*one' \\) -type f -exec cp {} Common_RNA_Fusion \\;")
source_dir= homefolderpath + "/Common_RNA_Fusion/input/"
os.system("find "+source_dir+" \( -mindepth 1 -maxdepth 1 -type d -o -name '*.py' \) -exec cp -r {} Common_RNA_Fusion \;")
os.chdir(crf_dir)
os.system("python3 Common_Fusion_Program.py")
os.chdir(path)

# coverconcat script
RNA_module.coverconcat()
# check if there is only one sample
if len(sample_list) == 1:
    onlyID = sample_list[0]
    # open Target_coverage.txt as a tsv file
    tcv = pd.read_csv("Target_coverage.txt", sep="\t", header = 0)
    # add the sample ID
    tcv['Sample'] = onlyID
    # remove the file and write the new one
    os.remove("Target_coverage.txt")
    tcv.to_csv("Target_coverage.txt", sep="\t", index=False)

RNA_module.fixcovfile("Target_coverage.txt")

# RNA QC script
pattern_files = glob.glob(os.path.join(path, pattern))

dfx = []

for fs in pattern_files:
    df = pd.read_csv(fs, sep="\t")
    dfx.append(df)
    qc_genstats = pd.concat(dfx)

genfile = dirname + "-genstats_QC.csv"
qc_genstats.to_csv(dirname+"/"+genfile, index=False)

coverage = pd.read_csv(path_cov, sep="\t", engine='python')
df = pd.merge(qc_genstats, coverage, on='Sample', how='inner')

# Scores

df.loc[df[tir] >= 20000000, 'TIP_SC'] = 1
df.loc[df[dups] <= 55, 'DUPS_SC'] = 1
df.loc[df[insertL] >= 75, 'INSERT_SC'] = 1
df.loc[df['TargetCoverage'] >= 30, 'COV_SC'] = 1

df.loc[df[tir] < 20000000, 'TIP_SC'] = 0
df.loc[df[dups] > 55, 'DUPS_SC'] = 0
df.loc[df[insertL] < 75, 'INSERT_SC'] = 0
df.loc[df['TargetCoverage'] < 30, 'COV_SC'] = 0

sc_cols = ['TIP_SC', 'DUPS_SC', 'INSERT_SC', 'COV_SC']
df['SCORE'] = df[sc_cols].sum(axis=1)

# Status
df.loc[df['SCORE'] >= 4, 'STATUS'] = 'PASS'
df.loc[df['SCORE'] == 3, 'STATUS'] = 'RECONSIDER'
df.loc[df['COV_SC'] == 0, 'STATUS'] = 'FAIL'
df.loc[df['INSERT_SC'] == 0, 'STATUS'] = 'FAIL'
df.loc[df['TIP_SC'] == 0, 'STATUS'] = 'FAIL'

# QC Comments
df.loc[(df['STATUS'] == 'RECONSIDER') & (df[tir] > 40000000) & (df[dups] < 60) & (df[dups] > 55) & (df['COV_SC'] == 1), 'REMARKS'] = 'ADEQUATE READS & COVERAGE'
df.loc[(df['STATUS'] == 'RECONSIDER') & (df[tir] > 40000000) & (df[dups] < 60) & (df[dups] > 55) & (df['COV_SC'] == 1), 'STATUS'] = 'PASSED AFTER RECONSIDERATION'
df.loc[(df['STATUS'] == 'RECONSIDER') & (df[tir] < 40000000) & (df[dups] < 60) & (df[dups] > 55) & (df['COV_SC'] == 1), 'REMARKS'] = 'LOW INPUT READS'
df.loc[(df['STATUS'] == 'RECONSIDER') & (df[tir] < 40000000) & (df[dups] < 60) & (df[dups] > 55) & (df['COV_SC'] == 1), 'STATUS'] = 'FAILED AFTER RECONSIDERATION'
df.loc[(df['STATUS'] == 'RECONSIDER') & (df[tir] > 40000000) & (df[dups] < 60) & (df[dups] > 55) & (df['COV_SC'] == 0), 'REMARKS'] = 'LOW COVERAGE'
df.loc[(df['STATUS'] == 'RECONSIDER') & (df[tir] > 40000000) & (df[dups] < 60) & (df[dups] > 55) & (df['COV_SC'] == 0), 'STATUS'] = 'FAILED AFTER RECONSIDERATION'

##LEVEL2, DUPS 60-65##----------------------------------------------------------
df.loc[(df['STATUS'] == 'RECONSIDER') & (df[tir] > 60000000) & (df[dups] < 65) & (df[dups] > 60) & (df['COV_SC'] == 1), 'REMARKS'] = 'ADEQUATE READS & COVERAGE'
df.loc[(df['STATUS'] == 'RECONSIDER') & (df[tir] > 60000000) & (df[dups] < 65) & (df[dups] > 60) & (df['COV_SC'] == 1), 'STATUS'] = 'PASSED AFTER RECONSIDERATION'
df.loc[(df['STATUS'] == 'RECONSIDER') & (df[tir] < 60000000) & (df[dups] < 65) & (df[dups] > 60) & (df['COV_SC'] == 1), 'REMARKS'] = 'LOW INPUT READS'
df.loc[(df['STATUS'] == 'RECONSIDER') & (df[tir] < 60000000) & (df[dups] < 65) & (df[dups] > 60) & (df['COV_SC'] == 1), 'STATUS'] = 'FAILED AFTER RECONSIDERATION'
df.loc[(df['STATUS'] == 'RECONSIDER') & (df[tir] > 60000000) & (df[dups] < 65) & (df[dups] > 60) & (df['COV_SC'] == 0), 'STATUS'] = 'FAILED AFTER RECONSIDERATION'
df.loc[(df['STATUS'] == 'RECONSIDER') & (df[tir] > 60000000) & (df[dups] < 65) & (df[dups] > 60) & (df['COV_SC'] == 0), 'REMARKS'] = 'LOW COVERAGE'

##LEVEL3, DUPS 65-70##----------------------------------------------------------
df.loc[(df['STATUS'] == 'RECONSIDER') & (df[tir] > 80000000) & (df[dups] < 70) & (df[dups] > 65) & (df['COV_SC'] == 1), 'REMARKS'] = 'ADEQUATE READS & COVERAGE'
df.loc[(df['STATUS'] == 'RECONSIDER') & (df[tir] > 80000000) & (df[dups] < 70) & (df[dups] > 65) & (df['COV_SC'] == 1), 'STATUS'] = 'PASSED AFTER RECONSIDERATION'
df.loc[(df['STATUS'] == 'RECONSIDER') & (df[tir] < 80000000) & (df[dups] < 70) & (df[dups] > 65) & (df['COV_SC'] == 1), 'REMARKS'] = 'LOW INPUT READS'
df.loc[(df['STATUS'] == 'RECONSIDER') & (df[tir] < 80000000) & (df[dups] < 70) & (df[dups] > 65) & (df['COV_SC'] == 1), 'STATUS'] = 'FAILED AFTER RECONSIDERATION'
df.loc[(df['STATUS'] == 'RECONSIDER') & (df[tir] > 80000000) & (df[dups] < 70) & (df[dups] > 65) & (df['COV_SC'] == 0), 'REMARKS'] = 'LOW COVERAGE'
df.loc[(df['STATUS'] == 'RECONSIDER') & (df[tir] > 80000000) & (df[dups] < 70) & (df[dups] > 65) & (df['COV_SC'] == 0), 'STATUS'] = 'FAILED AFTER RECONSIDERATION'

##FINAL REMARKS##==============================================================#
##FINAL ANNOTATION
df.loc[df['INSERT_SC'] == 0, 'REMARKS'] = 'INSERT LENGTH LESS THAN 75'
df.loc[df['COV_SC'] == 0, 'REMARKS'] = 'LOW COVERAGE'
df.loc[df['TIP_SC'] == 0, 'REMARKS'] = 'READS LESS THAN 20M'
df.loc[df['SCORE'] == 4, 'REMARKS'] = 'PASSED ALL QC CHECKPOINTS'
df.loc[(df['TIP_SC'] == 0) & (df['COV_SC'] == 0), 'REMARKS'] = 'LOW READS & COVERAGE'
df.loc[(df['TIP_SC'] == 0) & (df['INSERT_SC'] == 0), 'REMARKS'] = 'LOW READS & INSERT LENGTH LESS THAN 75'
df.loc[(df['COV_SC'] == 0) & (df['INSERT_SC'] == 0), 'REMARKS'] = 'LOW COVERAGE & INSERT LENGTH LESS THAN 75'
df.loc[df['SCORE'] == 0, 'REMARKS'] = 'FAILED ALL QC CHECKPOINTS'
df.loc[(df['SCORE'] == 3) & (df[dups] > 70), 'REMARKS'] = 'DUPS > 70'
df.loc[df['REMARKS'] == 'DUPS > 70', 'STATUS'] = 'FAILED AFTER RECONSIDERATION'

t_sams = len(df)
passed = (df['STATUS']=='PASS').sum()
failed = (df['STATUS']=='FAIL').sum()
recon_p = (df['STATUS']=='PASSED AFTER RECONSIDERATION').sum()
recon_f = (df['STATUS']=='FAILED AFTER RECONSIDERATION').sum()

data = [["PASSED", passed], ["FAILED", failed], ["FAILED AFTER RECON", recon_f], ["PASSED AFTER RECON", recon_p], ["TOTAL SAMPLES", t_sams]]

# Output
newname = dirname + "RAW-RNA-QC.csv"
df.to_csv(dirname+"/"+newname, index=False)

df_mod = df.drop(['TIP_SC', 'DUPS_SC', 'INSERT_SC', 'COV_SC', 'SCORE'], axis=1)

modname = dirname + "RNA-QC.xlsx"
df_mod = df_mod[['Sample', 'DRAGEN mapping_mqc-generalstats-dragen_mapping-Total_input_reads', 'DRAGEN mapping_mqc-generalstats-dragen_mapping-Number_of_duplicate_marked_reads_pct', 'DRAGEN mapping_mqc-generalstats-dragen_mapping-Insert_length_median', 'TargetCoverage', 'STATUS', 'REMARKS', 'DRAGEN mapping_mqc-generalstats-dragen_mapping-Average_sequenced_coverage_over_genome', 'DRAGEN mapping_mqc-generalstats-dragen_mapping-Reads_with_mate_sequenced_pct', 'DRAGEN mapping_mqc-generalstats-dragen_mapping-QC_failed_reads_pct', 'DRAGEN mapping_mqc-generalstats-dragen_mapping-Mapped_reads_pct', 'DRAGEN mapping_mqc-generalstats-dragen_mapping-Unmapped_reads_pct', 'DRAGEN mapping_mqc-generalstats-dragen_mapping-Number_of_unique_mapped_reads_excl_duplicate_marked_reads_pct', 'DRAGEN mapping_mqc-generalstats-dragen_mapping-Properly_paired_reads_pct', 'DRAGEN mapping_mqc-generalstats-dragen_mapping-Not_properly_paired_reads_discordant_pct', 'DRAGEN mapping_mqc-generalstats-dragen_mapping-Paired_reads_mapped_to_different_chromosomes_MAPQ_10_pct', 'DRAGEN mapping_mqc-generalstats-dragen_mapping-Q30_bases_pct', 'DRAGEN mapping_mqc-generalstats-dragen_mapping-Total_alignments', 'DRAGEN mapping_mqc-generalstats-dragen_mapping-Secondary_alignments_pct', 'DRAGEN coverage_mqc-generalstats-dragen_coverage-Aligned_reads', 'DRAGEN coverage_mqc-generalstats-dragen_coverage-Aligned_bases', 'DRAGEN coverage_mqc-generalstats-dragen_coverage-Average_alignment_coverage_over_genome', 'DRAGEN coverage_mqc-generalstats-dragen_coverage-Uniformity_of_coverage_PCT_0_2_mean_over_genome', 'DRAGEN coverage_mqc-generalstats-dragen_coverage-Mean_Median_autosomal_coverage_ratio_over_genome', 'DRAGEN coverage_mqc-generalstats-dragen_coverage-PCT_of_genome_with_coverage_1x_inf', 'DRAGEN coverage_mqc-generalstats-dragen_coverage-PCT_of_genome_with_coverage_20x_inf', 'DRAGEN coverage_mqc-generalstats-dragen_coverage-PCT_of_genome_with_coverage_50x_inf', 'DRAGEN coverage_mqc-generalstats-dragen_coverage-PCT_of_genome_with_coverage_100x_inf']]
df_mod.to_excel(dirname+"/"+modname, index=False)

# New RNA QC output

newqcname = dirname + "RNA-QC-upd.xlsx"
df_mod2 = pd.DataFrame()
df_mod2['Basespace_Sample_ID'] = df_mod['Sample']
# DNA parameters are all NA
df_mod2['Total_Sequences'] = 'NA'
df_mod2['PCT_DUP_AR'] = 'NA'
df_mod2['PCT_TC_50X'] = 'NA'
df_mod2['MTCD'] = 'NA'
df_mod2['UCT_02'] = 'NA'
df_mod2['UE_Base'] = 'NA'
df_mod2['UE_Read'] = 'NA'
df_mod2['FLM'] = 'NA'
df_mod2['PCT_Q30'] = 'NA'
df_mod2['Total_size_GB'] = 'NA'
df_mod2['DNA_QC_Status'] = 'NA'
df_mod2['DNA_QC_Comment'] = 'NA'
# RNA parameters from the data
df_mod2['RNA_Total_input_reads'] = df_mod['DRAGEN mapping_mqc-generalstats-dragen_mapping-Total_input_reads']
df_mod2['RNA_PCT_DUP_AR'] = df_mod['DRAGEN mapping_mqc-generalstats-dragen_mapping-Number_of_duplicate_marked_reads_pct']
df_mod2['RNA_Insert_length_median'] = df_mod['DRAGEN mapping_mqc-generalstats-dragen_mapping-Insert_length_median']
df_mod2['RNA_TargetCoverage'] = df_mod['TargetCoverage']
df_mod2['RNA_Total_size_GB'] = df_mod['DRAGEN mapping_mqc-generalstats-dragen_mapping-Total_input_reads'] * 150/1e9
df_mod2['QC_RNA_Status'] = df_mod['STATUS']
df_mod2['RNA_Comment'] = df_mod['REMARKS']

# Output to excel file
df_mod2.to_excel(dirname+"/"+newqcname, index=False)
