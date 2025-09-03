import pandas as pd
import os
import sys
import subprocess
import datetime
import time
import glob

# load module
homefolderpath = "/home/ubuntu/Programs/NGS3Pipeline/RNA/"
sys.path.append(homefolderpath)
import codes.RNA_module as RNA_module

# running details
wdir = os.getcwd()
foldername = sys.argv[1]
batch_ID = sys.argv[2]

# open the first csv at the current location
csv = glob.glob(os.path.join(wdir,"*.csv"))
maindf = pd.read_csv(csv[0], header=0)

bsinstall = RNA_module.bsinstall
terminal_emulator = "gnome-terminal"
# fastq work
# check for fastq files in wdir
fastqfiles = glob.glob(os.path.join(wdir, "*.fastq.gz"))
if not fastqfiles:
    filestodl = maindf["Raw_Name_R1_FastQ"].tolist() + maindf["Raw_Name_R2_FastQ"].tolist()
    RNA_module.awsdownload(foldername , filestodl)


# rename all fastq files
fastqlist = glob.glob(os.path.join(wdir, "*.fastq.gz"))

#for file in fastqlist:
    #RNA_module.renaming(file)

projID = str(maindf["Project_ID"][0])

os.system(bsinstall + "bs upload dataset --project="+ projID +" *R1_001.fastq.gz *R2_001.fastq.gz")
time.sleep(30)
# copy all scripts to the current directory
os.system("cp " + homefolderpath + "/codes/RNA_QC_P.py " + wdir)
os.system("cp " + homefolderpath + "/codes/RNA_DS_P.py " + wdir)
os.system("cp " + homefolderpath + "/codes/RNAshellv2.sh " + wdir)


projectname = maindf["Project_name"][0]
capkit = maindf["Capturing_Kit"].dropna().unique().tolist()
bedid = str(maindf["bed_num_ID"][0])
testname = str(maindf["Test_Name"][0])
sampleList = maindf["Sample_ID"].tolist()

# separate samples by sample type
samplesCT = []
for sample in sampleList:
    if "-CT" in sample:
        samplesCT.append(sample)
samplesST8 = []
for sample in sampleList:
    if "-ST8" in sample:
        samplesST8.append(sample)

BsampleListCT = " ".join(samplesCT)
BsampleListST8 = " ".join(samplesST8)

bscmdCT = ""
bscmdST8 = ""

#construct bscmd
current_date = datetime.datetime.now().strftime("%d_%b_%y")
formatted_date = os.path.splitext(current_date)[0]
# analysis names
appsessionnameCT = formatted_date + "_" + batch_ID + "_" + testname + "_CT_RNA"
appsessionnameST8 = formatted_date + "_" + batch_ID + "_" + testname + "_ST8_RNA"

# check sample type
if "CT" in capkit:
    sampletype = "CT"
    bscmdCT = bsinstall + "bs launch application -n \"DRAGEN RNA Pipeline\" --app-version 3.6.3 -o project-id:" + projID + " -o app-session-name:"+ appsessionnameCT +" -l "+ appsessionnameCT +" output_format:BAM -o coverage_list.coverage_bed_id:" + bedid + " -o sample-id:$bsidsCT -o ht-ref:hg19-altaware-cnv-anchor.v8 -o gene_fusion:1 -o quantification_checkbox:1 -o commandline-disclaimer:true"
else:
    print("No CT samples found.")
    bscmdCT = ""
if "ST8" in capkit:
    sampletype = "ST8"
    bscmdST8 = bsinstall + "bs launch application -n \"DRAGEN RNA Pipeline\" --app-version 3.6.3 -o project-id:" + projID + " -o app-session-name:"+ appsessionnameST8 +" -l "+ appsessionnameST8 +" output_format:BAM -o coverage_list.coverage_bed_id:" + bedid + " -o sample-id:$bsidsST8 -o ht-ref:hg19-altaware-cnv-anchor.v8 -o gene_fusion:1 -o quantification_checkbox:1 -o commandline-disclaimer:true"
else:
    print("No ST8 samples found.")
    bscmdST8 = ""

# create a list of sample names separated by spaces
BsampleList = " ".join(sampleList)
typelist = maindf["Capturing_Kit"].unique().tolist()
# open the shell file
shellfile = wdir + "/RNAshellv2.sh"
with open(shellfile, 'r') as file :
    filedata = file.read()

# replace the placeholders in the shell file
#filedata = filedata.replace('{{samplenames}}' , sampleList.strip("[]").replace("'","").replace(",",""))
filedata = filedata.replace('{{samplenamesCT}}' , BsampleListCT.strip("[]").replace("'","").replace(",",""))
filedata = filedata.replace('{{samplenamesST8}}' , BsampleListST8.strip("[]").replace("'","").replace(",",""))
filedata = filedata.replace('{{location}}', wdir)
filedata = filedata.replace('{{bscmdct}}', bscmdCT)
filedata = filedata.replace('{{bscmdst8}}', bscmdST8)
with open(shellfile, 'w') as file:
    file.write(filedata)

os.system("chmod +x " + shellfile)
os.system("bash " + shellfile)

os.system("fastqc *.fastq.gz")

stat_input = input("Is the analysis complete on Basespace? (y/n): ").strip().lower()

if stat_input == 'y':
    # run the QC pipelines
    # basemount in the current directory
    bs_file = open(os.path.join(wdir, "basespace.sh"), "w")
    bs_file.write("basemount basespace")
    bs_file.close()
    os.system("chmod 777 " + os.path.join(wdir, "basespace.sh"))
    os.chdir(wdir)

    p = subprocess.Popen(['bash', os.path.join(wdir, 'basespace.sh')],
        stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
        universal_newlines=True)
    output, error = p.communicate(input='yes\n')
    bspath = os.getcwd() + "/basespace"

    # print list.txt file
    listfile = os.path.join(wdir, "list.txt")
    with open(listfile, "w") as f:
        f.write("Sample_ID\n")
        for sample in sampleList:
            f.write(sample + "\n")
    # check if the basespace directory is created
    if not os.path.exists(bspath + "/Projects"):
        print("Error: Basespace directory not created, downstream cannot be run.")
        sys.exit(1)

    # run RNA_QC.py

    os.system("python3 RNA_QC_P.py " + projectname + " " + batch_ID + " " + bspath)
    time.sleep(5)
    os.system("python3 RNA_DS_P.py " + projID + " " + projectname + " " + bspath)
    time.sleep(5)
    p = subprocess.Popen(['bash', os.path.join(wdir, 'basespace.sh')],
        stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
        universal_newlines=True)
    output, error = p.communicate(input='yes\n')
    RNA_module.movesjsf(os.getcwd(), "SJ_Files", "SF_Files")



else:
    print("Exiting the script. RNA downstream may be run separately later.")
    sys.exit(1)

exit()
