# RNA_Pipeline
RNA_Pipeline_03/09/25

Documentation

Core files:
codes/RNA_module.py
ushellpipe.py
codes/RNA_QC_P.py
codes/RNA_DS_P.py
codes/RNA_Vis_V3_P.py

Installation:
Edit RNA_Module.py to set the path to the RNA pipeline directory.

Path edits in the following scripts:

codes/RNA_module.py

ushellpipe.py

codes/RNA_QC_P.py

Common_RNA_Fusion/CRF_script.sh

Ideogram and gene track file paths need to be changed in RNA_Vis_V3, as well as the hg19.fasta path
eg.

fasta_file = '/home/ubuntu/Freebayes_Test_Runs/Common_files/hg19.fasta'

ideogram_file = '/home/ubuntu/Programs/IGV_Linux_2.9.2/igv/genomes/hg19/cytoBand.txt'

gene_track = '/home/ubuntu/Programs/IGV_Linux_2.9.2/igv/genomes/hg19/ncbiRefGene.txt'

Execution:
python3 ushellpipe.py [arg1] [arg2]
arg1: AWS bucket ID
arg2: Batch ID (for analysis name and QC report folder)

Roadmap:
1. Open the csv file for project ID, capturing kit, sample list
2. Upload all samples to basespace in the specified project
3. Create the shell file
    3a. Check biosample IDs per capturing kit
    3b. initiate analysis per capturing kit
4. Check analysis status on bs, enter status as user input
5. If yes, mount basespace
6. Copy files from basespace to local directory, run QC
7. Create fusion files and HTML reports

Note for importing module:

Some systems don't accept import code.RNA_Module.
In order to get around this, add codes/ to the home directory path, remove codes. from the import statement across the scripts.
This may require editing the Filter database location.
