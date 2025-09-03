#!/bin/bash

# import home path
homefolderpath="/home/ubuntu/Programs/NGS3Pipeline/RNA/"

# Source directory containing files and folders
source_dir="$homefolderpath/Common_RNA_Fusion/input/"

# Destination directory
destination_dir="./Common_RNA_Fusion/"

# Create destination directory if it doesn't exist
mkdir -p "$destination_dir"

# Copy files ending with "preliminary" and "one" and one folder to the destination directory
find . \( -name "*preliminary" -o -name "*one" \) -type f -exec cp {} "$destination_dir" \;
find "$source_dir" \( -mindepth 1 -maxdepth 1 -type d -o -name "*.py" \) -exec cp -r {} "$destination_dir" \;

echo "Files and folder copied successfully to $destination_dir"

cd Common_RNA_Fusion
python3 Common_Fusion_Program.py
