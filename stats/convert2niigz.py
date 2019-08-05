"""
-----------------------------------------------------------------------------------------
convert2niigz.py
-----------------------------------------------------------------------------------------
Goal of the script:
Convert PAR/REC to nifti (nii.gz) format
-----------------------------------------------------------------------------------------
Input(s):
sys.argv[1]: raw data directory
-----------------------------------------------------------------------------------------
Output(s):
nifti files
-----------------------------------------------------------------------------------------
To run:
ssh -Y compute-01
module load collections/default
cd /data1/projects/fMRI-course/Spinoza_Course/
python convert2niigz.py [your raw data path]
-----------------------------------------------------------------------------------------
Written by Martin Szinte (martin.szinte@gmail.com)
-----------------------------------------------------------------------------------------
"""

# imports modules
import sys
import os
import glob
import nibabel as nb
opj = os.path.join

# define subject folder
input_folder = sys.argv[1]

# create output folder
output_folder = opj(input_folder,'nifti')
try: os.makedirs(output_folder)
except: pass

# get PAR REC file list
list_par_files = glob.glob(opj(input_folder,'*.PAR'))

# convert files
print('convert files to nifti')

# using dcm2niix
cmd_txt = "dcm2niix -b y -o {out} -z y {in_folder} ".format(out = output_folder, in_folder = input_folder)
os.system(cmd_txt)
