"""
-----------------------------------------------------------------------------------------
pre_fit.py
-----------------------------------------------------------------------------------------
Goal of the script:
SG filter, PSC, AVG runs and combine data of both hemisphere
-----------------------------------------------------------------------------------------
Input(s):
sys.argv[1]: subject name
sys.argv[2]: start voxel index
sys.argv[3]: end voxel index
sys.argv[4]: data file path
sys.argv[5]: main directory
-----------------------------------------------------------------------------------------
Output(s):
# Preprocessed timeseries files
-----------------------------------------------------------------------------------------
To run:
python pre_fit/pre_fit.py
-----------------------------------------------------------------------------------------
"""

# Stop warnings
# -------------
import warnings
warnings.filterwarnings("ignore")

# General imports
# ---------------
import json
import sys
import os
import glob
import ipdb
import platform
import numpy as np
opj = os.path.join
deb = ipdb.set_trace

# MRI analysis imports
# --------------------
import nibabel as nb
from scipy.signal import savgol_filter

with open('settings.json') as f:
    json_s = f.read()
    analysis_info = json.loads(json_s)
trans_cmd = 'rsync -avuz --progress'

# Define cluster/server specific parameters
# -----------------------------------------
base_dir = analysis_info['base_folder'] 

# Copy files in pp_data folder
# ----------------------------
for sub_name in analysis_info['subject_list'] :

    dest_folder = "{base_dir}/pp_data/{sub}".format(base_dir = base_dir, sub = sub_name)
    try: os.makedirs(dest_folder)
    except: pass

    orig_folder = "{base_dir}/deriv_data/fmriprep/fmriprep/{sub}/func".format(base_dir = base_dir, sub = sub_name)
            
    for acq in analysis_info['acq']:
        for run in analysis_info['runs']:
            orig_file = "{orig_fold}/{sub}_task-AttendStim_{acq}_{run}_space-T1w_desc-preproc_bold.nii.gz".format(
                                            orig_fold = orig_folder, sub = sub_name, acq = acq, run = run)
            dest_file = "{dest_fold}/{sub}_task-AttendStim_{acq}_{run}_fmriprep.nii.gz".format(
                                            dest_fold = dest_folder, sub = sub_name, acq = acq, run = run)

            os.system("{cmd} {orig} {dest}".format(cmd = trans_cmd, orig = orig_file, dest = dest_file))
            

# SG + PSC + AVG
# --------------

for sub_name in analysis_info['subject_list'] :
    
    for acq in analysis_info['acq']:
        # SG + PSC
        # --------
        print("{sub}-{acq}: sg + psc".format(sub = sub_name, acq = acq))
    
        file_list = sorted(glob.glob("{base_dir}/pp_data/{sub}/*{acq}*fmriprep.nii.gz".format(
                                            base_dir = base_dir, sub = sub_name, acq = acq)))
        
        for file in file_list:
        
            print(file)
            # load
            img = nb.load(file)
            pp_data = img.get_fdata()
        
            # sg filter
            pp_data_filt = savgol_filter( x = pp_data,
                                          window_length = analysis_info['sg_filt_window_length'],
                                          polyorder = analysis_info['sg_filt_polyorder'],
                                          deriv = analysis_info['sg_filt_deriv'],
                                          axis = 3, 
                                          mode = 'nearest')
            
            pp_data_filt_mean = pp_data_filt.mean(axis=3)
            pp_data_filt_mean = np.repeat(pp_data_filt_mean[:, :, :, np.newaxis], analysis_info['TR'], axis=3)
            pp_data_sg = pp_data - pp_data_filt + pp_data_filt_mean
            
            # percent signal change            
            pp_data_sg_median = np.median(pp_data_sg, axis=3)
            pp_data_sg_median = np.repeat(pp_data_sg_median[:, :, :, np.newaxis], analysis_info['TR'], axis=3)
            pp_data_sg_psc = 100.0 * (pp_data_sg - pp_data_sg_median)/pp_data_sg_median
            
            # save
            new_file = "{file}_sg_psc.nii.gz".format(file = file[:-7])
            new_img = nb.Nifti1Image(dataobj = pp_data_sg_psc, affine = img.affine, header = img.header)
            new_img.to_filename(new_file)
    
        # AVERAGE RUNS
        # ------------
        print("{sub}-{acq}: average runs".format(sub = sub_name, acq = acq))
        file_list = sorted(glob.glob("{base_dir}/pp_data/{sub}/*{acq}*_sg_psc.nii.gz".format(base_dir = base_dir, sub = sub_name, acq = acq)))
        img = nb.load(file_list[0])
        pp_data_sg_psc_avg = np.zeros(img.shape)
    
        for file in file_list:
            print(file)
            # load
            pp_data_sg_psc = []
            pp_data_sg_psc_img = nb.load(file)
            pp_data_sg_psc = pp_data_sg_psc_img.get_fdata()
            pp_data_sg_psc_avg += pp_data_sg_psc/len(file_list)
        
        # save
        new_file = "{file}_fmriprep_sg_psc_avg.nii.gz".format(file = file[:-30])
        new_img = nb.Nifti1Image(dataobj = pp_data_sg_psc_avg, affine = img.affine, header = img.header)
        new_img.to_filename(new_file)