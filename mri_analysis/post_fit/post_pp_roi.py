"""
-----------------------------------------------------------------------------------------
post_pp_roi.py
-----------------------------------------------------------------------------------------
Goal of the script:
Create roi-masks and save all data in hdf5 format
-----------------------------------------------------------------------------------------
Input(s):
sys.argv[1]: subject name (e.g. 'sub-001')
-----------------------------------------------------------------------------------------
Output(s):
hdf5 files
-----------------------------------------------------------------------------------------
To run:
cd /home/szinte/projects/pRF_analysis/
python post_pp_roi.py sub-002
-----------------------------------------------------------------------------------------
"""

# General imports
# ---------------
import os
import shutil
import sys
import json
import numpy as np
import ipdb
import matplotlib.pyplot as plt
opj = os.path.join
deb = ipdb.set_trace

# MRI imports
# -----------
import nibabel as nb
import cortex

# Function import
# ---------------
from utils import set_pycortex_config_file

# Get inputs
# ----------
subject = sys.argv[1]
acq = sys.argv[2]

# Define analysis parameters
# --------------------------
with open('settings.json') as f:
    json_s = f.read()
    analysis_info = json.loads(json_s)

# Define folders and settings
# ---------------------------
base_dir = analysis_info['base_dir_local']
xfm_name = "identity.{acq}".format(acq = acq)
rois = analysis_info['rois']
gm_sampler = analysis_info['gm_sampler']

# Set pycortex db and colormaps
# -----------------------------
set_pycortex_config_file(base_dir)

# Create ROI masks
# ----------------
# deb()
# mask = cortex.get_cortical_mask(subject = subject, xfmname = xfm_name, type='cortical')
deb()
roi_masks = cortex.utils.get_roi_masks(subject = subject, xfmname = xfm_name, gm_sampler = gm_sampler, roi_list = 'V1',return_dict=True)
# index_volume, index_labels = cortex.utils.get_roi_masks(subject = subject, xfmname = xfm_name, gm_sampler = gm_sampler, roi_list = 'V1',return_dict=True)
deb()
#

# # Create HDF5 files
# # -----------------

# # define file name of data to include in hdf5
# loo_avg_files                   =   sorted(glob.glob(opj(analysis_info['aeneas_project_directory'],\
#                                                         'pp', sub_id, 'av_gaze*','loo','*.nii.gz')))                    # get averaged loo files

# loo_prf_files                   =   sorted(glob.glob(opj(analysis_info['aeneas_project_directory'],\
#                                                         'pp', sub_id, 'deriv_gazeAll','all','loo_*.nii.gz')))           # get loo pRF analysis including pos/neg pRF  

# loo_prf_files_left              =   sorted(glob.glob(opj(analysis_info['aeneas_project_directory'],\
#                                                         'pp', sub_id, 'deriv_gazeLeft','all','loo_*.nii.gz')))          # get loo pRF analysis including pos/neg pRF  

# loo_prf_files_right             =   sorted(glob.glob(opj(analysis_info['aeneas_project_directory'],\
#                                                         'pp', sub_id, 'deriv_gazeRight','all','loo_*.nii.gz')))         # get loo pRF analysis including pos/neg pRF  

# files                           =   [ (loo_avg_files,       'loo_avg'),      (loo_prf_files,        'loo_prf'), 
#                                       (loo_prf_files_left,  'loo_prf_left'), (loo_prf_files_right,  'loo_prf_right')]

# for roi_mask in analysis_info['pycortex_roi_mask_types']:
#     for roi in analysis_info['rois']:
#         mask_files                      =   sorted(glob.glob(opj(masks_folder,roi_mask,\
#                                                         '{roi}*.nii.gz'.format(roi = roi))))                            # get roi mask file name
#         h5_file                         =   opj(h5_folder,'{roi}.h5'.format(\
#                                                         roi = roi))                                                     # define h5 file

#         try: os.system('rm '+ h5_file)                                                                                  # delete old file
#         except: pass
        
#         for file, fName in files:
#             _                           =   mask_nii_2_hdf5(                                                            # add masked loo data
#                                             in_files        =   file,                                                   # absolute path to functional nifti-files.
#                                             mask_files      =   mask_files,                                             # absolute path to mask nifti-files
#                                             hdf5_file       =   h5_file,                                                # absolute path to hdf5 file
#                                             folder_alias    =   fName)                                                  # name of the to-be-created folder


