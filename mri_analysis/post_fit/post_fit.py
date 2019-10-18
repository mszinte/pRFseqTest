"""
-----------------------------------------------------------------------------------------
pp_roi.py
-----------------------------------------------------------------------------------------
Goal of the script:
Combine fit files
Compute pRF derivatives
-----------------------------------------------------------------------------------------
Input(s):
sys.argv[1]: subject number
sys.argv[2]: acquisition ('acq-2p5mm','acq-2mm')
sys.argv[3]: fit model ('gauss','css')
-----------------------------------------------------------------------------------------
Output(s):
None
-----------------------------------------------------------------------------------------
To run:
python post_fit/post_fit.py sub-01 acq-2p5mm gauss
python post_fit/post_fit.py sub-01 acq-2mm gauss
-----------------------------------------------------------------------------------------
"""

# Stop warnings
# -------------
import warnings
warnings.filterwarnings("ignore")

# General imports
# ---------------
import os
import sys
import json
import glob
import numpy as np
import matplotlib.pyplot as pl
import ipdb
import platform
opj = os.path.join
deb = ipdb.set_trace

# MRI imports
# -----------
import nibabel as nb

# Functions import
# ----------------
from utils import set_pycortex_config_file, convert_fit_results, draw_cortex_vertex

# Get inputs
# ----------
subject = sys.argv[1]
acq = sys.argv[2]
fit_model = sys.argv[3]

if fit_model == 'gauss': fit_val = 6
elif fit_model == 'css': fit_val = 7

# Define analysis parameters
# --------------------------
with open('settings.json') as f:
    json_s = f.read()
    analysis_info = json.loads(json_s)

# Define folder
# -------------
base_dir = analysis_info['base_dir']
deriv_dir = opj(base_dir,'pp_data',subject,fit_model,'deriv')

# Copy data from scratchw to scratch
# ----------------------------------
os.system("rsync -az --no-g --no-p --progress {scratchw}/ {scratch}".format(
                scratch = analysis_info['base_dir'],
                scratchw  = analysis_info['base_dir_westmere']))

# Check if all slices are present
# -------------------------------

# Original data to analyse
data_file  =  "{base_dir}/pp_data/{sub}/func/{sub}_task-AttendStim_{acq}_fmriprep_sg_psc_avg.nii.gz".format(
                                base_dir = base_dir,
                                sub = subject,
                                acq = acq)
img_data = nb.load(data_file)
data = img_data.get_fdata()

mask_file  =  "{base_dir}/pp_data/{sub}/func/{sub}_task-AttendStim_{acq}_fmriprep_mask_avg.nii.gz".format(
                                base_dir = base_dir,
                                sub = subject,
                                acq = acq)

img_mask = nb.load(mask_file)
mask = img_mask.get_fdata()
slices = np.arange(mask.shape[2])[mask.mean(axis=(0,1))>0]

est_files = []
miss_files_nb = 0
for slice_nb in slices:
    
    est_file = "{base_dir}/pp_data/{subject}/{fit_model}/fit/{subject}_task-AttendStim_{acq}_est_z_{slice_nb}.nii.gz".format(
                                base_dir = base_dir,
                                subject = subject,
                                fit_model = fit_model,
                                acq = acq,
                                slice_nb = slice_nb
                                )
    if os.path.isfile(est_file):
        if os.path.getsize(est_file) == 0:
            num_miss_part += 1 
        else:
            est_files.append(est_file)
    else:
        miss_files_nb += 1

if miss_files_nb != 0:
    sys.exit('%i missing files, analysis stopped'%miss_files_nb)

# Combine and save estimates
# --------------------------
print('combining est files')
ests = np.zeros((data.shape[0],data.shape[1],data.shape[2],fit_val))
for est_file in est_files:
    img_est = nb.load(est_file)
    est = img_est.get_fdata()
    ests = ests + est

# Save estimates data
estfn = "{base_dir}/pp_data/{subject}/{fit_model}/fit/{subject}_task-AttendStim_{acq}_est.nii.gz".format(
                                base_dir = base_dir,
                                subject = subject,
                                fit_model = fit_model,
                                acq = acq)

new_img = nb.Nifti1Image(dataobj = ests, affine = img_data.affine, header = img_data.header)
new_img.to_filename(estfn)

# Compute derived measures from prfs
# ----------------------------------
print('extracting pRF derivatives')
convert_fit_results(est_fn = estfn,
                    output_dir = deriv_dir,
                    stim_width = analysis_info['stim_width'],
                    stim_height = analysis_info['stim_height'],
                    fit_model = fit_model,
                    acq = acq)