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
import ipdb
import platform
opj = os.path.join
deb = ipdb.set_trace

# MRI imports
# -----------
import nibabel as nb
import cortex
from cortex.fmriprep import *
from nilearn import image

# Functions import
# ----------------
from utils import set_pycortex_config_file, convert_fit_results

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

# Define folder
# -------------
fmriprep_dir = "{base_dir}/deriv_data/fmriprep/".format(base_dir = base_dir)
fs_dir = "{base_dir}/deriv_data/fmriprep/freesurfer/".format(base_dir = base_dir)
xfm_name = 'identity.t1w'
cortex_dir = "{base_dir}/pp_data/cortex/db/{subject}".format(base_dir = base_dir, subject = subject)

# Set pycortex db and colormaps
# -----------------------------
set_pycortex_config_file(base_dir)

# Add participant to pycortex db
# ------------------------------
if os.path.isdir(cortex_dir) == False:
    print('import subject in pycortex')
    cortex.fmriprep.import_subj(subject = subject[-2:], source_dir = fmriprep_dir, sname = subject)
    cortex.freesurfer.import_flat(subject = subject, patch = 'full', freesurfer_subject_dir = fs_dir, sname = subject)

# Add transform to pycortex db
# ------------------------------
t1w = cortex.db.get_anat(subject)
transform = cortex.xfm.Transform(np.identity(4), t1w)
transform.save(subject, xfm_name, 'magnet')
    
# Resample to t1w
# ---------------
for mask_dir in ['all','pos','neg']:

    print('resample deriv mat: {mask_dir}'.format(mask_dir = mask_dir))
    
    # Load data
    deriv_mat_file = "{deriv_dir}/{mask_dir}/prf_deriv_{acq}_{mask_dir}.nii.gz".format(deriv_dir = deriv_dir, acq = acq, mask_dir = mask_dir)
    img_deriv_mat = nb.load(deriv_mat_file)

    # interpolate to t1w cortex
    deriv_rs = image.resample_to_img(source_img = img_deriv_mat, target_img = t1w, interpolation = 'nearest')
    
    # save resample file
    deriv_mat_rs_file = "{deriv_dir}/{mask_dir}/prf_deriv_{acq}_{mask_dir}_rs.nii.gz".format(deriv_dir = deriv_dir, acq = acq, mask_dir = mask_dir)
    new_img = nb.Nifti1Image(dataobj = deriv_rs.get_data(), affine = deriv_rs.affine, header = deriv_rs.header)
    new_img.to_filename(deriv_mat_rs_file)
    
# Resample tc file
# ----------------

print('resample tc file')

# Load data
tc_file = "{base_dir}/pp_data/{subject}/func/{subject}_task-AttendStim_{acq}_fmriprep_sg_psc_avg.nii.gz".format(base_dir = base_dir, subject = subject, acq = acq)
img_tc = nb.load(tc_file)

# interpolate to t1w cortex
tc_rs = image.resample_to_img(source_img = img_tc, target_img = t1w, interpolation = 'nearest')

# save resample file
tc_rs_file = "{base_dir}/pp_data/{subject}/func/{subject}_task-AttendStim_{acq}_fmriprep_sg_psc_avg_rs.nii.gz".format(base_dir = base_dir, subject = subject, acq = acq)
new_img = nb.Nifti1Image(dataobj = tc_rs.get_data(), affine = tc_rs.affine, header = tc_rs.header)
new_img.to_filename(tc_rs_file)