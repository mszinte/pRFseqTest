"""
-----------------------------------------------------------------------------------------
roi_to_hdf5.py
-----------------------------------------------------------------------------------------
Goal of the script:
Create roi-masks and save all data in hdf5 format
-----------------------------------------------------------------------------------------
Input(s):
sys.argv[1]: subject name (e.g. 'sub-01')
sys.argv[2]: acquisition ('acq-2p5mm','acq-2mm')
sys.argv[3]: fit model (e.g. 'gauss')
-----------------------------------------------------------------------------------------
Output(s):
hdf5 files
-----------------------------------------------------------------------------------------
To run:
cd /home/mszinte/projects/pRFseqTest/mri_analysis/
python post_fit/roi_to_hdf5.py sub-01 acq-2p5mm gauss
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
from utils import set_pycortex_config_file, mask_nifti_2_hdf5

# Get inputs
# ----------
subject = sys.argv[1]
acq = sys.argv[2]
fit_model = sys.argv[3]

# Define analysis parameters
# --------------------------
with open('settings.json') as f:
    json_s = f.read()
    analysis_info = json.loads(json_s)

# Define folders and settings
# ---------------------------
base_dir = analysis_info['base_dir']
prf_signs = analysis_info['prf_signs']
cortex_dir = "{base_dir}/pp_data/cortex/db/{subject}".format(base_dir = base_dir, subject = subject)
rois_mask_dir = "{base_dir}/pp_data/{subject}/{fit_model}/masks/".format(base_dir = base_dir, subject = subject, fit_model = fit_model)
h5_dir = "{base_dir}/pp_data/{subject}/{fit_model}/h5/".format(base_dir = base_dir, subject = subject, fit_model = fit_model)
deriv_dir = opj(base_dir,'pp_data',subject,fit_model,'deriv')
xfm_name = "identity.{acq}".format(acq = acq)
rois = analysis_info['rois']
cortical_mask = analysis_info['cortical_mask']

# Set pycortex db and colormaps
# -----------------------------
set_pycortex_config_file(base_dir)

# Create ROI masks
# ----------------
ref_file = "{cortex_dir}/transforms/{xfm_name}/reference.nii.gz".format(cortex_dir = cortex_dir, xfm_name = xfm_name)
ref_img = nb.load(ref_file)
try: os.makedirs(rois_mask_dir)
except: pass

for roi in rois:
	print('creating {roi} {cortical_mask} mask'.format(roi = roi, cortical_mask = cortical_mask))
	roi_mask = cortex.utils.get_roi_masks(subject = subject, xfmname = xfm_name, gm_sampler = cortical_mask, roi_list = roi, return_dict = True)
	roi_mask_img = nb.Nifti1Image(dataobj = roi_mask[roi].transpose((2,1,0)), affine = ref_img.affine, header = ref_img.header)
	roi_mask_file = "{rois_mask_dir}/{roi}_{cortical_mask}_{acq}.nii.gz".format(rois_mask_dir = rois_mask_dir, roi = roi, cortical_mask = cortical_mask, acq = acq)
	roi_mask_img.to_filename(roi_mask_file)

# Create HDF5 files
# -----------------
try: os.makedirs(h5_dir)
except: pass

for roi in rois:
	print('creating {roi} h5 files'.format(roi = roi))

	h5_file = "{h5_dir}{roi}_{acq}.h5".format(h5_dir = h5_dir, roi = roi, acq = acq)
	
	try: os.system('rm {h5_file}'.format(h5_file = h5_file))
	except: pass

	mask_file = "{rois_mask_dir}/{roi}_{cortical_mask}_{acq}.nii.gz".format(rois_mask_dir = rois_mask_dir, roi = roi, cortical_mask = cortical_mask, acq = acq)
	
	for prf_sign in prf_signs:

		deriv_file = "{deriv_dir}/{prf_sign}/prf_deriv_{acq}_{prf_sign}.nii.gz".format(deriv_dir = deriv_dir,acq = acq, prf_sign = prf_sign)

		mask_nifti_2_hdf5(	in_file = deriv_file,
							mask_file = mask_file,
							hdf5_file = h5_file,
							folder_alias = prf_sign)

