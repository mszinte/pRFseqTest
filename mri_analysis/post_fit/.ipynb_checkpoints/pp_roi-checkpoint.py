"""
-----------------------------------------------------------------------------------------
pp_roi.py
-----------------------------------------------------------------------------------------
Goal of the script:
Plot pycortex maps
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
python post_fit/post_fit.py sub-01 gauss 2500
python post_fit/post_fit.py sub-02 gauss 2500
> to run localy
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
import numpy as np
import ipdb
import platform
import matplotlib.pyplot as plt
opj = os.path.join
deb = ipdb.set_trace

# MRI imports
# -----------
import nibabel as nb
import cortex
from cortex.fmriprep import *
from nipype.interfaces import fsl, freesurfer
from nilearn import image

# Functions import
# ----------------
from utils import set_pycortex_config_file, draw_cortex_vertex

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

# Define folder
# -------------
base_dir = analysis_info['base_dir_local']
deriv_dir = opj(base_dir,'pp_data',subject,fit_model,'deriv')
cortex_dir = "{base_dir}/pp_data/cortex/db/{subject}".format(base_dir = base_dir, subject = subject)
fs_dir = "{base_dir}/deriv_data/fmriprep/freesurfer/".format(base_dir = base_dir)
fmriprep_dir = "{base_dir}/deriv_data/fmriprep/".format(base_dir = base_dir)
bids_dir = "{base_dir}/bids_data/".format(base_dir = base_dir)
xfm_name = 'identity.t1w'

# Set pycortex db and colormaps
# -----------------------------
set_pycortex_config_file(base_dir)

# Add participant to pycortex db
# ------------------------------
if os.path.isdir(cortex_dir) == False:
    cortex.fmriprep.import_subj(subject = subject[-2:], source_dir = fmriprep_dir, sname = subject)
    cortex.freesurfer.import_flat(subject = subject, patch = 'full', freesurfer_subject_dir = fs_dir, sname = subject)
    
    t1w = cortex.db.get_anat(subject)
    transform = cortex.xfm.Transform(np.identity(4), t1w)
    transform.save(subject, xfm_name, 'magnet')

# Draw pycortex flatmaps
# ----------------------
print('draw pycortex flatmaps')
sign_idx, rsq_idx, ecc_idx, polar_real_idx, polar_imag_idx , size_idx, \
            non_lin_idx, amp_idx, baseline_idx, cov_idx, x_idx, y_idx = 0,1,2,3,4,5,6,7,8,9,10,11

cmap_neg_pos = 'RdBu_r'
cmap_polar = 'Retinotopy_RYBCR'
col_offset = 14/16
polar_col_steps = [4.0, 8.0, 16.0, 255.0]

cmap_uni = 'Reds'
cmap_ecc_size = 'Spectral'

for mask_dir in ['all','pos','neg']:

    # Create figure folders
    maps_names = []
    exec('fig_roi_dir_{mask_dir} = opj(base_dir,"pp_data",subject,fit_model,"figs","flatmaps","{mask_dir}")'.format(mask_dir = mask_dir))
    try: exec('os.makedirs(fig_roi_dir_{mask_dir})'.format(mask_dir=mask_dir))
    except: pass

    # Load data
    deriv_mat=[]
    deriv_mat_file = "{deriv_dir}/{mask_dir}/prf_deriv_{acq}_{mask_dir}.nii.gz".format(deriv_dir = deriv_dir,acq = acq, mask_dir = mask_dir)
    img_deriv_mat = nb.load(deriv_mat_file)
    deriv_mat = img_deriv_mat.get_data()

    # interpolate to t1w
    t1w = cortex.db.get_anat(subject)
    deriv_rs = image.resample_to_img(source_img = img_deriv_mat, target_img = t1w, interpolation = 'nearest')
    deriv_mat_rs = deriv_rs.get_data()

    # R-square
    rsq_data = deriv_mat_rs[...,rsq_idx]
    alpha = rsq_data
    param_rsq = {'data': rsq_data, 'cmap': cmap_uni, 'alpha': alpha, 'vmin': 0,'vmax': 1,'cbar': 'discrete', 'description': 'pRF rsquare'}
    maps_names.append('rsq')

    # Polar angle
    pol_comp_num = deriv_mat_rs[...,polar_real_idx] + 1j * deriv_mat_rs[...,polar_imag_idx]
    polar_ang = np.angle(pol_comp_num)
    ang_norm = (polar_ang + np.pi) / (np.pi * 2.0)
    ang_norm = np.fmod(ang_norm + col_offset,1)

    for cmap_steps in polar_col_steps:
        param_polar = {'data': ang_norm, 'cmap': cmap_polar, 'alpha': alpha, 'vmin': 0, 'vmax': 1, 'cmap_steps': cmap_steps,
                       'cbar': 'polar', 'col_offset': col_offset, 'description': 'pRF polar:{:3.0f} steps'.format(cmap_steps)}
        exec('param_polar_{csteps} = param_polar'.format(csteps = int(cmap_steps)))
        exec('maps_names.append("polar_{csteps}")'.format(csteps = int(cmap_steps)))

    # Eccentricity
    ecc_data = deriv_mat_rs[...,ecc_idx]
    param_ecc = {'data': ecc_data, 'cmap': cmap_ecc_size, 'alpha': alpha, 'vmin': 0, 'vmax': 15,'cbar': 'ecc', 'description': 'pRF eccentricity'}
    maps_names.append('ecc')

    # Sign
    sign_data = deriv_mat_rs[...,sign_idx]
    param_sign = {'data': sign_data, 'cmap': cmap_neg_pos, 'alpha': alpha, 'vmin': -1, 'vmax': 1, 'cbar': 'discrete', 'description': 'pRF sign'}
    maps_names.append('sign')

    # Size
    size_data = deriv_mat_rs[...,size_idx]
    param_size = {'data': size_data, 'cmap': cmap_ecc_size, 'alpha': alpha, 'vmin': 0, 'vmax': 8, 'cbar': 'discrete', 'description': 'pRF size'}
    maps_names.append('size')

    # Coverage
    cov_data = deriv_mat_rs[...,cov_idx]
    param_cov = {'data': cov_data, 'cmap': cmap_uni, 'alpha': alpha,'vmin': 0, 'vmax': 1, 'cbar': 'discrete', 'description': 'pRF coverage'}
    maps_names.append('cov')

    # Draw figures
    volumes = {}
    for maps_name in maps_names:
        roi_name = '{maps_name}_{mask_dir}'.format(maps_name = maps_name, mask_dir = mask_dir)

        roi_param = {   'subject': subject,
                        'xfmname': xfm_name,
                        'add_roi': False, 
                        'roi_name': roi_name,
                        'curv_brightness': 0.7, 
                        'curv_contrast': 0.3}

        exec('param_{maps_name}.update(roi_param)'.format(maps_name = maps_name))
        exec('volume_{maps_name} = draw_cortex_vertex(**param_{maps_name})'.format(maps_name=maps_name))
        exec('plt.savefig(opj(fig_roi_dir_{mask_dir}, "{maps_name}_{acq}_{mask_dir}.pdf"),facecolor="w")'.format(mask_dir=mask_dir,maps_name = maps_name, acq = acq))
        plt.close()    
        print('new_volume = { param_{maps_name}["description"]: volume_{maps_name}}'.format(maps_name=maps_name))
        deb()
        volumes.update(new_volume)
        
        
    
        
    'First Dataset': volume1,
    'Second Dataset': volume2,
    'Third Dataset': volume3,
}

# create viewer
cortex.webgl.show(data=volumes)
    