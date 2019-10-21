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

# Functions import
# ----------------
from utils import draw_cortex_vertex

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

for mask_dir in ['all','pos','neg']:

    print('draw pycortex flatmaps: {mask_dir}'.format(mask_dir = mask_dir))
    # Create figure folders
    maps_names = []
    exec('flatmaps_dir_{mask_dir} = opj(base_dir,"pp_data",subject,fit_model,"pycortex_outputs","flatmaps","{mask_dir}")'.format(mask_dir = mask_dir))
    exec('dataset_dir_{mask_dir} = opj(base_dir,"pp_data",subject,fit_model,"pycortex_outputs","dataset","{mask_dir}")'.format(mask_dir = mask_dir))
    exec('webviewer_dir_{mask_dir} = opj(base_dir,"pp_data",subject,fit_model,"pycortex_outputs","webviewer","{mask_dir}")'.format(mask_dir = mask_dir))
    try: 
        exec('os.makedirs(flatmaps_dir_{mask_dir})'.format(mask_dir=mask_dir))
        exec('os.makedirs(dataset_dir_{mask_dir})'.format(mask_dir=mask_dir))
        exec('os.makedirs(webviewer_dir_{mask_dir})'.format(mask_dir=mask_dir))
    except: 
        pass

    # Load data
    deriv_mat_rs_file = "{deriv_dir}/{mask_dir}/prf_deriv_{acq}_{mask_dir}_rs.nii.gz".format(deriv_dir = deriv_dir,acq = acq, mask_dir = mask_dir)
    img_deriv_mat_rs = nb.load(deriv_mat_rs_file)
    deriv_mat_rs = img_deriv_mat_rs.get_data()
    
    # R-square
    rsq_data = deriv_mat_rs[...,rsq_idx]
    alpha = rsq_data
    param_rsq = {'data': rsq_data, 'cmap': cmap_uni, 'alpha': alpha, 'vmin': 0,'vmax': 1,'cbar': 'discrete', 'description': 'pRF rsquare', 'curv_brightness': 1, 'curv_contrast': 0.1}
    maps_names.append('rsq')

    # Polar angle
    pol_comp_num = deriv_mat_rs[...,polar_real_idx] + 1j * deriv_mat_rs[...,polar_imag_idx]
    polar_ang = np.angle(pol_comp_num)
    ang_norm = (polar_ang + np.pi) / (np.pi * 2.0)
    ang_norm = np.fmod(ang_norm + col_offset,1)

    for cmap_steps in polar_col_steps:
        param_polar = {'data': ang_norm, 'cmap': cmap_polar, 'alpha': alpha, 'vmin': 0, 'vmax': 1, 'cmap_steps': cmap_steps,
                       'cbar': 'polar', 'col_offset': col_offset, 'description': 'pRF polar:{:3.0f} steps'.format(cmap_steps), 'curv_brightness': 0.1, 'curv_contrast': 0.25}
        exec('param_polar_{csteps} = param_polar'.format(csteps = int(cmap_steps)))
        exec('maps_names.append("polar_{csteps}")'.format(csteps = int(cmap_steps)))

    # Eccentricity
    ecc_data = deriv_mat_rs[...,ecc_idx]
    param_ecc = {'data': ecc_data, 'cmap': cmap_ecc_size, 'alpha': alpha, 'vmin': 0, 'vmax': 15,'cbar': 'ecc', 'description': 'pRF eccentricity', 'curv_brightness': 1, 'curv_contrast': 0.1}
    maps_names.append('ecc')

    # Sign
    sign_data = deriv_mat_rs[...,sign_idx]
    param_sign = {'data': sign_data, 'cmap': cmap_neg_pos, 'alpha': alpha, 'vmin': -1, 'vmax': 1, 'cbar': 'discrete', 'description': 'pRF sign', 'curv_brightness': 1, 'curv_contrast': 0.1}
    maps_names.append('sign')

    # Size
    size_data = deriv_mat_rs[...,size_idx]
    param_size = {'data': size_data, 'cmap': cmap_ecc_size, 'alpha': alpha, 'vmin': 0, 'vmax': 8, 'cbar': 'discrete', 'description': 'pRF size', 'curv_brightness': 1, 'curv_contrast': 0.1}
    maps_names.append('size')

    # Coverage
    cov_data = deriv_mat_rs[...,cov_idx]
    param_cov = {'data': cov_data, 'cmap': cmap_uni, 'alpha': alpha,'vmin': 0, 'vmax': 1, 'cbar': 'discrete', 'description': 'pRF coverage', 'curv_brightness': 1, 'curv_contrast': 0.1}
    maps_names.append('cov')

    # Draw flatmaps
    volumes = {}
    for maps_name in maps_names:
        roi_name = '{maps_name}_{mask_dir}'.format(maps_name = maps_name, mask_dir = mask_dir)

        roi_param = {'subject': subject, 'xfmname': xfm_name, 'add_roi': False, 'roi_name': roi_name}

        exec('param_{maps_name}.update(roi_param)'.format(maps_name = maps_name))
        exec('volume_{maps_name} = draw_cortex_vertex(**param_{maps_name})'.format(maps_name=maps_name))
        exec('plt.savefig(opj(flatmaps_dir_{mask_dir}, "{maps_name}_{acq}_{mask_dir}.pdf"),facecolor="w")'.format(mask_dir = mask_dir,maps_name = maps_name, acq = acq))
        plt.close()
        exec('vol_description = param_{maps_name}["description"]'.format(maps_name = maps_name))
        exec('volume = volume_{maps_name}'.format(maps_name = maps_name))
        volumes.update({vol_description:volume})
    
    print('save pycortex dataset: {mask_dir}'.format(mask_dir = mask_dir))
    exec('dataset_{mask_dir} = cortex.Dataset(data = volumes)'.format(mask_dir = mask_dir))
    exec('dataset_{mask_dir}.save(opj(dataset_dir_{mask_dir}, "{acq}_{mask_dir}.hdf"))'.format(acq = acq, mask_dir = mask_dir))
    
    print('save pycortex webviewer: {mask_dir}'.format(mask_dir = mask_dir))
    exec('cortex.webgl.make_static(outpath = webviewer_dir_{mask_dir}, data = volumes)'.format(acq = acq, mask_dir = mask_dir))