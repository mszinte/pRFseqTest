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

# Functions import
# ----------------
from utils import set_pycortex_config_file, draw_cortex_vertex

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
base_dir = analysis_info['base_dir_local']
deriv_dir = opj(base_dir,'pp_data',subject,fit_model,'deriv')
cortex_dir = "{base_dir}/pp_data/cortex/db/{subject}".format(base_dir = base_dir, subject = subject)
fs_dir = "{base_dir}/deriv_data/fmriprep/freesurfer/".format(base_dir = base_dir)
fmriprep_dir = "{base_dir}/deriv_data/fmriprep/fmriprep/".format(base_dir = base_dir)
xfm_name = "AttendStim_{acq}".format(acq = acq)
xfm_dir = "{cortex_dir}/transforms/{xfm_name}".format(cortex_dir = cortex_dir, xfm_name = xfm_name)

# Set pycortex db and colormaps
# -----------------------------
set_pycortex_config_file(base_dir)

# Add participant to pycortex db
# ------------------------------

if os.path.isdir(cortex_dir) == False:
    print('add subject to pycortex database')

    cortex.freesurfer.import_subj(subject = subject, sname = subject, freesurfer_subject_dir = fs_dir)
    cortex.freesurfer.import_flat(subject = subject, patch='full', freesurfer_subject_dir = fs_dir)

# Create pycortex xfm
# -------------------
if os.path.isdir(xfm_dir) == False:

    deb()
    data_file  =  "{base_dir}/pp_data/{sub}/func/{sub}_task-AttendStim_{acq}_fmriprep_sg_psc_avg.nii.gz".format(
                                base_dir =base_dir,
                                sub = subject,
                                acq = acq)

    try:os.makedirs(opj(base_dir,'pp_data',subject,'reg'))
    except:pass

    # get mean file
    mean_file  =  "{base_dir}/pp_data/{sub}/reg/{sub}_task-AttendStim_{acq}_mean.nii.gz".format(
                                base_dir = base_dir, sub = subject, acq = acq)

    mean_cmd = fsl.maths.MeanImage( output_type = 'NIFTI_GZ', in_file = data_file,
                                    dimension = 'T', out_file = mean_file)    
    mean_cmd.run()

    # get reg mat
    reg_cmd = freesurfer.BBRegister(init = 'header', contrast_type = 't2', 
                                    subjects_dir = fs_dir, source_file = mean_file,
                                    subject_id = subject, out_fsl_file = True)
    reg_cmd.run()

    # create xfm transform
    xfm_file = "{mean_file}_bbreg_{subject}.mat".format(mean_file = mean_file[:-7],subject = subject)
    t1_file = "{fmriprep_dir}{subject}/anat/{subject}_desc-preproc_T1w.nii.gz".format(    
                                    fmriprep_dir = fmriprep_dir, subject = subject)

    xfm = cortex.xfm.Transform.from_fsl(xfm = xfm_file, func_nii = mean_file, anat_nii = t1_file)
    xfm.save(subject = subject, name = xfm_name)


# Draw pycortex flatmaps
# ----------------------
print('draw pycortex flatmaps')
cmap_neg_pos = 'RdBu_r'
cmap_polar = 'hsv'
col_offset = 1/14.0
polar_col_steps = [4.0, 8.0, 16.0, 255.0]
cmap_ecc_size = 'Spectral'
cmap_pos = 'Reds'
sign_idx, rsq_idx, ecc_idx, polar_real_idx, polar_imag_idx , size_idx, \
            non_lin_idx, amp_idx, baseline_idx, cov_idx, x_idx, y_idx = 0,1,2,3,4,5,6,7,8,9,10,11

for mask_dir in ['all','pos','neg']:

    # Create figure folders
    vertex_names = []
    all_vertex   = []
    exec('fig_roi_dir_{mask_dir} = opj(base_dir,"pp_data",subject,fit_model,"figs","roi","{mask_dir}")'.format(mask_dir=mask_dir))
    try: exec('os.makedirs(fig_roi_dir_{mask_dir})'.format(mask_dir=mask_dir))
    except: pass

    # Load data
    deriv_mat=[]
    deriv_matfn = "{deriv_dir}/{mask_dir}/prf_deriv_{acq}_{mask_dir}.nii.gz".format(deriv_dir = deriv_dir,acq = acq, mask_dir = mask_dir)
    img_deriv_mat = nb.load(deriv_matfn)
    deriv_mat = img_deriv_mat.get_fdata()

    # R-square
    rsq_data = deriv_mat[...,rsq_idx]
    alpha = rsq_data
    param_rsq = {'subject': 'fsaverage', 'data': rsq_data.T, 'cmap': cmap_pos, 'alpha': alpha.T, 'vmin': 0,'vmax': 1,'cbar': 'discrete'}
    vertex_names.append('rsq')

    # Polar angle
    pol_comp_num = deriv_mat[...,polar_real_idx] + 1j * deriv_mat[...,polar_imag_idx]
    polar_ang = np.angle(pol_comp_num)
    ang_norm = (polar_ang + np.pi) / (np.pi * 2.0)

    for cmap_steps in polar_col_steps:
        param_polar = {'data': ang_norm.T, 'cmap': cmap_polar, 'alpha': alpha.T, 'vmin': 0, 'vmax': 1, 'cmap_steps': cmap_steps,\
                       'curv_brightness': 0.05, 'curv_contrast': 0.1, 'cbar': 'polar', 'col_offset': col_offset}
        exec('param_polar_{csteps} = param_polar'.format(csteps = int(cmap_steps)))
        exec('vertex_names.append("polar_{csteps}")'.format(csteps = int(cmap_steps)))

    # Eccentricity
    ecc_data = deriv_mat[...,ecc_idx]
    param_ecc = {'data': ecc_data.T, 'cmap': cmap_ecc_size, 'alpha': alpha.T, 'vmin': 0, 'vmax': 8,'cbar': 'ecc'}
    vertex_names.append('ecc')

    # Sign
    sign_data = deriv_mat[...,sign_idx]
    param_sign = {'data': sign_data.T, 'cmap': cmap_neg_pos, 'alpha': alpha.T, 'vmin': -1, 'vmax': 1, 'cbar': 'discrete'}
    vertex_names.append('sign')

    # Size
    size_data = deriv_mat[...,size_idx]
    param_size = {'data': size_data.T, 'cmap': cmap_ecc_size, 'alpha': alpha.T, 'vmin': 0, 'vmax': 8, 'cbar': 'discrete'}
    vertex_names.append('size')

    # Amplitude
    amp_data = deriv_mat[...,amp_idx]
    param_amp = {'data': amp_data.T, 'cmap': cmap_neg_pos, 'alpha': alpha.T, 'vmin': -0.005, 'vmax': 0.005, 'cbar': 'discrete'}
    vertex_names.append('amp')

    # Baseline
    baseline_data = deriv_mat[...,baseline_idx]
    param_baseline = {'data': baseline_data.T, 'cmap': cmap_neg_pos, 'alpha': alpha.T, 'vmin': -10, 'vmax': 10,\
                      'curv_brightness': 0.05, 'curv_contrast': 0.1,'cbar': 'discrete'}
    vertex_names.append('baseline')

    # Non-linearity
    non_lin_data = deriv_mat[...,non_lin_idx]
    param_non_lin = {'data': non_lin_data.T, 'cmap': cmap_pos, 'alpha': alpha.T, 'vmin': 0, 'vmax': 1.5, 'cbar': 'discrete'}
    vertex_names.append('non_lin')

    # Coverage
    cov_data = deriv_mat[...,cov_idx]
    param_cov = {'data': cov_data.T, 'cmap': cmap_pos, 'alpha': alpha.T,'vmin': 0, 'vmax': 1, 'cbar': 'discrete'}
    vertex_names.append('cov')

    # Draw figures
    for vertex_name in vertex_names:
        roi_name = '{vertex_name}_{mask_dir}'.format(vertex_name = vertex_name, mask_dir = mask_dir)

        roi_param = {   'subject': subject,
                        'xfmname': xfm_name,
                        'add_roi': False,
                        'roi_name': roi_name}

        exec('param_{vertex_name}.update(roi_param)'.format(vertex_name = vertex_name))
        exec('vertex_rgb = draw_cortex_vertex(**param_{vertex_name})'.format(vertex_name=vertex_name))
        exec('plt.savefig(opj(fig_roi_dir_{mask_dir}, "{vertex_name}_{acq}_{mask_dir}.pdf"),facecolor="w")'.format(mask_dir=mask_dir,vertex_name = vertex_name,acq = acq))
        
        
    plt.close()