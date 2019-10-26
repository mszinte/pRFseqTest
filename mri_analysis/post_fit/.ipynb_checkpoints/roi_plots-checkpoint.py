"""
-----------------------------------------------------------------------------------------
roi_plots.py
-----------------------------------------------------------------------------------------
Goal of the script:
Draw roi plots (maps, ecc vs. params, laterality, time course)
-----------------------------------------------------------------------------------------
Input(s):
sys.argv[1]: subject number
sys.argv[2]: acquisition ('acq-2p5mm','acq-2mm')
sys.argv[3]: fit model ('gauss','css')
sys.argv[4]: save svg (1 = YES, 0 = NO)
-----------------------------------------------------------------------------------------
Output(s):
None
-----------------------------------------------------------------------------------------
To run:
cd /home/mszinte/projects/pRFseqTest/mri_analysis/
python post_fit/roi_plots.py sub-01 acq-2p5mm gauss 0
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
import h5py
import scipy.io
opj = os.path.join
deb = ipdb.set_trace

# MRI imports
# -----------
import nibabel as nb
import cortex

# Functions import
# ----------------
from plot_class import PlotOperator
from utils import set_pycortex_config_file

# Bokeh imports
# ---------------
from bokeh.io import output_notebook, show, save, output_file, export_png, export_svgs
from bokeh.layouts import row, column, gridplot

# Popeye imports
import popeye.utilities as utils
from popeye.visual_stimulus import VisualStimulus
import popeye.css_neg as css
import popeye.og_neg as og

# Check system
# ------------
sys.exit('Drawing Flatmaps only works with Python 3. Aborting.') if sys.version_info[0] > 3 else None

# Get inputs
# ----------
subject = sys.argv[1]
acq = sys.argv[2]
fit_model = sys.argv[3]
save_svg = int(sys.argv[4])

# Define analysis parameters
# --------------------------
with open('settings.json') as f:
    json_s = f.read()
    analysis_info = json.loads(json_s)

# Define folders and settings
# ---------------------------
base_dir = analysis_info['base_dir']
prf_signs = analysis_info['prf_signs']
rois = analysis_info['rois']
deriv_dir = "{base_dir}/pp_data/{subject}/{fit_model}/deriv".format(base_dir = base_dir, subject = subject, fit_model = fit_model)
h5_dir = "{base_dir}/pp_data/{subject}/{fit_model}/h5".format(base_dir = base_dir, subject = subject, fit_model = fit_model)
bokeh_dir = "{base_dir}/pp_data/{subject}/{fit_model}/bokeh_outputs".format(base_dir = base_dir, subject = subject, fit_model = fit_model)
try: os.makedirs(bokeh_dir)
except: pass

# Create stimulus design
# ----------------------
visual_dm_file = scipy.io.loadmat("{base_dir}/pp_data/visual_dm/vis_design.mat".format(base_dir = base_dir))
visual_dm = visual_dm_file['stim']
stimulus = VisualStimulus(  stim_arr = visual_dm,
                            viewing_distance = analysis_info["screen_distance"],
                            screen_width = analysis_info["screen_width"],
                            scale_factor = 1/10.0,
                            tr_length = analysis_info["TR"],
                            dtype = np.short)

if fit_model == 'gauss':
    fit_func = og.GaussianFit
    model_func = og.GaussianModel(  stimulus = stimulus,
                                    hrf_model = utils.spm_hrf)
    step_r2 = [0,100/3.0,250/3.0,100]
    list_r2_level = ['High','Low']
    step_params = [0,100/3.0,200/3.0,100]
    list_params_level = ['High','Low']
    list_params = ['ecc','amp','size','cov']
    num_plots = len(list_params)*len(step_params)*(len(step_r2)-1)
    
elif fit_model == 'css':
    fit_func = css.CompressiveSpatialSummationFit
    model_func = css.CompressiveSpatialSummationModel(  stimulus = stimulus,
                                                        hrf_model = utils.spm_hrf)
    step_r2 = [0,100/3.0,200/3.0,100]
    list_r2_level = ['Low','High']
    step_params = [0,100/3.0,200/3.0,100]
    list_params_level = ['Low','High']
    list_params = ['ecc','amp','size','cov','non_lin']

    num_plots = len(list_params)*len(step_params)*(len(step_r2)-1)
    
model_func.hrf_delay = 0

# Draw main analysis figure
# -------------------------
print('creating bokeh plots')
sign_idx, rsq_idx, ecc_idx, polar_real_idx, polar_imag_idx , size_idx, \
            non_lin_idx, amp_idx, baseline_idx, cov_idx, x_idx, y_idx = 0,1,2,3,4,5,6,7,8,9,10,11

# Draw main analysis figure
# -------------------------
print('creating bokeh plots')
sign_idx, rsq_idx, ecc_idx, polar_real_idx, polar_imag_idx , size_idx, \
            non_lin_idx, amp_idx, baseline_idx, cov_idx, x_idx, y_idx = 0,1,2,3,4,5,6,7,8,9,10,11


# Initialize data
# ---------------
for roi_num, roi in enumerate(rois):
    for prf_sign in prf_signs:
        
        tc_mat = []

        # create html folder
        exec("html_dir_{prf_sign} = '{bokeh_dir}/html/{prf_sign}'".format(bokeh_dir = bokeh_dir, prf_sign = prf_sign))
        try: exec('os.makedirs(html_dir_{prf_sign})'.format(prf_sign = prf_sign))
        except: pass

        # create svg folder
        if save_svg == 1:
            exec("svg_dir_{prf_sign} = '{bokeh_dir}/svg/{prf_sign}'".format(bokeh_dir = bokeh_dir, prf_sign = prf_sign))
            try: os.makedirs('os.makedirs(svg_dir_{prf_sign})'.format(prf_sign = prf_sign))
            except: pass
            exec('svg_folder = svg_dir_{prf_sign}'.format(prf_sign = prf_sign))

        # load h5 file
        h5_file = h5py.File("{h5_dir}/{roi}_{acq}.h5".format(h5_dir = h5_dir, roi = roi, acq = acq))

        # load deriv data
        deriv_data_name = "prf_deriv_{acq}_{prf_sign}".format(acq = acq, prf_sign = prf_sign)
        deriv_data = h5_file['{folder_alias}/{data_name}'.format(folder_alias = prf_sign, data_name = deriv_data_name)]
        
        # load time course data
        tc_data_name = "{subject}_task-AttendStim_{acq}_fmriprep_sg_psc_avg".format(subject = subject, acq = acq)
        tc_data = h5_file['{folder_alias}/{data_name}'.format(folder_alias = prf_sign, data_name = tc_data_name)]

        # load coordinates data
        coord_data_name = "coord"
        coord_data = h5_file['{folder_alias}/{data_name}'.format(folder_alias = prf_sign, data_name = coord_data_name)]
        
        # threshold data

        voxel_num_ini = deriv_data.shape[0]
        voxel_num = 0
        if voxel_num_ini > 0:
            # take out nan
            deriv_data = deriv_data[~np.isnan(deriv_data[:,rsq_idx]),:]
            tc_data = tc_data[~np.isnan(deriv_data[:,rsq_idx]),:]
            coord_data = coord_data[~np.isnan(deriv_data[:,rsq_idx]),:]
            
            # threshold on eccentricity and size
            deriv_data_th = deriv_data
            size_th_down = deriv_data_th[:,size_idx] >= analysis_info['size_th'][0]
            size_th_up = deriv_data_th[:,size_idx] <= analysis_info['size_th'][1]
            ecc_th_down = deriv_data_th[:,ecc_idx] >= analysis_info['ecc_th'][0]
            ecc_th_up = deriv_data_th[:,ecc_idx] <= analysis_info['ecc_th'][1]
            all_th = np.array((size_th_down,size_th_up,ecc_th_down,ecc_th_up)) 

            deriv_data = deriv_data[np.logical_and.reduce(all_th),:]
            tc_data = tc_data[np.logical_and.reduce(all_th),:]
            coord_data = coord_data[np.logical_and.reduce(all_th),:]
            voxel_num = deriv_data.shape[0]
        
        # Draw and save figures
        # ---------------------
        if voxel_num > 0:
            # randomize order and take 10 %
            new_order = np.random.permutation(voxel_num)
            deriv_data = deriv_data[new_order,:]
            tc_data = tc_data[new_order,:]
            coord_data = coord_data[new_order,:]
            
            print("drawing {roi}_{prf_sign} figures, n = {voxel_num}".format(roi = roi, voxel_num = voxel_num, prf_sign = prf_sign)) 

            data_source = { 'sign': deriv_data[:,sign_idx],
                            'rsq': deriv_data[:,rsq_idx],
                            'ecc': deriv_data[:,ecc_idx],
                            'sigma': deriv_data[:,size_idx],
                            'non_lin': deriv_data[:,non_lin_idx],
                            'beta': deriv_data[:,amp_idx],
                            'baseline': deriv_data[:,baseline_idx],
                            'cov': deriv_data[:,cov_idx],
                            'x': deriv_data[:,x_idx],
                            'y': deriv_data[:,y_idx],
                            'colors_ref': deriv_data[:,ecc_idx]}

            param_all = {   'roi_t': roi, 
                            'p_width': 400, 
                            'p_height': 400, 
                            'min_border_large': 10, 
                            'min_border_small': 5,
                            'bg_color': tuple([229,229,229]), 
                            'stim_color': tuple([250,250,250]), 
                            'hist_fill_color': tuple([255,255,255]),
                            'hist_line_color': tuple([0,0,0]), 
                            'stim_radius': analysis_info['stim_width']/2,
                            'cmap': 'Spectral',
                            'cmap_steps': 8,
                            'col_offset': 0,
                            'vmin': 0,
                            'vmax': 8,
                            'leg_xy_max_ratio': 1.8,
                            'dataMat': deriv_data,
                            'data_source': data_source,
                            'hist_range': (0,0.5),
                            'hist_steps': 0.5,
                            'h_hist_bins': 16,
                            'link_x': False,
                            'link_y': False,
                            'save_svg': save_svg
                            }

            if save_svg == 1:
                param_all.update({'svg_folder': svg_folder})

            plotter = PlotOperator(**param_all)
            
            # pRFmap
            title = '{roi} {prf_sign}: pRF map (n = {voxel_num})'.format(roi = roi, voxel_num = voxel_num, prf_sign = prf_sign)
            params_pRFmap = param_all

            params_pRFmap.update({   
                            'x_range': (-15, 15),
                            'y_range': (-15, 15),
                            'x_label': 'Horizontal coordinate (dva)',
                            'y_label': 'Vertical coordinate (dva)',
                            'x_source_label': 'x',
                            'y_source_label': 'y',
                            'x_tick_steps': 5,
                            'y_tick_steps': 5,
                            'v_hist_bins': 12,
                            'h_hist_bins': 12,
                            'main_fig_title': title})
            if save_svg == 1:
                params_pRFmap.update({'svg_filename': '{roi}_{prf_sign}_pRFmap'.format(prf_sign = prf_sign, roi= roi)})

            f_pRFmap, old_main_fig1 = plotter.draw_figure(  parameters = params_pRFmap, 
                                                            plot = 'map')

            
            # pRF cov
            title = '{roi} {prf_sign}: pRF density map'.format(roi = roi, prf_sign = prf_sign)
            params_pRFcov = param_all
            params_pRFcov.update({   
                            'x_range': (-15, 15), 
                            'y_range': (-15, 15),
                            'x_label': 'Horizontal coordinate (dva)',
                            'y_label': 'Vertical coordinate (dva)',
                            'x_tick_steps': 5,
                            'y_tick_steps': 5,
                            'smooth_factor': 15,
                            'cmap': 'viridis',
                            'cmap_steps': 10,
                            'col_offset': 0,
                            'vmin': 0,
                            'vmax': 1,
                            'cb_tick_steps': 0.2,
                            'condition': 'cov',
                            'cb_label': 'pRF coverage (norm.)',
                            'link_x': True,
                            'link_y': True})
            if save_svg == 1:
                params_pRFcov.update({'svg_filename': '{roi}_{prf_sign}_pRFcov'.format(prf_sign = prf_sign, roi= roi)})
                        
            params_pRFcov.update({'main_fig_title':   title})

            f_pRFcov, old_main_fig1 = plotter.draw_figure(parameters = params_pRFcov, 
                                           plot = 'cov',
                                           old_main_fig = old_main_fig1)

            all_fig1 = gridplot([ [f_pRFmap,f_pRFcov]], toolbar_location = 'right')

            exec('output_file_html = opj(html_dir_{prf_sign},"{roi}_{prf_sign}_pRFmap.html")'.format(prf_sign = prf_sign,roi = roi))
            html_title = "Subject: {subject} | ROI: {roi} | Sign: {prf_sign} | Voxel num: {voxel_num} | Figures: pRF maps parameters and density".format(
                                                subject = subject,roi = roi,prf_sign = prf_sign,voxel_num = voxel_num)
            output_file(output_file_html, title = html_title)
            save(all_fig1)
        else:
            print("drawing {roi}_{prf_sign} figures not possible: n = {voxel_num}".format(roi = roi,voxel_num = voxel_num,prf_sign = prf_sign)) 

            # pRFecc
            # old_main_fig = []
                    #     f_pRFecc = []
                    #     if fit_model == 'gauss':
                    #         type_comp_list = ['Size','R2','Amplitude','Coverage','Baseline']
                    #     elif fit_model == 'css':
                    #         type_comp_list = ['Size','R2','Non-Linearity','Amplitude','Coverage','Baseline']

                    #     for numData, type_comp in enumerate(type_comp_list):

                    #         params_pRFecc = param_all
                    #         params_pRFecc.update(   
                    #                    {    'x_range':          (0, 10),
                    #                         'x_label':          'Eccentricity (dva)',
                    #                         'x_tick_steps':     2,
                    #                         'x_source_label':   'ecc',
                    #                         'draw_reg':         False,
                    #                         'link_x':           True})
                    #         if save_svg == 1:
                    #             params_pRFecc.update(
                    #                         {   'svg_subfolder':    '{roi}_{hemi}_{prf_sign}_pRFecc'.format(prf_sign=prf_sign,hemi=hemi,roi= roi_text)})

                    #         if type_comp == 'Size':
                    #             params_pRFecc.update(
                    #                         {   'y_range':          (0, 10),
                    #                             'y_label':          'Size (dva)',
                    #                             'y_source_label':   'sigma',
                    #                             'y_tick_steps':     2,
                    #                             'v_hist_bins':      20,
                    #                             'draw_reg':         True})
                    #             if save_svg == 1:
                    #                 params_pRFecc.update(
                    #                         {   'svg_filename':     '{roi}_{hemi}_{prf_sign}_pRFecc_size'.format(prf_sign=prf_sign,hemi=hemi,roi= roi_text)})

                    #         elif type_comp == 'R2':
                    #             params_pRFecc.update(
                    #                         {   'y_range':          (0, 1),
                    #                             'y_label':          'R2',
                    #                             'y_source_label':   'rsq',
                    #                             'y_tick_steps':     0.2,
                    #                             'v_hist_bins':      15})
                    #             if save_svg == 1:
                    #                 params_pRFecc.update(
                    #                         {   'svg_filename':     '{roi}_{hemi}_{prf_sign}_pRFecc_rsq'.format(prf_sign=prf_sign,hemi=hemi,roi= roi_text)})

                    #         elif type_comp == 'Non-Linearity':
                    #             params_pRFecc.update(
                    #                         {   'y_range':          (0, 1.5),
                    #                             'y_label':          'Non-linearity',
                    #                             'y_source_label':   'non_lin',
                    #                             'y_tick_steps':     0.25,
                    #                             'v_hist_bins':      18})
                    #             if save_svg == 1:
                    #                 params_pRFecc.update(
                    #                         {   'svg_filename':     '{roi}_{hemi}_{prf_sign}_pRFecc_non_lin'.format(prf_sign=prf_sign,hemi=hemi,roi= roi_text)})

                    #         elif type_comp == 'Amplitude':
                    #             params_pRFecc.update(
                    #                         {   'y_range':          (-0.02, 0.02),
                    #                             'y_label':          'Amplitude',
                    #                             "y_source_label":   'beta',
                    #                             'y_tick_steps':     0.01,
                    #                             'v_hist_bins':      16})
                    #             if save_svg == 1:
                    #                 params_pRFecc.update(
                    #                         {   'svg_filename':     '{roi}_{hemi}_{prf_sign}_pRFecc_amp'.format(prf_sign=prf_sign,hemi=hemi,roi= roi_text)})

                    #         elif type_comp == 'Coverage':
                    #             params_pRFecc.update(
                    #                         {   'y_range':          (0, 1),
                    #                             'y_label':          'pRF coverage (%)',
                    #                             'y_source_label':   'cov',
                    #                             'y_tick_steps':     0.2,
                    #                             'v_hist_bins':      15})
                    #             if save_svg == 1:
                    #                 params_pRFecc.update(
                    #                         {   'svg_filename':     '{roi}_{hemi}_{prf_sign}_pRFecc_cov'.format(prf_sign=prf_sign,hemi=hemi,roi= roi_text)})

                    #         elif type_comp == 'Baseline':
                    #             params_pRFecc.update(
                    #                         {   'y_range':          (-200, 200),
                    #                             'y_label':          'Baseline',
                    #                             'y_source_label':   'baseline',
                    #                             'y_tick_steps':     25,
                    #                             'v_hist_bins':      16})
                    #             if save_svg == 1:
                    #                 params_pRFecc.update(
                    #                         {   'svg_filename':     '{roi}_{hemi}_{prf_sign}_pRFecc_baseline'.format(prf_sign=prf_sign,hemi=hemi,roi= roi_text)})

                    #         title = '{roi}{hemi} {prf_sign}: Eccentricity vs. {type_comp}'.format(roi = roi_text, hemi = hemi, type_comp = type_comp, prf_sign = prf_sign)
                    #         params_pRFecc.update({'main_fig_title':   title})

                    #         out1, old_main_fig  = plotter.draw_figure(  parameters = params_pRFecc,
                    #                                                     plot = 'ecc',
                    #                                                     old_main_fig = old_main_fig)
                    #         f_pRFecc.append(out1)


                    
                    #     # pRF lat
                    #     # -------
                    #     params_pRFlat = param_all
                    #     params_pRFlat.update(
                    #                 {   'p_width':          500, 
                    #                     'p_height':         500,
                    #                     'dataMat':          data,
                    #                     'x_range':          (-2.6, 2.6), 
                    #                     'y_range':          (-2.6, 2.6),
                    #                     'vmin':             0,
                    #                     'vmax':             0.2,
                    #                     'weighted_data':    True,
                    #                     'main_fig_title':   title,
                    #                     'cmap':             'hsv',
                    #                     'cmap_steps':       16,
                    #                     'ang_bins':         36,
                    #                     'hemi':             hemi})
                    #     if save_svg == 1:
                    #         params_pRFlat.update(
                    #                 {   'svg_filename':     '{roi}_{hemi}_{prf_sign}_pRFlat'.format(prf_sign=prf_sign,hemi=hemi,roi= roi_text)})
                    #     title = '{roi}{hemi} {prf_sign}: pRF laterality histogram'.format(roi = roi_text, hemi = hemi, type_comp = type_comp, prf_sign = prf_sign)
                    #     params_pRFcov.update({'main_fig_title':   title})
                    #     f_pRFlat = plotter.draw_figure(parameters = params_pRFlat, plot = 'lat')

                    #     # pRF tc
                    #     # -------
                    #     params_pRFtc = param_all
                    #     params_pRFtc.update(
                    #                 {   'p_width':          500, 
                    #                     'p_height':         500,
                    #                     'x_range_map':      (-12,12),
                    #                     'y_range_map':      (-12,12),
                    #                     'x_label_map':      'Horizontal coord. (dva)',
                    #                     'y_label_map':      'Vertical coord. (dva)',
                    #                     'x_tick_map':       4,
                    #                     'y_tick_map':       4,
                    #                     'x_range_tc':       (0,analysis_info['n_timepoints_per_run']),
                    #                     'x_label_tc':       'Time (s)',
                    #                     'y_label_tc':       'BOLD signal change (%)',
                    #                     'x_tick_tc':        50,
                    #                     'tr_dur':           analysis_info['TR'],
                    #                     'model_line_color': tuple([254,51,10]),
                    #                     'model_fill_color': tuple([254,51,10])
                    #                 })

                    #     if save_svg == 1:
                    #         params_pRFtc.update(
                    #                     {   'svg_subfolder':    '{roi}_{hemi}_{prf_sign}_pRFtc'.format(prf_sign=prf_sign,hemi=hemi,roi= roi_text)})

                    #     f_pRFtc = []

                    #     # get index matrices
                    #     prct_r2 = np.nanpercentile(a = data[:,rsq_idx],q = step_r2)
                    #     idx_r2 = []
                    #     for r2_step in np.arange(0,len(step_r2)-1,2):

                    #         idx_r2.append(np.logical_and(data[:,rsq_idx]>=prct_r2[r2_step],data[:,rsq_idx]<=prct_r2[r2_step+1]))


                    #     # get vertex number to draw
                    #     for param in list_params:
                    #         exec('prct_{param} = np.nanpercentile(a = data[:,{param}_idx],q = step_params)'.format(param = param))
                    #         exec('idx_{param} = []'.format(param = param))

                    #         for param_step in np.arange(0,len(step_params)-1,2):
                    #             exec('idx_{param}.append(np.logical_and(data[:,{param}_idx]>=prct_{param}[param_step],data[:,{param}_idx]<=prct_{param}[param_step+1]))'.format(param = param))

                    #         exec('num_{param} = []'.format(param = param))
                    #         exec('num_{param} = []'.format(param = param))
                    #         for r2_step in np.arange(0,2,1):

                    #             for param_step in np.arange(0,2,1):

                    #                 exec('mat = np.where(np.logical_and(idx_r2[r2_step],idx_{param}[param_step]))'.format(param = param))
                    #                 if mat[0].size == 0:
                    #                     exec('num_{param}.append(-1)'.format(param = param))
                    #                 else:
                    #                     exec('num_{param}.append(mat[0][np.random.choice(len(mat[0]))])'.format(param = param))

                    #     for param in list_params:
                    #         exec('num_vertex = num_{param}'.format(param = param))

                    #         for r2_level in list_r2_level:
                    #             if r2_level == 'Low':
                    #                 num_vertex2draw = num_vertex[0:2]
                    #             elif r2_level == 'High':
                    #                 num_vertex2draw = num_vertex[2:4]

                    #             params_pRFtc.update(   
                    #                             {   'params':               param,
                    #                                 'r2_level':             r2_level,
                    #                                 'deriv_mat':            data,
                    #                                 'tc_mat':               tc_mat,
                    #                                 'num_vertex':           num_vertex2draw,
                    #                                 'fit_model':            fit_model,
                    #                                 'model_func':           model_func,
                    #                                 'prf_sign':             prf_sign,
                    #                                 'title':                '{roi} {hemi} {sign}'.format(sign = prf_sign,hemi=hemi,roi = roi_text)
                    #                             })
                    #             if save_svg == 1:
                    #                 params_pRFtc.update(
                    #                                 {   'svg_filename':     '{roi}_{hemi}_{prf_sign}_pRFtc_rsq_{r2_level}_{param}'.format(prf_sign=prf_sign,hemi=hemi,roi= roi_text,r2_level = r2_level, param = param)})

                    #             out4,main_fig4  = plotter.draw_figure(    parameters =    params_pRFtc,
                    #                                             plot =          'tc')

                    #             f_pRFtc.append(out4)

                    #     # save files
                    #     if fit_model == 'gauss':
                    #         all_f1 = gridplot([ [f_pRFecc[0],f_pRFecc[1],f_pRFecc[2]],
                    #                             [f_pRFecc[3],f_pRFecc[4],None]],toolbar_location='right')


                    #         all_f4 = gridplot([ [f_pRFtc[0],f_pRFtc[1]],
                    #                             [f_pRFtc[2],f_pRFtc[3]],
                    #                             [f_pRFtc[4],f_pRFtc[5]],
                    #                             [f_pRFtc[6],f_pRFtc[7]]],toolbar_location = None)
                    #     elif fit_model == 'css':
                    #         all_f1 = gridplot([ [f_pRFecc[0],f_pRFecc[1],f_pRFecc[2]],
                    #                             [f_pRFecc[3],f_pRFecc[4],f_pRFecc[5]]],toolbar_location='right')
                    #         all_f4 = gridplot([ [f_pRFtc[0],f_pRFtc[1]],
                    #                             [f_pRFtc[2],f_pRFtc[3]],
                    #                             [f_pRFtc[2],f_pRFtc[3]],
                    #                             [f_pRFtc[4],f_pRFtc[5]]])

                    #     exec('output_file_html = opj(fig_bokeh_dir_{prf_sign}_{hemi},"{roi_text}_{hemi}_{prf_sign}_pRFecc.html")'.format(prf_sign = prf_sign,roi_text = roi_text, hemi = hemi))
                    #     output_file(output_file_html, title='Subject: %s | ROI: %s | Hemisphere: %s | Sign: %s | Vertex: %i | Figures: pRF parameters and density'%(subject,roi_text,hemi,prf_sign,vertex))
                    #     save(all_f1)

                    #     all_f2 = gridplot([ [f_pRFmap,f_pRFcov[0]]],toolbar_location='right')
                    #     exec('output_file_html = opj(fig_bokeh_dir_{prf_sign}_{hemi},"{roi_text}_{hemi}_{prf_sign}_pRFmap.html")'.format(prf_sign = prf_sign,roi_text = roi_text, hemi = hemi))
                    #     output_file(output_file_html, title='Subject: %s | ROI: %s | Hemisphere: %s | Sign: %s | Vertex: %i | Figures:pRF maps parameters and density'%(subject,roi_text,hemi,prf_sign,vertex))
                    #     save(all_f2)

                    #     all_f3 = gridplot([[f_pRFlat[0]]],toolbar_location = 'right')
                    #     exec('output_file_html = opj(fig_bokeh_dir_{prf_sign}_{hemi},"{roi_text}_{hemi}_{prf_sign}_pRFlat.html")'.format(prf_sign = prf_sign,roi_text = roi_text, hemi = hemi))
                    #     output_file(output_file_html, title='Subject: %s | ROI: %s | Hemisphere: %s | Sign: %s | Vertex: %i | Figures:pRF laterality histogram'%(subject,roi_text,hemi,prf_sign,vertex))
                    #     save(all_f3)

                    #     exec('output_file_html = opj(fig_bokeh_dir_{prf_sign}_{hemi},"{roi_text}_{hemi}_{prf_sign}_pRFtc.html")'.format(prf_sign = prf_sign,roi_text = roi_text, hemi = hemi))
                    #     output_file(output_file_html, title='Subject: %s | ROI: %s | Hemisphere: %s | Sign: %s | Figures: pRF time course'%(subject,roi_text,hemi,prf_sign))
                    #     save(all_f4)

                    # 
