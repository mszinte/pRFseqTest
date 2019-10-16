import os
import nibabel as nb
import numpy as np
import matplotlib.pyplot as pl
import scipy as sp
from scipy.signal import savgol_filter
from skimage.transform import rotate
from math import *
import cortex


def set_pycortex_config_file(data_folder):

    # Import necessary modules
    import os
    import cortex
    import ipdb
    from pathlib import Path

    # Define the new database and colormaps folder
    pycortex_db_folder = data_folder + '/pp_data/cortex/db/'
    pycortex_cm_folder = data_folder + '/pp_data/cortex/colormaps/'

    # Get pycortex config file location
    pycortex_config_file  = cortex.options.usercfg

    # Create name of new config file that will be written
    new_pycortex_config_file = pycortex_config_file[:-4] + '_new.cfg'

    # Create the new config file
    Path(new_pycortex_config_file).touch()

    # Open the config file in read mode and the newly created one in write mode.
    # Loop over every line in the original file and copy it into the new one.
    # For the lines containing either 'filestore' or 'colormap', it will
    # change the saved folder path to the newly created one above (e.g. pycortex_db_folder)
    with open(pycortex_config_file, 'r') as fileIn:
        with open(new_pycortex_config_file, 'w') as fileOut:

            for line in fileIn:

                if 'filestore' in line:
                    newline = 'filestore=' + pycortex_db_folder
                    fileOut.write(newline)
                    newline = '\n'

                elif 'Colormaps' in line:
                    newline = 'Colormaps=' + pycortex_cm_folder
                    fileOut.write(newline)
                    newline = '\n'

                else:
                    newline = line

                fileOut.write(newline)

    
    # Renames the original config file als '_old' and the newly created one to the original name
    os.rename(pycortex_config_file, pycortex_config_file[:-4] + '_old.cfg')
    os.rename(new_pycortex_config_file, pycortex_config_file)
    return None



def convert_fit_results(est_fn,
                        output_dir,
                        stim_width,
                        stim_height,
                        fit_model):
    """
    Convert pRF fitting value in different parameters for following analysis
   
    Parameters
    ----------
    est_fn: absolute paths to estimates file
    output_dir: absolute path to directory into which to put the resulting files
    stim_width: stimulus width in deg
    stim_heigth: stimulus height in deg
    fit_model: fit model ('gauss','css')

    Returns
    -------
    prf_deriv_all: derivative of pRF analysis for all pRF voxels
    prf_deriv_neg: derivative of pRF analysis for all negative pRF voxels
    prf_deriv_pos: derivative of pRF analysis for all positive pRF voxels

    stucture output:
    columns: 1->size of input
    row00 : sign
    row01 : R2
    row02 : eccentricity in deg
    row03 : polar angle real component in deg
    row04 : polar angle imaginary component in deg
    row05 : size in deg
    row06 : non-linerity or nans
    row07 : amplitude
    row08 : baseline
    row09 : coverage
    row10 : x
    row11 : y
    ['prf_sign','prf_rsq','prf_ecc','prf_polar_real','prf_polar_imag','prf_size','prf_non_lin','prf_amp','prf_baseline','prf_cov','prf_x','prf_y']
    """

    # Imports
    # -------
    # General imports
    import os
    import nibabel as nb
    import glob
    import numpy as np
    import ipdb
    deb = ipdb.set_trace

    
    # Popeye imports
    from popeye.spinach import generate_og_receptive_fields

    # Create folders
    # --------------
    try:
        os.makedirs(os.path.join(output_dir,'all'))
        os.makedirs(os.path.join(output_dir,'pos'))
        os.makedirs(os.path.join(output_dir,'neg'))
    except:
        pass


    # Get data details
    # ----------------
    est = []
    img_est = nb.load(est_fn)
    est = img_est.get_fdata()

    # Compute derived measures from prfs
    # ----------------------------------
    # get data index
    if fit_model == 'gauss':
        x_idx, y_idx, sigma_idx, beta_idx, baseline_idx, rsq_idx = 0, 1, 2, 3, 4, 5
    elif fit_model == 'css':
        x_idx, y_idx, sigma_idx, non_lin_idx, beta_idx, baseline_idx, rsq_idx = 0, 1, 2, 3, 4, 5, 6


    # change to nan empty voxels
    est[est[:,:,:,rsq_idx] == 0] = np.nan

    # pRF sign
    prf_sign_all = np.sign((est[:,:,:,beta_idx]))
    pos_mask = prf_sign_all > 0.0
    neg_mask = prf_sign_all < 0.0
    all_mask = pos_mask | neg_mask
    
    # r-square
    prf_rsq_all = est[:,:,:,rsq_idx]

    # pRF eccentricity
    prf_ecc_all = np.nan_to_num(np.sqrt(est[:,:,:,x_idx]**2 + est[:,:,:,y_idx]**2))

    # pRF polar angle
    complex_polar = est[:,:,:,x_idx] + 1j * est[:,:,:,y_idx]
    normed_polar = complex_polar / np.abs(complex_polar)
    prf_polar_real_all = np.real(normed_polar)
    prf_polar_imag_all = np.imag(normed_polar)
    
    # pRF size
    prf_size_all = est[:,:,:,sigma_idx].astype(np.float64)
    prf_size_all[prf_size_all<1e-4] = 1e-4

    # pRF non-linearity
    if fit_model == 'gauss':
        prf_non_lin_all = np.zeros((prf_size_all.shape))*np.nan
    elif fit_model == 'css':
        prf_non_lin_all = est[:,:,:,non_lin_idx]

    # pRF amplitude
    prf_amp_all = est[:,:,:,beta_idx]
    
    # pRF baseline
    prf_baseline_all = est[:,:,:,baseline_idx]

    # pRF coverage
    deg_x, deg_y = np.meshgrid(np.linspace(-30, 30, 121), np.linspace(-30, 30, 121))         # define prfs in visual space
    
    flat_est = est.reshape((-1, est.shape[-1])).astype(np.float64)

    rfs = generate_og_receptive_fields( flat_est[:,x_idx],
                                        flat_est[:,y_idx],
                                        flat_est[:,sigma_idx],
                                        flat_est[:,beta_idx].T*0+1,
                                        deg_x,
                                        deg_y)
    if fit_model == 'css':
        rfs = rfs ** flat_est[:,non_lin_idx]
    

    total_prf_content = rfs.reshape((-1, flat_est.shape[0])).sum(axis=0)
    log_x = np.logical_and(deg_x >= -stim_width/2.0, deg_x <= stim_width/2.0)
    log_y = np.logical_and(deg_y >= -stim_height/2.0, deg_y <= stim_height/2.0)
    stim_vignet = np.logical_and(log_x,log_y)
    prf_cov_all = rfs[stim_vignet, :].sum(axis=0) / total_prf_content
    prf_cov_all = prf_cov_all.reshape(prf_baseline_all.shape)
    
    # pRF x
    prf_x_all = est[:,:,:,x_idx]

    # pRF y
    prf_y_all = est[:,:,:,y_idx]

    # Save results
    prf_deriv_all = np.zeros((est.shape[0],est.shape[1],est.shape[2],12))*np.nan
    prf_deriv_pos = np.zeros((est.shape[0],est.shape[1],est.shape[2],12))*np.nan
    prf_deriv_neg = np.zeros((est.shape[0],est.shape[1],est.shape[2],12))*np.nan
    output_list = ['prf_sign','prf_rsq','prf_ecc','prf_polar_real','prf_polar_imag','prf_size','prf_non_lin','prf_amp','prf_baseline','prf_cov','prf_x','prf_y']

    for mask_dir in ['all','pos','neg']:
        print('saving: %s'%('os.path.join(output_dir,"{mask_dir}","prf_deriv_{mask_dir}.nii.gz")'.format(mask_dir = mask_dir)))
        for output_num, output_type in enumerate(output_list):
            exec('{output_type}_{mask_dir} = np.copy({output_type}_all)'.format(mask_dir = mask_dir, output_type = output_type))
            exec('{output_type}_{mask_dir}[~{mask_dir}_mask] = np.nan'.format(mask_dir = mask_dir, output_type = output_type))

            exec('prf_deriv_{mask_dir}[...,{output_num}]  = {output_type}_{mask_dir}'.format(mask_dir = mask_dir, output_type = output_type, output_num = output_num))

        exec('prf_deriv_{mask_dir} = prf_deriv_{mask_dir}.astype(np.float32)'.format(mask_dir = mask_dir))
        exec('new_img = nb.Nifti1Image(dataobj = prf_deriv_{mask_dir}, affine = img_est.affine, header = img_est.header)'.format(mask_dir = mask_dir))
        exec('new_img.to_filename(os.path.join(output_dir,"{mask_dir}","prf_deriv_{mask_dir}.nii.gz"))'.format(mask_dir = mask_dir))

    return None
    
def draw_cortex_vertex(subject,xfmname,data,cmap,vmin,vmax,cbar = 'discrete',cmap_steps = 255,\
                        alpha = None,depth = 1,thick = 1,height = 1024,sampler = 'nearest',\
                        with_curvature = True,with_labels = False,with_colorbar = False,\
                        with_borders = False,curv_brightness = 0.95,curv_contrast = 0.05,add_roi = False,\
                        roi_name = 'empty',col_offset = 0, zoom_roi = None, zoom_hem = None, zoom_margin = 0.0):
    """
    Plot brain data onto a previously saved flatmap.
    Parameters
    ----------
    subject             : subject id (e.g. 'sub-001')
    xfmname             : xfm transform
    data                : the data you would like to plot on a flatmap
    cmap                : colormap that shoudl be used for plotting
    vmin                : minimal value iqn colormap
    vmax                : maximal value in colormap
    cbar                : color bar layout
    cmap_steps          : number of colormap bins
    alpha               : alpha map
    depth               : Value between 0 and 1 for how deep to sample the surface for the flatmap (0 = gray/white matter boundary, 1 = pial surface)
    thick               : Number of layers through the cortical sheet to sample. Only applies for pixelwise = True
    height              : Height of the image to render. Automatically scales the width for the aspect of the subject's flatmap
    sampler             : Name of sampling function used to sample underlying volume data. Options include 'trilinear', 'nearest', 'lanczos'
    with_curvature      : Display the rois, labels, colorbar, annotated flatmap borders, or cross-hatch dropout?
    with_labels         : Display labels?
    with_colorbar       : Display pycortex' colorbar?
    with_borders        : Display borders?
    curv_brightness     : Mean brightness of background. 0 = black, 1 = white, intermediate values are corresponding grayscale values.
    curv_contrast       : Contrast of curvature. 1 = maximal contrast (black/white), 0 = no contrast (solid color for curvature equal to curvature_brightness).
    add_roi             : add roi -image- to overlay.svg
    roi_name            : roi name
    col_offset          : colormap offset between 0 and 1
    zoom_roi            : name of the roi on which to zoom on
    zoom_hem            : hemifield fo the roi zoom
    zoom_margin         : margin in mm around the zoom
    Returns
    -------
    vertex_rgb - pycortex vertex file
    """
    
    import cortex
    import numpy as np
    import matplotlib.pyplot as plt
    import matplotlib.colors as colors
    from matplotlib import cm
    import matplotlib as mpl
    import ipdb
    deb = ipdb.set_trace

    # define colormap
    base = cortex.utils.get_cmap(cmap)
    val = np.fmod(np.linspace(0+col_offset, 1+col_offset,cmap_steps+1,endpoint=False),1.0)
    colmap = colors.LinearSegmentedColormap.from_list(  'my_colmap',
                                                        base(val),
                                                        N = cmap_steps)

    # convert data to RGB
    vrange = float(vmax) - float(vmin)
    norm_data = ((data-float(vmin))/vrange)*cmap_steps
    mat = colmap(norm_data.astype(int))

    deb()

    vertex_rgb = cortex.VolumeRGB(  red = mat[...,0],
                                    green = mat[...,1],
                                    blue = mat[...,2],
                                    alpha = alpha,
                                    subject = subject,
                                    xfmname = xfmname
                                    )
    deb()
    
    vertex_rgb_fig = cortex.quickshow(  braindata = vertex_rgb,
                                        depth = depth,
                                        thick = thick,
                                        height = height,
                                        sampler = sampler,
                                        with_curvature = with_curvature,
                                        with_labels = with_labels,
                                        with_colorbar = with_colorbar,
                                        with_borders = with_borders,
                                        curvature_brightness = curv_brightness,
                                        curvature_contrast = curv_contrast)
    
    deb()
    zoom_plot = []
    if zoom_roi != None:
        roi_verts = cortex.get_roi_verts(subject, zoom_roi)[zoom_roi]
        roi_map = cortex.Vertex.empty(subject)
        roi_map.data[roi_verts] = 1

        (lflatpts, lpolys), (rflatpts, rpolys) = cortex.db.get_surf(subject, "flat", nudge=True)
        sel_pts = dict(left=lflatpts, right=rflatpts)[zoom_hem]
        roi_pts = sel_pts[np.nonzero(getattr(roi_map, zoom_hem))[0],:2]
        
        xmin, ymin = roi_pts.min(0) - zoom_margin
        xmax, ymax = roi_pts.max(0) + zoom_margin
        zoom_width = roi_pts.max(0)[0]- roi_pts.min(0)[0]
        zoom_height = roi_pts.max(0)[1]- roi_pts.min(0)[1]
        sqr_zoom_size = np.round(np.max([zoom_width,zoom_height]))*1.1
        zoom_ctr_x = np.mean([roi_pts.max(0)[0],roi_pts.min(0)[0]])
        zoom_ctr_y = np.mean([roi_pts.max(0)[1],roi_pts.min(0)[1]])
        mat = [zoom_ctr_x-sqr_zoom_size/2.0,zoom_ctr_x+sqr_zoom_size/2.0,zoom_ctr_y-sqr_zoom_size/2.0,zoom_ctr_y+sqr_zoom_size/2.0]
        zoom_plot = pl.axis(mat)
    
    else:
        
        if cbar == 'polar':

            # Polar angle color bar
            colorbar_location = [0.5, 0.07, 0.8, 0.2]
            n = 200
            cbar_axis = vertex_rgb_fig.add_axes(colorbar_location, projection='polar')
            norm = mpl.colors.Normalize(0, 2*np.pi)

            # Plot a color mesh on the polar plot
            # with the color set by the angle
            t = np.linspace(2*np.pi,0,n)
            r = np.linspace(1,0,2)
            rg, tg = np.meshgrid(r,t)
            c = tg
            im = cbar_axis.pcolormesh(t, r, c.T,norm= norm, cmap = colmap)
            cbar_axis.set_theta_zero_location("W")
            cbar_axis.set_yticklabels([])
            cbar_axis.set_xticklabels([])
            cbar_axis.spines['polar'].set_visible(False)

        elif cbar == 'ecc':
        
            # Ecc color bar
            colorbar_location = [0.5, 0.07, 0.8, 0.2]
            n = 200
            cbar_axis = vertex_rgb_fig.add_axes(colorbar_location, projection='polar')

            t = np.linspace(0,2*np.pi, n)
            r = np.linspace(0,1, n)
            rg, tg = np.meshgrid(r,t)
            c = tg
            
            im = cbar_axis.pcolormesh(t, r, c, norm = mpl.colors.Normalize(0, 2*np.pi), cmap = colmap)
            cbar_axis.tick_params(pad = 1,labelsize = 15)
            cbar_axis.spines['polar'].set_visible(False)
            
            # superimpose new axis for dva labeling
            box = cbar_axis.get_position()
            cbar_axis.set_yticklabels([])
            cbar_axis.set_xticklabels([])
            axl = vertex_rgb_fig.add_axes(  [1.8*box.xmin,
                                        0.5*(box.ymin+box.ymax),
                                        box.width/600,
                                        box.height*0.5],
                                        axisbg = None)
            axl.spines['top'].set_visible(False)
            axl.spines['right'].set_visible(False)
            axl.spines['bottom'].set_visible(False)
            axl.yaxis.set_ticks_position('right')
            axl.xaxis.set_ticks_position('none')
            axl.set_xticklabels([])
            axl.set_yticklabels(np.linspace(vmin,vmax,3),size = 'x-large')
            axl.set_ylabel('$dva$\t\t', rotation = 0, size = 'x-large')
            axl.yaxis.set_label_coords(box.xmax+30,0.4)
            axl.patch.set_alpha(0.5)

        elif cbar == 'discrete':

            # Discrete color bars
            # -------------------
            colorbar_location= [0.9, 0.05, 0.03, 0.25]
            cmaplist = [colmap(i) for i in range(colmap.N)]

            # define the bins and normalize
            bounds = np.linspace(vmin, vmax, cmap_steps + 1)
            bounds_label = np.linspace(vmin, vmax, 3)
            norm = mpl.colors.BoundaryNorm(bounds, colmap.N)
            
            cbar_axis = vertex_rgb_fig.add_axes(colorbar_location)
            cb = mpl.colorbar.ColorbarBase(cbar_axis,cmap = colmap,norm = norm,ticks = bounds_label,boundaries = bounds)

        # add to overalt
        if add_roi == True:
            cortex.utils.add_roi(   data = vertex_rgb,
                                    name = roi_name,
                                    open_inkscape = False,
                                    add_path = False,
                                    depth = depth,
                                    thick = thick,
                                    sampler = sampler,
                                    with_curvature = with_curvature,
                                    with_colorbar = with_colorbar,
                                    with_borders = with_borders,
                                    curvature_brightness = curv_brightness,
                                    curvature_contrast = curv_contrast)

    return vertex_rgb
