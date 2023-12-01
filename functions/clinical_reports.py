"""
Generate the clinical reports for ZBRAINS

zbrains MUST exist on the global environment

USAGE:
    
    outdir='<PATH to output directory>/derivatives/z-brains'
    zbrains='<PATH to repository>/z-brains'

    # Run the function per norm ['z', 'normative'] to create one PDF per norm
    clinical_reports(norm='z', 
                     sub='sub-PX069', 
                     ses='ses-01', 
                     outdir=outdir, 
                     cmap='cmo.balance', 
                     cmap_asymmetry='PRGn',
                     Range=(-3,3),
                     thr=1.96, 
                     thr_alpha=0.3, 
                     color_bar='left', 
                     den='0p5', 
                     tmp='/tmp', 
                     sex='F', 
                     age='40',
                     analysis=['asymmetry', 'regional'], 
                     feature=['ADC', 'FA', 'T1map', 'flair', 'thickness', 'volume']
                     )

"""
# libraries
from xhtml2pdf import pisa
import os
import numpy as np
import glob
import nibabel as nib
from brainspace.datasets import load_mask
from brainspace.plotting import plot_hemispheres
from brainspace.mesh.mesh_io import read_surface
from brainspace.plotting.surface_plotting import plot_surf
from scipy.spatial.transform import Rotation as R
from brainspace.vtk_interface import wrap_vtk, serial_connect
from vtk import vtkPolyDataNormals
from pyvirtualdisplay import Display
import pandas as pd
import cmocean
cmaps = cmocean.cm.cmap_d

# -----------------------------------------------------------------------------
# Funtions for plotting
def plot_surfs(surfaces, values, filename=None, views =None, size=None, zoom=1.75, color_bar='bottom', share='both',
                Color_range=(-2, 2), cmap="cmo.balance", save=False, text=None, embed=False, interactive=True,
                transparent_bg=False, nan_color=(0, 0, 0, 1)):
    '''
    surfaces = [hipp_mid_l, hipp_unf_l, hipp_unf_r, hipp_mid_r]  Can be 1 or more
    views = ['dorsal', 'lateral', 'lateral', 'dorsal'] Can be 1 or more
    '''
    # Append values to surfaces
    my_surfs = {}
    array_names = []
    for i, value in enumerate(surfaces):
        surf_key = f's{i+1}'  # Assuming you want keys like 's1', 's2', ...
        surfaces[i].append_array(values[i], name=f'surf{i+1}')
        my_surfs[surf_key] = surfaces[i]
        array_names.append(f'surf{i+1}')
    
    # Set the layout as the list of keys from my_surfs
    layout = [list(my_surfs.keys())]

    # Size
    if size == None: size=(200*len(surfaces),350)
    
    p = plot_surf(my_surfs, layout, array_name=array_names, view=views, color_bar=color_bar,
                  color_range=Color_range, share=share, label_text=text, cmap=cmap,
                  nan_color=nan_color, zoom=zoom, background=(1, 1, 1),
                  size=size, embed_nb=embed, interactive=interactive, scale=(1,1),
                  transparent_bg=transparent_bg, screenshot=save, filename=filename,
                  return_plotter=False, suppress_warnings=False)

    return p

# -----------------------------------------------------------------------------
# Funtions for report
def convert_html_to_pdf(source_html, output_filename):
    # open output file for writing (truncated binary)
    result_file = open(output_filename, "w+b")

    # convert HTML to PDFt
    pisa_status = pisa.CreatePDF(
            source_html,                # the HTML to convert
            dest=result_file)           # file handle to recieve result

    # close output file
    result_file.close()                 # close output file

    # return True on success and False on errors
    return pisa_status.err

def report_header_template(zbrains='', analysis='', sub='', ses='', age='', sex=''):
    # Header
    report_header = (
        # Micapipe banner
        f'<img id=\"top\" src=\"{zbrains}/data/zbrains_banner.png\" alt=\"zbrains\">'

        # Title
        '<p style="margin-bottom:0;font-family:gill sans, sans-serif;text-align:center;font-size:20px;'
        f'color:#505050"> Clinical Summary &nbsp; | &nbsp; <b> {analysis} analysis </b> </p>'

        # Subject's ID - Session - Basic demographics
        '<p style="margin-bottom:0;margin-top:-100px;font-family:gill sans,sans-serif;text-align:center;font-size:14px;'
        f'color:#505050"> <b>Subject</b>: {sub}, <b>Session</b>: {ses}, <b>Sex</b>: {sex}, <b>Age</b>: {age}, '
    )

    return report_header

def report_colors(analysis='regional'): 
    report = ('<hr>')
    if analysis != 'regional':
        report += (
            '<p style="margin-bottom:0;margin-top:0;font-family:gill sans,sans-serif;text-align:center;font-size:14px;'
            'color:#7532a8"> <b> Purple </b> = <b>right</b> MORE THAN <b>left</b> </p>'
            '<p style="margin-bottom:0;margin-top:0;font-family:gill sans,sans-serif;text-align:center;font-size:14px;'
            'color:#32a852"> <b> Green </b> = <b>left</b> MORE THAN <b>right</b> </p>'
            )
    else:
        report += (
            '<p style="margin-bottom:0;margin-top:0;font-family:gill sans,sans-serif;text-align:center;font-size:14px;'
            'color:#b31b2c"> <b> Red </b> = INCREASED compared to controls </p>'
            '<p style="margin-bottom:0;margin-top:0;font-family:gill sans,sans-serif;text-align:center;font-size:14px;'
            'color:#13365d"> <b> Blue </b> = DECREASED compared to controls </p>'
            )
    return report


# Header template
def feature_header_template(feature='', extra=''):
    # Module header:
    feature_header_template = (
        '<p style="border:0px solid #666;padding-top:10px;padding-left:5px;background-color:#eee;font-family:Helvetica, '
        'sans-serif;font-size:14px;text-align:center;color:#5d5070">'
        f'<b> Feature: {feature} {extra}</b> </p>'
    )
    return feature_header_template

def report_1x2_table(fig1='', fig2='', height=250):
    # QC summary table
    report_fig_table = (
        '<table style="border:1px solid white;width:100%">'
            # Row 1: fig1, fig2
            f'<tr><td style=padding-top:4px;padding-bottom:4px;padding-left:3px;padding-right:3px;text-align:center><img style="height:{height}px;margin-top:-100px;" src="{fig1}"></td>'
            f'<td style=padding-top:4px;padding-bottom:4px;padding-left:3px;padding-right:3px;text-align:center><img style="height:{height}px;margin-top:-100px;" src="{fig2}"></td></tr>'
        '</table>'
    )
    return report_fig_table

def report_cortex(feature='', smooth_ctx=5, analysis='', norm='z', cmap='cmo.balance', Range=(-2,2), thr=None, 
                  color_bar='bottom', thr_alpha=0.5, subjdir=None, ses=None, sub=None,  cmap_asymmetry='PRGn', tmp='/tmp'):
    
    # Load native mid CORTICAL surface
    inf_lh = read_surface(f'{zbrains}/data/fsLR-32k.L.inflated.surf.gii', itype='gii')
    inf_rh = read_surface(f'{zbrains}/data/fsLR-32k.R.inflated.surf.gii', itype='gii')
    mask = load_mask(join=True)
    
    # Path to the feature data
    file_lh = f'{subjdir}/norm-{norm}/cortex/{sub}_{ses}_hemi-L_surf-fsLR-32k_feature-{feature}_smooth-{smooth_ctx}mm_analysis-{analysis}.func.gii'
    file_rh = f'{subjdir}/norm-{norm}/cortex/{sub}_{ses}_hemi-R_surf-fsLR-32k_feature-{feature}_smooth-{smooth_ctx}mm_analysis-{analysis}.func.gii'

    # Make figure
    def make_png(data_lh, data_rh, Range=(-2,2)):
        
        # threshold
        if thr != None:
            data_lh[np.abs(data_lh)<thr]*=thr_alpha
            data_rh[np.abs(data_rh)<thr]*=thr_alpha
            thr_str=f' | threshold={thr}'
        else:
            thr_str=''
        
        # outname of the plot
        out_png=f'{tmp}/{sub}_{ses}_cortex_feature-{feature}_smooth-{smooth_ctx}mm_analysis-{analysis}{thr}.png'

        # Replace values in f with NaN where mask_32k is False
        feat_map = np.hstack(np.concatenate((data_lh, data_rh), axis=0))
        feat_map[mask == False] = np.nan
        
        # Plot the mean FEATURE on fsLR-32k
        plot_hemispheres(inf_lh, inf_rh, array_name=feat_map, layout_style='grid',
                         label_text={'left': ['Lateral', 'Medial'], 'top': ['Left', 'Right']},
                         cmap=cmap, color_bar=color_bar, color_range=Range,
                         share=None, transparent_bg=False, nan_color=(0, 0, 0, 1),
                         zoom=1.25, size=(600, 600), embed_nb=True, offscreen=True,
                         interactive=False, screenshot=True, filename=out_png)
        
        # Plot cortex with different views
        plot_surfs(surfaces=[inf_lh, inf_rh, inf_rh, inf_lh],
                   values=[data_lh, data_rh, data_rh, data_lh],
                   views=['dorsal', 'dorsal', 'ventral', 'ventral'],
                   text={'bottom':['Left Superior','Right Superior','Right Inferior','Left Inferior']},
                   cmap=cmap, color_bar=color_bar, Color_range=Range,
                   share='both', transparent_bg=False, nan_color=(0, 0, 0, 1),
                   zoom=2.95, size=(550,350), embed=True,
                   interactive=False, save=True, filename=out_png.replace('.png', '_SI.png'))
        
        # create html chunck
        png_block = report_1x2_table(fig1=out_png, fig2=out_png.replace('.png', '_SI.png'), height=300)
        
        return(png_block, thr_str)
    
    # If the FEATURE is not found it will print a warning
    if not os.path.exists(file_lh):
        print(f"\n[WARNING] cortex file was not found: \n{file_lh}")
        png_block = '<p style="margin-bottom:0;margin-top:0;font-family:gill sans,sans-serif;text-align:center;font-size:14px;color:#ffb311"> <b> [WARNING] </b>cortex file was not found </p>'
        thr_str=''
    else:
        # Load the LEFT feature data
        data_lh = nib.load(file_lh).darrays[0].data
        # Load the RIGHT feature data
        if analysis != 'asymmetry':
            data_rh =nib.load(file_rh).darrays[0].data
        else:
            data_rh =nib.load(file_lh).darrays[0].data
            # LEFT is possitive RIGHT is negative
            data_lh[(data_lh)<0] = 0
            data_rh[(data_rh)>0] = 0
            cmap=cmap_asymmetry
            Range=(-1.5,1.5)
        # Create HTML png chunck
        png_block, thr_str = make_png(data_lh, data_rh, Range=Range)
        
    # Create cortical chunck
    ctx_block = (f'<p style="margin-bottom:0;margin-top:0;font-family:gill sans,sans-serif;'
          f'text-align:left;font-size:14px;color:#5d5070"> <b>Cortical {feature}</b> | {smooth_ctx}mm smooth | {analysis} analysis | norm-{norm} {thr_str}</p>')
    ctx_block += png_block
    
    return(ctx_block)

def subcorticalvertices(subcortical_values=None):
    """
    Taken from the ENIGMA toolbox https://github.com/MICA-MNI/ENIGMA/blob/master/enigmatoolbox/utils/parcellation.py#L353
    Map one value per subcortical area to surface vertices (author: @saratheriver)

    Parameters
    ----------
    subcortical_values : 1D ndarray
        Shape = (16,), order of subcortical structure must be = 
        L_accumbens, L_amygdala, L_caudate, L_hippocampus, L_pallidun, L_putamen, L_thalamus, L_ventricles, 
        R_accumbens, R_amygdala, R_caudate, R_hippocampus, R_pallidun, R_putamen, R_thalamus, R_ventricles

    zbrains output (the fisrt colum is SubjID):
        'Laccumb', 'Lamyg', 'Lcaud', 'Lhippo', 'Lpal', 'Lput','Lthal', 
        'Raccumb', 'Ramyg', 'Rcaud', 'Rhippo', 'Rpal', 'Rput', 'Rthal'
    
    Returns
    -------
    data : 1D ndarray
        Transformed data, shape = (51278,)
    """
    numvertices = [867, 1419, 3012, 3784, 1446, 4003, 3726, 7653, 838, 1457, 3208, 3742, 1373, 3871, 3699, 7180]
    data = []
    if isinstance(subcortical_values, np.ndarray):
        for ii in range(16):
            data.append(np.tile(subcortical_values[ii], (numvertices[ii], 1)))
        data = np.vstack(data).flatten()
    return data


def report_subcortex(feature='', analysis='', cmap='cmo.balance', norm='z', thr=None, Range=(-2,2), 
                     thr_alpha=0.5, subjdir=None, ses=None, sub=None, cmap_asymmetry='PRGn', tmp='/tmp'):

    # Read subcortical surfaces
    surf_lh = read_surface(f'{zbrains}/data/sctx.L.surf.gii', itype='gii')
    surf_rh = read_surface(f'{zbrains}/data/sctx.R.surf.gii', itype='gii')
    
    # List CSV file
    sctx_file = f'{subjdir}/norm-{norm}/subcortex/{sub}_{ses}_feature-{feature}_analysis-{analysis}.csv'
    
    # Make figure
    def make_png(array_values, cmap='cmo.balance', Range=(-2,2)):
        # threshold
        if thr != None:
            array_values[np.abs(array_values)<thr]*=thr_alpha
            thr_str=f' | threshold={thr}'
        else:
            thr_str=''
        
        # Array of data
        array_16 = np.full(16, np.nan)
        array_16[0:7] = array_values[0:7]
        if analysis != 'asymmetry':
            array_16[8:15] = array_values[7:]
        else:
            array_16[8:15] = array_values[0:7]
            # LEFT [0:7] is positive RiGHT [8:15] is negative
            array_16[0:7][(array_16[0:7])<0] = 0
            array_16[8:15][(array_16[8:15])>0] = 0
            cmap=cmap_asymmetry
            Range=(-1.5,1.5)
        
        # Map array to vertices of surfaces
        feat_map = subcorticalvertices(array_16)
            
        # Map array to left and right surfaces
        data_lh = feat_map[0:25910]
        data_rh = feat_map[25910:]
        
        # Plot subcortical structures
        out_png=f'{tmp}/{sub}_{ses}_subcortex_feature-{feature}_analysis-{analysis}{thr}.png'
        plot_surfs(surfaces=[surf_lh, surf_lh, surf_rh, surf_rh],
                   values=[data_lh, data_lh, data_rh, data_rh], 
                   views=['lateral', 'medial', 'lateral', 'medial'],
                   text={'left':['left'], 'right':['right']},
                   cmap=cmap, color_bar='bottom', Color_range=Range,
                   share='both', transparent_bg=True, nan_color=(0, 0, 0, 0),
                   zoom=1.4, size=(900, 250), embed=True,
                   interactive=False, save=True, filename=out_png)
        png_block = (f'<p style="text-align:center;margin-left=0px;"> <a href="{out_png}" target="_blank">'
        f'<img style="height:150px;margin-top:-100px;" src="{out_png}"> </a> </p>')
        return(png_block, thr_str)
    
    # If the FEATURE is not found it will print a warning
    if not os.path.exists(sctx_file):
        print(f"\n[WARNING] subcortex file was not found: \n{sctx_file}")
        png_block='<p style="margin-bottom:0;margin-top:0;font-family:gill sans,sans-serif;text-align:center;font-size:14px;color:#ffb311"> <b> [WARNING] </b>subcortex file was not found </p>'
        thr_str=''
    else:
        # Load CSV data
        sctx_data = pd.read_csv(sctx_file, sep=',')
        # Get the values into an array
        array_values = sctx_data.iloc[:, 1:].values.reshape(-1)
        # Create figure
        png_block, thr_str = make_png(array_values, Range=Range)
    
    # Create subcortical chunck
    sctx_block = (f'<p style="margin-bottom:0;margin-top:0;font-family:gill sans,sans-serif;'
          f'text-align:left;font-size:14px;color:#5d5070"> <b>Subcortical {feature}</b> | {analysis} analysis | norm-{norm} {thr_str}</p>')
    sctx_block += png_block

    return(sctx_block)


def report_hippocampus(feature='', analysis='', smooth_hipp=2, den='0p5', norm='z', thr=None, Range=(-2,2), cmap_asymmetry='PRGn',
                       cmap='cmo.balance', color_bar='bottom', thr_alpha=0.5, subjdir=None, ses=None, sub=None, tmp='/tmp'):
    
    # Load the cannonical HIPPOCAMPAL surfaces    
    hipp_mid_r=read_surface(f'{zbrains}/data/tpl-avg_space-canonical_den-{den}mm_label-hipp_midthickness.surf.gii', itype='gii')
    hipp_unf_r=read_surface(f'{zbrains}/data/tpl-avg_space-unfold_den-{den}mm_label-hipp_midthickness.surf.gii', itype='gii')
    hipp_mid_l=read_surface(f'{zbrains}/data/tpl-avg_space-canonical_den-{den}mm_label-hipp_midthickness.surf.gii', itype='gii')
    hipp_unf_l=read_surface(f'{zbrains}/data/tpl-avg_space-unfold_den-{den}mm_label-hipp_midthickness.surf.gii', itype='gii')

    # Flip Right to left surface
    vflip = np.ones(hipp_mid_l.Points.shape)
    vflip[:,0] = -1
    hipp_mid_l.Points = hipp_mid_l.Points*vflip
    hipp_unf_l.Points = hipp_unf_l.Points*vflip

    # Rotate surfaces because Reinder didn't acept my pull request

    # Rotate around Y-axis 270
    rY270 = R.from_rotvec(3*np.pi/2 * np.array([0, 1, 0]))
    hipp_unf_r.Points = rY270.apply(hipp_unf_r.Points)
    
    # Rotate around X-axis 90
    rY90 = R.from_rotvec(np.pi/2 * np.array([0, 1, 0]))
    hipp_unf_l.Points = rY90.apply(hipp_unf_l.Points)

    # Rotate around Z-axis 180 
    rZ = R.from_rotvec(np.pi * np.array([0, 0, 1]))
    hipp_unf_r.Points = rZ.apply(hipp_unf_r.Points)

    # Right Antero-posterior lateral
    hipp_lat_r = read_surface(f'{zbrains}/data/tpl-avg_space-canonical_den-{den}mm_label-hipp_midthickness.surf.gii', itype='gii')
    hipp_lat_r.Points = rY270.apply(hipp_lat_r.Points)

    # Left Antero-posterior lateral
    hipp_lat_l = read_surface(f'{zbrains}/data/tpl-avg_space-canonical_den-{den}mm_label-hipp_midthickness.surf.gii', itype='gii')
    hipp_lat_l.Points = rY90.apply(hipp_mid_l.Points)
    
    # Load the feature data
    file_lh=f'{subjdir}/norm-{norm}/hippocampus/{sub}_{ses}_hemi-L_den-{den}mm_feature-{feature}_smooth-{smooth_hipp}mm_analysis-{analysis}.func.gii'
    file_rh=f'{subjdir}/norm-{norm}/hippocampus/{sub}_{ses}_hemi-R_den-{den}mm_feature-{feature}_smooth-{smooth_hipp}mm_analysis-{analysis}.func.gii'
    
    # Make figure
    def make_png(feat_lh, feat_rh, Range=(-2,2)):
        
        # threshold
        if thr != None:
            feat_lh[np.abs(feat_lh)<thr]*=thr_alpha
            feat_rh[np.abs(feat_rh)<thr]*=thr_alpha
            thr_str=f' | threshold={thr}'
        else:
            thr_str=''
        out_png=f'{tmp}/{sub}_{ses}_hippocampus_feature-{feature}_smooth-{smooth_hipp}mm_analysis-{analysis}{thr}.png'
        
        # Plot the hippocampal data
        plot_surfs(surfaces=[hipp_lat_l, hipp_mid_l, hipp_unf_l, hipp_unf_r, hipp_mid_r, hipp_lat_r], 
                        values=[feat_lh, feat_lh, feat_lh, feat_rh, feat_rh, feat_rh], 
                        views=['dorsal', 'dorsal', 'lateral', 'lateral', 'dorsal', 'dorsal'], color_bar=color_bar,
                        zoom=1.75, cmap=cmap, Color_range=Range, interactive=False, save=True, filename=out_png)
        
        png_block = (f'<p style="text-align:center;margin-left=0px;"> <a href="{out_png}" target="_blank">'
                     f'<img style="height:175px;margin-top:-100px;" src="{out_png}"> </a> </p>' )
    
        return(png_block, thr_str)
    
    # If the FEATURE is not found it will print a warning
    if not os.path.exists(file_lh):
        print(f"\n[WARNING] hippocampal file was not found: \n{file_lh}")
        png_block = '<p style="margin-bottom:0;margin-top:0;font-family:gill sans,sans-serif;text-align:center;font-size:14px;color:#ffb311"> <b> [WARNING] </b>hippocampal file was not found </p>'
        thr_str=''
    else:
        # Load the LEFT feature data
        feat_lh = nib.load(file_lh).darrays[0].data
        # Load the RIGHT feature data
        if analysis != 'asymmetry':
            feat_rh = nib.load(file_rh).darrays[0].data
        else:
            feat_rh =nib.load(file_lh).darrays[0].data
            # left is positive RIGHT negative
            feat_lh[(feat_lh)<0] = 0
            feat_rh[(feat_rh)>0] = 0
            cmap=cmap_asymmetry
            Range=(-1.5,1.5)
        # create HTML chunck
        png_block, thr_str = make_png(feat_lh, feat_rh, Range=Range)
    
    # Create cortical chunck
    hipp_block = ('<p style="margin-bottom:0;margin-top:0;font-family:gill sans,sans-serif;'
          'text-align:left;font-size:14px;color:#5d5070"> <b>Hippocampal {feature}</b> | {smooth_hipp}mm smooth | density {den}mm | {analysis} analysis | norm-{norm} {thr_str}</p>'
          ).format(analysis=analysis, feature=feature, smooth_hipp=smooth_hipp, norm=norm, thr_str=thr_str, den=den.replace('p','.'))
    hipp_block += png_block
    
    return(hipp_block)


# -----------------------------------------------------------------------------
def clinical_reports(norm='z', sub=None, ses=None, outdir=None, cmap='cmo.balance', Range=(-2,2), 
                     thr=1.96, color_bar='bottom', den='0p5', smooth_ctx=5, smooth_hipp=2, sex='unknown', age='unknown',
                     thr_alpha=0.3, analysis=None, feature=None, cmap_asymmetry='PRGn', tmp='/tmp'):
    ''' Zbrains: Clinical report generator
    
        Global Parameters
        -----------------
        zbrains  :  str
        MUST exist and is the PATH to the z-brain repository <path_to/z-brains>
        
        Parameters
        ----------
        outdir  : str
            Output directory for processed z-brains derivatives.
        sub     : str
            BIDS_ID identification
        ses     : str
            OPTIONAL flag that indicates the session name (if omitted will manage as SINGLE session)
        norm    : str
            Default='z' MUST be one of the models used ['z', 'normative']
        features     : [] list or None
            Default=None it will run all the features unless a list of features is provided
        analysis     : ['regional', 'asymmetry'] list or None
            Default=None it will run all the analyses unless a list of analyses is provided
        cmap    : str
            Default='cmo.balance'. Colormap for plotting the regional changes.
        cmap_asymmetry
            Default='PRGn'. Colormap for plotting the regional changes.
        Range   : tuple
            Default=(-2,2). Range of values for the regional plot
        color_bar: str
            Default='right' position of the color bar for hippocampal and cortical plots
        den     : str
            Default='0p5' Density of the hippocampal maps ['2', '0p5']
        thr     : float
            Default=1.96. Threshold for the maps
        thr_alpha : float
            Default=0.3 Alpha channel for the colormap when thr is applied 0=white, 1=full color
        smooth_ctx     : int
            Default=5 Smoothing for the cortical data in mm.
        smooth_hipp     : int
            default=2 Smoothing for the hippocampal data in mm.
        age     : int
            Defaul='unknown'. OPTIONAL Covariate for the normative model. 
        sex     : str
            Defaul='unknown'. OPTIONAL MUST be [F,M] Covariate for the normative model.
        tmp     : str
            Default='/tmp'. Working directory for temporary files
    
      NOTE: clinical_reports does not take 'volume' as an feature.
            if it exist it will be removed. To process morphometry only 'thickness' is allowed
            
      McGill University, MNI, MICA-lab, Created 24 November 2023 'Zbrains Hackathon', MICA collective
      https://github.com/MICA-MNI/micapipe
      http://mica-mni.github.io/
    '''  

    # subject directory
    subjdir=f"{outdir}/{sub}/{ses}"
    # List all files of norm
    subses_files = glob.glob(f"{subjdir}/norm-{norm}/*/*")
    
    # If ANALYSIS is empty run all avaliable analysis
    if analysis == None:
        # get unique ANALYSIS
        analysis = sorted(list(set([file.split('.')[0].split('analysis-')[1].split('_')[0] for file in subses_files])))
    
    # If FEATURES is empty run all avaliable features
    if feature == None:
        # get unique FEATURES
        feature = sorted(list(set([file.split('feature-')[1].split('_')[0] for file in subses_files])))
    print(f'{sub}_{ses} features:\n{feature}')
    # Remove volume from the features
    if 'volume' in feature: feature.remove('volume')
    
    # Initialize the report
    static_report = ''
    
    # Display for headless plotting
    dsize = (900, 750)
    display = Display(visible=0, size=dsize)
    display.start()
    
    for A in analysis:
        for feat in feature:
            if feat == 'thickness':
                feat_sctx='volume'
            else:
                feat_sctx=feat
            print(f'Running summary of norm-{norm}, analysis-{A}, feature-{feat}')
            # -------------------------------------------------------------------------
            # PAGE 1: Analysis raw
            # Generate the report Title
            static_report += report_header_template(zbrains=zbrains, analysis=A, sub=sub, ses=ses, age=age, sex=sex)
            static_report += feature_header_template(feat)
            
            # Report cortex
            static_report +=  report_cortex(feature=feat, analysis=A, cmap=cmap, thr=None, norm=norm, Range=Range,
                                            color_bar=color_bar, smooth_ctx=smooth_ctx, 
                                            subjdir=subjdir, ses=ses, sub=sub, cmap_asymmetry=cmap_asymmetry, tmp=tmp)
            # Report subcortical
            static_report +=  report_subcortex(feature=feat_sctx, analysis=A, cmap=cmap, thr=None, norm=norm, Range=Range,
                                               subjdir=subjdir, ses=ses, sub=sub, cmap_asymmetry=cmap_asymmetry, tmp=tmp)
            # Report hippocampis
            static_report +=  report_hippocampus(feature=feat, analysis=A, cmap=cmap, thr=None, norm=norm, Range=Range,
                                                 color_bar=color_bar, den=den, smooth_hipp=smooth_hipp, 
                                                 subjdir=subjdir, ses=ses, sub=sub, cmap_asymmetry=cmap_asymmetry, tmp=tmp)
            static_report += report_colors(analysis=A)
            
            # page break
            static_report += '<div style="page-break-after: always;"></div>'
            
            # -------------------------------------------------------------------------
            # New page title
            static_report += report_header_template(zbrains=zbrains, analysis=A, sub=sub, ses=ses, age=age, sex=sex)
            static_report += feature_header_template(feat, extra=f'| Threshold: {thr}')
            
            # THRESHOLDED: Report cortex
            static_report +=  report_cortex(feature=feat, analysis=A, cmap=cmap, thr=thr, norm=norm, Range=Range,
                                            color_bar=color_bar, smooth_ctx=smooth_ctx, thr_alpha=thr_alpha, 
                                            subjdir=subjdir, ses=ses, sub=sub, cmap_asymmetry=cmap_asymmetry, tmp=tmp)
            # THRESHOLDED: Report subcortical 
            static_report +=  report_subcortex(feature=feat_sctx, analysis=A, cmap=cmap, thr=thr, norm=norm, Range=Range,
                                               thr_alpha=thr_alpha, subjdir=subjdir, ses=ses, sub=sub, cmap_asymmetry=cmap_asymmetry, tmp=tmp)
            # THRESHOLDED: Report hippocampis
            static_report +=  report_hippocampus(feature=feat, analysis=A, cmap=cmap, thr=thr, norm=norm, Range=Range, 
                                                 thr_alpha=thr_alpha, color_bar=color_bar, den=den, smooth_hipp=smooth_hipp, 
                                                 subjdir=subjdir, ses=ses, sub=sub, cmap_asymmetry=cmap_asymmetry, tmp=tmp)
            static_report += report_colors(analysis=A)
            
            # page break
            static_report += '<div style="page-break-after: always;"></div>'
    
    # Stop display for headless plotting
    display.stop()
    
    # Report file name
    file_pdf=f'{subjdir}/{sub}_{ses}_norm-{norm}_summary-report.pdf'
    print(f"Out file: {file_pdf}")
    
    # Create the HTML file
    convert_html_to_pdf(static_report, file_pdf)
