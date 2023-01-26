import os
import sys
import pandas as pd
import numpy as np
import nibabel as nib
import glob
from enigmatoolbox.utils.useful import zscore_matrix
from enigmatoolbox.datasets import load_mask
from enigmatoolbox.plotting import plot_cortical, plot_subcortical, plot_hippocampal
from enigmatoolbox.datasets import getaffine
import warnings

warnings.simplefilter('ignore')

import argparse
# defined command line options
# this also generates --help and error handling
CLI=argparse.ArgumentParser()
CLI.add_argument("subject", nargs=1, type=str, default="")
CLI.add_argument("session", nargs=1, type=str, default="")
CLI.add_argument("out", nargs=1, type=str, default="")
CLI.add_argument("hipDir", nargs=1, type=str, default="")
CLI.add_argument("demo", nargs=1, type=str, default="")
CLI.add_argument("thr", nargs=1, type=float, default=1.96)
CLI.add_argument(
  "--featList_ctx",  # name on the CLI - drop the `--` for positional/required parameters
  nargs="*",
  type=str,
  default=['flair', 'qt1', 'adc', 'thickness'])
CLI.add_argument(
  "--featList_sctx",
  nargs="*",
  type=str,  # any type/callable can be used here
  default=[])
CLI.add_argument(
  "--featList_hipp",
  nargs="*",
  type=str,  # any type/callable can be used here
  default=[])

# parse the command line
args = CLI.parse_args()

# access CLI options
subject = args.subject[0]
session = args.session[0]
out = args.out[0]
hipDir = args.hipDir[0]
demo = args.demo[0]
thr = args.thr[0]
featList_ctx = args.featList_ctx
featList_sctx = args.featList_sctx
featList_hipp = args.featList_hipp

def zscore(mtrx, TBL):
    grpk = TBL['grp']
    data_z = zscore_matrix(mtrx, grpk, 'HC')
    return data_z

def matrix(area, path, filename, TBL):
    # Get number of vertices/structures
    if area == "ctx":
        nvert = 64984
    elif area == "sctx":
        nvert = 14
    elif area == "hipp":
        nvert = 14524

    # Subject IDs
    sub = TBL['ID']
    (nSub,) = sub.shape

    # Initialize matrix
    mtrx = np.empty((nSub, nvert))
    mtrx[:] = np.nan

    # Cortical feature matrix
    if area == "ctx":
        for i in range(nSub):
            try:
                dpath = out + "/" + sub[i] + "/" + session + "/" + path
                d = []
                for _, h in enumerate(['lh', 'rh']):
                    d = np.append(d, nib.load(dpath + sub[i] + '_' + session +
                                            filename.format(h)).get_fdata().squeeze())
                mtrx[i, :] = d
            except:
                print(sub[i] + " : NO " + filename)
                pass
    
    # Subcortical feature matrix
    elif area == "sctx":
        for i in range(nSub):
            try:
                dpath = out + "/" + sub[i] + "/" + session + "/" + path
                mtrx[i] = np.loadtxt(dpath + sub[i] + "_" + session + filename, 
                                    delimiter=",", skiprows=1, usecols=range(1,15))
            except:
                print(sub[i] + " : NO " + filename)
                pass

    # Hippocampal feature matrix
    elif area == "hipp" and "_thickness" in filename:
        for i in range(nSub):
            try:
                dpath = hipDir + "/" + sub[i] + "/" + path
                d = []
                for _, h in enumerate(['L', 'R']):
                    d = np.append(d, nib.load(dpath + sub[i] +
                                            filename.format(h)).darrays[0].data)
                mtrx[i] = d
            except:
                print(sub[i] + " : NO " + filename)
                pass


    elif area == "hipp" and "_thickness" not in filename:
        for i in range(nSub):
            try:
                dpath = out + "/" + sub[i] + "/" + session + "/" + path
                d = []
                for _, h in enumerate(['lh', 'rh']):
                    d = np.append(d, nib.load(dpath + sub[i] + '_' + session +
                                            filename.format(h)).darrays[0].data)
                mtrx[i] = d
            except:
                print(sub[i] + " : NO " + filename)
                pass
    
    # Compute asymmetry
    mtrx = (mtrx[:, 0:int(nvert / 2)] - mtrx[:, int(nvert / 2):nvert]) / (
                (mtrx[:, 0:int(nvert / 2)] + mtrx[:, int(nvert / 2):nvert]) / 2)

    # Create dataframe
    mtrx = pd.DataFrame(mtrx)

    # Compute z-score data
    zmtrx = zscore(mtrx, TBL)
    return zmtrx

if __name__ == "__main__":
    # Load demographics
    TBL = pd.read_excel(demo, engine='openpyxl', skiprows=[0,2]).dropna(how="all").dropna(axis=1, how="all")
    #TBL = TBL.drop(TBL[TBL.incl == 0].index).reset_index()

    #=============================================================================
    # Load, concatenate, and save CORTICAL features (T2-FLAIR, qT1, AD, thickness)
    #=============================================================================
    mv_c = []
    mvdic_c = {}
    if "flair" in featList_ctx:
        t2z_c = matrix("ctx", '/anat/surfaces/flair/', '_space-conte69-32k_desc-{}_flair_10mm.mgh', TBL)
        mv_c.append(t2z_c.values)
        mvdic_c.update({'flair': t2z_c})
    if "qt1" in featList_ctx:
        qt1z_c = matrix("ctx", '/anat/surfaces/qt1/', '_space-conte69-32k_desc-{}_qt1_10mm.mgh', TBL)
        mv_c.append(qt1z_c.values)
        mvdic_c.update({'qt1': qt1z_c})
    if "adc" in featList_ctx:
        adcz_c = matrix("ctx", '/dwi/surfaces/', '_space-conte69-32k_desc-{}_model-DTI_map-ADC_10mm.mgh', TBL)
        mv_c.append(adcz_c.values)
        mvdic_c.update({'adc': adcz_c})
    if "thickness" in featList_ctx:
        ctz_c = matrix("ctx", '/anat/surfaces/morphology/', '_space-conte69-32k_desc-{}_thickness_10mm.mgh', TBL)
        ctz_c = ctz_c * -1
        mv_c.append(ctz_c.values)
        mvdic_c.update({'thickness': ctz_c})

    if mv_c:
        mv_c = np.nanmean(mv_c, axis=0)

        # Save unthrehsolded map, then filter out noise or non-significant findings and save
        m = load_mask(surface_name="conte69", join=True).astype(int)
        mv_c_orig = np.tile(mv_c, (1, 2)) * m
        mv_c_unthr = mv_c_orig
        np.savetxt(os.path.join(os.path.dirname(out), "analysis", "asymmetry", "allSubjects_unthr_ctx-z.csv"),
                   mv_c_orig, delimiter=",")
                   
        mv_c[(mv_c > -thr) & (mv_c < thr)] = 0
        mv_c = np.tile(mv_c, (1, 2)) * m
        np.savetxt(os.path.join(os.path.dirname(out), "analysis", "asymmetry", "allSubjects_ctx-z.csv"),
                   mv_c, delimiter=",")

        # Save individual features (unthresholded and thresholded)
        for (_, mvName) in enumerate(mvdic_c):
            mv = mvdic_c[mvName]
            np.savetxt(os.path.join(os.path.dirname(out), "analysis", "asymmetry", "allSubjects_unthr_ctx-{}.csv".
                                    format(str(mvName))), mv, delimiter=",")
            mv[(mv > -thr) & (mv < thr)] = 0
            mv = np.tile(mv, (1, 2)) * m
            np.savetxt(os.path.join(os.path.dirname(out), "analysis", "asymmetry", "allSubjects_ctx-{}.csv".
                                    format(str(mvName))), mv, delimiter=",")
    
    
    
    #================================================================================
    # Load, concatenate, and save SUBCORTICAL features (T2-FLAIR, qT1, AD, thickness)
    #================================================================================
    mv_s = []
    mvdic_s = {}
    if "flair" in featList_sctx:
        t2z_s = matrix("sctx", '/anat/surfaces/flair/', '_space-flair_subcortical-intensities.csv', TBL)
        mv_s.append(t2z_s.values)
        mvdic_s.update({'flair': t2z_s})
    if "qt1" in featList_sctx:
        qt1z_s = matrix("sctx", '/anat/surfaces/qt1/', '_space-qt1_subcortical-intensities.csv', TBL)
        mv_s.append(qt1z_s.values)
        mvdic_s.update({'qt1': qt1z_s})
    if "adc" in featList_sctx:
        adcz_s = matrix("sctx", '/dwi/surfaces/', '_space-dwi_subcortical-ADC.csv', TBL)
        mv_s.append(adcz_s.values)
        mvdic_s.update({'adc': adcz_s})
    if "thickness" in featList_sctx:
        ctz_s = matrix("sctx", '/anat/surfaces/morphology/sctx_volume/', '_sctx_volume.csv', TBL)
        ctz_s = ctz_s * -1
        mv_s.append(ctz_s.values)
        mvdic_s.update({'thickness': ctz_s})

    if mv_s:
        mv_s = np.nanmean(mv_s, axis=0)
        mv_s_unthr = np.tile(mv_s, (1, 2))

        # Save unthrehsolded map, then filter out noise or non-significant findings and save
        np.savetxt(os.path.join(os.path.dirname(out), "analysis", "asymmetry", "allSubjects_unthr_sctx-z.csv"),
                   mv_s, delimiter=",")
        mv_s[(mv_s > -thr) & (mv_s < thr)] = 0
        mv_s = np.tile(mv_s, (1, 2))
        np.savetxt(os.path.join(os.path.dirname(out), "analysis", "asymmetry", "allSubjects_sctx-z.csv"),
                   mv_s, delimiter=",")

        # Save individual features (unthresholded and thresholded)
        for (_, mvName) in enumerate(mvdic_s):
            mv = mvdic_s[mvName]
            np.savetxt(os.path.join(os.path.dirname(out), "analysis", "asymmetry", "allSubjects_unthr_sctx-{}.csv".
                                    format(str(mvName))), mv, delimiter=",")
            mv[(mv > -thr) & (mv < thr)] = 0
            mv = np.tile(mv, (1, 2))
            np.savetxt(os.path.join(os.path.dirname(out), "analysis", "asymmetry", "allSubjects_sctx-{}.csv".
                                    format(str(mvName))), mv, delimiter=",")
    
    
    
    #=========================================================================================
    # Load, concatenate, and save HIPPOCAMPAL features (T2-FLAIR, qT1, AD, thickness)
    #=========================================================================================
    mv_h = []
    mvdic_h = {}
    if "flair" in featList_hipp:
        t2z_h = matrix("hipp", '/anat/surfaces/flair/', '_hemi-{}_space-flair_desc-flair_N4_den-0p5mm_label-hipp_midthickness_10mm.func.gii', TBL)
        mv_h.append(t2z_h.values)
        mvdic_h.update({'flair': t2z_h})
    if "qt1" in featList_hipp:
        qt1z_h = matrix("hipp", '/anat/surfaces/qt1/', '_hemi-{}_space-qt1_desc-qt1_den-0p5mm_label-hipp_midthickness_10mm.func.gii', TBL)
        mv_h.append(qt1z_h.values)
        mvdic_h.update({'qt1': qt1z_h})
    if "adc" in featList_hipp:
        adcz_h = matrix("hipp", '/dwi/surfaces/', '_hemi-{}_space-dwi_desc-dwi-ADC_den-0p5mm_label-hipp_midthickness_10mm.func.gii', TBL)
        mv_h.append(adcz_h.values)
        mvdic_h.update({'adc': adcz_h})
    if "thickness" in featList_hipp:
        ctz_h = matrix("hipp", '/surf/', '_hemi-{}_space-T1w_den-0p5mm_label-hipp_thickness.shape.gii', TBL)
        ctz_h = ctz_h * -1
        mv_h.append(ctz_h.values)
        mvdic_h.update({'thickness': ctz_h})

    if mv_h:
        mv_h = np.nanmean(mv_h, axis=0)
        mv_h_unthr = np.tile(mv_h, (1, 2))

        # Save unthrehsolded map, then filter out noise or non-significant findings and save
        np.savetxt(os.path.join(os.path.dirname(out), "analysis", "asymmetry", "allSubjects_unthr_hipp-z.csv"),
                   mv_h, delimiter=",")
        mv_h[(mv_h > -thr) & (mv_h < thr)] = 0
        mv_h = np.tile(mv_h, (1, 2))
        np.savetxt(os.path.join(os.path.dirname(out), "analysis", "asymmetry", "allSubjects_hipp-z.csv"),
                   mv_h, delimiter=",")

        # Save individual features (unthresholded and thresholded)
        for (_, mvName) in enumerate(mvdic_h):
            mv = mvdic_h[mvName]
            np.savetxt(os.path.join(os.path.dirname(out), "analysis", "asymmetry", "allSubjects_unthr_hipp-{}.csv".
                                    format(str(mvName))), mv, delimiter=",")
            mv[(mv > -thr) & (mv < thr)] = 0
            mv = np.tile(mv, (1, 2))
            np.savetxt(os.path.join(os.path.dirname(out), "analysis", "asymmetry", "allSubjects_hipp-{}.csv".
                                    format(str(mvName))), mv, delimiter=",")
    
    
    
    # ==================================================================
    # Plot and save multivariate cortical, subcortical, and hippocampal asymmetry
    # ==================================================================
    # Get row number of subject
    rn = (TBL.loc[TBL['ID'] == subject].index).tolist()

    # Save and plot multivariate z-score | cortical
    mv_c_half = mv_c_unthr[rn, :]
    fname = os.path.join(os.path.dirname(out), "analysis", "asymmetry", subject, session, subject + "_ctx-mz-unthr_{}.mgh")
    nib.freesurfer.mghformat.MGHImage(np.float32(mv_c_half[:, :mv_c_half.shape[1]//2].flatten()),
                                      getaffine('conte69', 'lh')).to_filename(fname.format('lh'))
    nib.freesurfer.mghformat.MGHImage(np.float32(mv_c_half[:, mv_c_half.shape[1]//2:].flatten()),
                                      getaffine('conte69', 'rh')).to_filename(fname.format('rh'))
    
    mv_c_half = mv_c[rn, :]
    fname = os.path.join(os.path.dirname(out), "analysis", "asymmetry", subject, session, subject + "_ctx-mz_{}.mgh")
    nib.freesurfer.mghformat.MGHImage(np.float32(mv_c_half[:, :mv_c_half.shape[1]//2].flatten()),
                                      getaffine('conte69', 'lh')).to_filename(fname.format('lh'))
    nib.freesurfer.mghformat.MGHImage(np.float32(mv_c_half[:, mv_c_half.shape[1]//2:].flatten()),
                                      getaffine('conte69', 'rh')).to_filename(fname.format('rh'))

    mv_c_half[:, :len(mv_c_half.flatten())//2] = np.nan
    fname = os.path.join(os.path.dirname(out), "analysis", "asymmetry", subject, session, subject + "_ctx-mz.png")
    plot_cortical(array_name=mv_c_half, surface_name='conte69', size=(800, 180), zoom=1.18, cmap='RdBu_r',
                  color_bar=True, color_range=(-3, 3), screenshot=True, view=['lateral', 'medial', 'medial', 'lateral'],
                  filename=fname)

    # Save and plot multivariate z-score | subcortical
    mv_s_half = mv_s_unthr[rn, :].flatten()
    np.savetxt(os.path.join(os.path.dirname(out), "analysis", "asymmetry", subject, session, subject + "_sctx-mz-unthr.txt"), 
               mv_s_half, delimiter=",")
    
    mv_s_half = mv_s[rn, :].flatten()
    np.savetxt(os.path.join(os.path.dirname(out), "analysis", "asymmetry", subject, session, subject + "_sctx-mz.txt"), 
               mv_s_half, delimiter=",")
    
    mv_s_half[:len(mv_s_half)//2] = np.nan
    fname = os.path.join(os.path.dirname(out), "analysis", "asymmetry", subject, session, subject + "_sctx-mz.png")
    plot_subcortical(array_name=mv_s_half, ventricles=False, size=(800, 180), zoom=1.18, cmap='RdBu_r',
                  color_bar=True, color_range=(-3, 3), screenshot=True, view=['lateral', 'medial', 'medial', 'lateral'],
                  filename=fname)

    # Save and plot multivariate z-score | hippocampal
    mv_h_half = mv_h_unthr[rn, :].flatten()
    mv_h_lh = nib.gifti.gifti.GiftiImage()
    mv_h_lh.add_gifti_data_array(nib.gifti.gifti.GiftiDataArray(data=mv_h_half[:len(mv_h_half)//2]))
    mv_h_rh = nib.gifti.gifti.GiftiImage()
    mv_h_rh.add_gifti_data_array(nib.gifti.gifti.GiftiDataArray(data=mv_h_half[len(mv_h_half)//2:]))
    fname = os.path.join(os.path.dirname(out), "analysis", "asymmetry", subject, session, subject + "_hipp-mz-unthr_{}.gii")
    nib.save(mv_h_lh, fname.format('lh'))
    nib.save(mv_h_rh, fname.format('rh'))

    mv_h_half = mv_h[rn, :].flatten()
    mv_h_lh = nib.gifti.gifti.GiftiImage()
    mv_h_lh.add_gifti_data_array(nib.gifti.gifti.GiftiDataArray(data=mv_h_half[:len(mv_h_half)//2]))
    mv_h_rh = nib.gifti.gifti.GiftiImage()
    mv_h_rh.add_gifti_data_array(nib.gifti.gifti.GiftiDataArray(data=mv_h_half[len(mv_h_half)//2:]))
    fname = os.path.join(os.path.dirname(out), "analysis", "asymmetry", subject, session, subject + "_hipp-mz_{}.gii")
    nib.save(mv_h_lh, fname.format('lh'))
    nib.save(mv_h_rh, fname.format('rh'))

    mv_h_half[:len(mv_h_half) // 2] = np.nan
    fname = os.path.join(os.path.dirname(out), "analysis", "asymmetry", subject, session, subject + "_hipp-mz.png")
    plot_hippocampal(array_name=mv_h_half, size=(800, 180), zoom=1.18, cmap='RdBu_r',
                     color_bar=True, color_range=(-3, 3), screenshot=True,
                     view=['lateral', 'medial', 'ventral', 'dorsal'], filename=fname)
    
    
    
    # ================================================================
    # Plot and save univariate cortical, subcortical, and hippocampal asymmetry
    # ================================================================
    # Plot univariate z-score | cortical
    for (_, mvName) in enumerate(mvdic_c):
        mv_tmp = np.loadtxt(os.path.join(os.path.dirname(out), "analysis", "asymmetry", "allSubjects_ctx-{}.csv".
                                    format(str(mvName))), delimiter=",")[rn, :]
        mv_tmp[:, :len(mv_c_half.flatten())//2] = np.nan
        plot_cortical(array_name=mv_tmp, surface_name='conte69', size=(800, 180), zoom=1.18, cmap='RdBu_r',
                  color_bar=True, color_range=(-3, 3), screenshot=True, view=['lateral', 'medial', 'medial', 'lateral'],
                  filename=os.path.join(os.path.dirname(out), "analysis", "asymmetry", subject, session, subject + "_ctx-{}.png".
                  format(str(mvName))))
        
        mv_tmp = np.loadtxt(os.path.join(os.path.dirname(out), "analysis", "asymmetry", "allSubjects_unthr_ctx-{}.csv".
                                    format(str(mvName))), delimiter=",")[rn, :]
        fname = os.path.join(os.path.dirname(out), "analysis", "asymmetry", subject, session, subject + "_ctx-{}-unthr_{}.mgh")
        nib.freesurfer.mghformat.MGHImage(np.float32(mv_tmp.flatten()),
                                          getaffine('conte69', 'lh')).to_filename(fname.format(str(mvName), 'lh'))
        nib.freesurfer.mghformat.MGHImage(np.float32(mv_tmp.flatten()),
                                          getaffine('conte69', 'rh')).to_filename(fname.format(str(mvName), 'rh'))


    # Plot univariate z-score | subcortical
    for (_, mvName) in enumerate(mvdic_s):
        mv_tmp = np.loadtxt(os.path.join(os.path.dirname(out), "analysis", "asymmetry", "allSubjects_sctx-{}.csv".
                                    format(str(mvName))), delimiter=",")[rn, :].flatten()
        mv_tmp[:len(mv_tmp)//2] = np.nan
        plot_subcortical(array_name=mv_tmp, ventricles=False, size=(800, 180), zoom=1.18, cmap='RdBu_r',
                  color_bar=True, color_range=(-3, 3), screenshot=True, view=['lateral', 'medial', 'medial', 'lateral'],
                  filename=os.path.join(os.path.dirname(out), "analysis", "asymmetry", subject,
                                        session, subject + "_sctx-{}.png".format(str(mvName))))
        
        mv_tmp = np.loadtxt(os.path.join(os.path.dirname(out), "analysis", "asymmetry", "allSubjects_unthr_sctx-{}.csv").
                                    format(str(mvName)), delimiter=",")[rn, :].flatten()
        np.savetxt(os.path.join(os.path.dirname(out), "analysis", "asymmetry", subject, session, subject + "_sctx-{}-unthr.txt").
                                    format(str(mvName)), mv_tmp, delimiter=",")
        

    # Plot univariate z-score | hippocampal
    for (_, mvName) in enumerate(mvdic_h):
        mv_tmp = np.loadtxt(os.path.join(os.path.dirname(out), "analysis", "asymmetry", "allSubjects_hipp-{}.csv".
                                    format(str(mvName))), delimiter=",")[rn, :].flatten()
        mv_tmp[:len(mv_tmp) // 2] = np.nan
        plot_hippocampal(array_name=mv_tmp, size=(800, 180), zoom=1.18, cmap='RdBu_r',
                     color_bar=True, color_range=(-3, 3), screenshot=True,
                     view=['lateral', 'medial', 'ventral', 'dorsal'],
                     filename=os.path.join(os.path.dirname(out), "analysis", "asymmetry", subject,
                                           session, subject + "_hipp-{}.png".format(str(mvName))))
        
        mv_tmp = np.loadtxt(os.path.join(os.path.dirname(out), "analysis", "asymmetry", "allSubjects_unthr_hipp-{}.csv".
                    format(str(mvName))), delimiter=",")[rn, :].flatten()
        mv_h_lh = nib.gifti.gifti.GiftiImage()
        mv_h_lh.add_gifti_data_array(nib.gifti.gifti.GiftiDataArray(data=mv_tmp))
        mv_h_rh = nib.gifti.gifti.GiftiImage()
        mv_h_rh.add_gifti_data_array(nib.gifti.gifti.GiftiDataArray(data=mv_tmp))
        fname = os.path.join(os.path.dirname(out), "analysis", "asymmetry", subject, session, subject + "_hipp-{}-unthr_{}.gii")
        nib.save(mv_h_lh, fname.format(str(mvName), 'lh'))
        nib.save(mv_h_rh, fname.format(str(mvName), 'rh'))
    
    
