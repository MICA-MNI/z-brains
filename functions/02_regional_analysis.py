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
  default=['flair', 'qt1', 'adc', 'thickness'])
CLI.add_argument(
  "--featList_hipp",
  nargs="*",
  type=str,  # any type/callable can be used here
  default=['flair', 'qt1', 'adc', 'thickness'])

# parse the command line
args = CLI.parse_args()

# access CLI options
subject = args.subject[0]
session = args.session[0]
out = args.out[0]
demo = args.demo[0]
thr = args.thr[0]
featList_ctx = args.featList_ctx
featList_sctx = args.featList_sctx
featList_hipp = args.featList_hipp

def zscore(mtrx, TBL):
    grpk = TBL['grp']
    data_z = zscore_matrix(mtrx, grpk, 'HC')
    return data_z

def matrix(area, filename, TBL):
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
                dpath = f"{out}/../z-brains/scene-nativepro/{sub[i]}/{session}/{sub[i]}_{session}"
                d = []
                for _, h in enumerate(['lh', 'rh']):
                    d = np.append(d, nib.load(dpath + filename.format(h)).darrays[0].data)
                mtrx[i, :] = d
            except:
                print(sub[i] + " : NO " + dpath + filename.format(h))
                pass
    
    # Subcortical feature matrix
    elif area == "sctx":
        for i in range(nSub):
            try:
                dpath = f"{out}/../z-brains/scene-nativepro/{sub[i]}/{session}/{sub[i]}_{session}"
                mtrx[i] = np.loadtxt(dpath + filename, 
                                    delimiter=",", skiprows=1, usecols=range(1,15))
            except:
                print(sub[i] + " : NO " + dpath + filename)
                pass

    # Hippocampal feature matrix
    elif area == "hipp":
        for i in range(nSub):
            try:
                dpath = f"{out}/../z-brains/scene-nativepro/{sub[i]}/{session}/{sub[i]}_{session}"
                d = []
                for _, h in enumerate(['lh', 'rh']):
                    d = np.append(d, nib.load(dpath + filename.format(h)).darrays[0].data)
                mtrx[i] = d
            except:
                print(sub[i] + " : NO " + dpath + filename.format(h))
                pass
    
    # Get vertexwise z-scores
    mtrx_full = mtrx
    mtrx_full = pd.DataFrame(mtrx_full)
    zmtrx_full = zscore(mtrx_full, TBL)
    
    return zmtrx_full

if __name__ == "__main__":
    # Load demographics
    TBL = pd.read_excel(demo, engine='openpyxl', skiprows=[0,2]).dropna(how="all").dropna(axis=1, how="all")
    #TBL = TBL.drop(TBL[TBL.incl == 0].index).reset_index()

    #=============================================================================
    # Load, concatenate, and save CORTICAL features (T2-FLAIR, qT1, AD, thickness)
    #=============================================================================
    mvFull_c = []
    mvFullDic_c = {}
    if "flair" in featList_ctx:
        t2zFull_c = matrix("ctx", "_space-conte69_hemi-{}_midthickness_desc-flair_10mm.func.gii", TBL)
        mvFull_c.append(t2zFull_c.values)
        mvFullDic_c.update({'flair': t2zFull_c})
    if "qt1" in featList_ctx:
        qt1zFull_c = matrix("ctx", "_space-conte69_hemi-{}_midthickness_desc-qt1_10mm.func.gii", TBL)
        mvFull_c.append(qt1zFull_c.values)
        mvFullDic_c.update({'qt1': qt1zFull_c})
    if "adc" in featList_ctx:
        adczFull_c = matrix("ctx", "_space-conte69_hemi-{}_midthickness_desc-ADC_10mm.func.gii", TBL)
        mvFull_c.append(adczFull_c.values)
        mvFullDic_c.update({'adc': adczFull_c})
    if "thickness" in featList_ctx:
        ctzFull_c = matrix("ctx", "_space-conte69_hemi-{}_desc-thickness_10mm.func.gii", TBL)
        ctzFull_c = ctzFull_c * -1
        mvFull_c.append(ctzFull_c.values)
        mvFullDic_c.update({'thickness': ctzFull_c})

    if mvFull_c:
        mvFull_c = np.nanmean(mvFull_c, axis=0)

        # Save unthrehsolded map, then filter out noise or non-significant findings and save
        m = load_mask(surface_name="conte69", join=True).astype(int)
        mv_c_orig = mvFull_c * m
        mvFull_c_unthr = mv_c_orig
        np.savetxt(os.path.join(os.path.dirname(out), "z-brains", "regional", "allSubjects_unthr_ctx-z.csv"),
                   mv_c_orig, delimiter=",")
        mvFull_c[(mvFull_c > -thr) & (mvFull_c < thr)] = 0
        np.savetxt(os.path.join(os.path.dirname(out), "z-brains", "regional", "allSubjects_ctx-z.csv"),
                   mvFull_c, delimiter=",")

        # Save individual features (unthresholded and thresholded)
        for (_, mvName) in enumerate(mvFullDic_c):
            mv = mvFullDic_c[mvName]
            np.savetxt(os.path.join(os.path.dirname(out), "z-brains", "regional", "allSubjects_unthr_ctx-{}.csv".
                                    format(str(mvName))), mv, delimiter=",")
            mv[(mv > -thr) & (mv < thr)] = 0
            np.savetxt(os.path.join(os.path.dirname(out), "z-brains", "regional", "allSubjects_ctx-{}.csv".
                                    format(str(mvName))), mv, delimiter=",")




    #================================================================================
    # Load, concatenate, and save SUBCORTICAL features (T2-FLAIR, qT1, AD, thickness)
    #================================================================================
    mvFull_s = []
    mvFullDic_s = {}
    if "flair" in featList_sctx:
        t2zFull_s = matrix("sctx", '_space-nativepro_desc-subcortical-flair.csv', TBL)
        mvFull_s.append(t2zFull_s.values)
        mvFullDic_s.update({'flair': t2zFull_s})
    if "qt1" in featList_sctx:
        qt1zFull_s = matrix("sctx", '_space-nativepro_desc-subcortical-qt1.csv', TBL)
        mvFull_s.append(qt1zFull_s.values)
        mvFullDic_s.update({'qt1': qt1zFull_s})
    if "adc" in featList_sctx:
        adczFull_s = matrix("sctx", '_space-nativepro_desc-subcortical-ADC.csv', TBL)
        mvFull_s.append(adczFull_s.values)
        mvFullDic_s.update({'adc': adczFull_s})
    if "thickness" in featList_sctx:
        ctzFull_s = matrix("sctx", '_space-nativepro_desc-subcortical-volume.csv', TBL)
        ctzFull_s = ctzFull_s * -1
        mvFull_s.append(ctzFull_s.values)
        mvFullDic_s.update({'thickness': ctzFull_s})

    if mvFull_s:
        mvFull_s = np.nanmean(mvFull_s, axis=0)
        mvFull_s_unthr = mvFull_s

        # Save unthrehsolded map, then filter out noise or non-significant findings and save
        np.savetxt(os.path.join(os.path.dirname(out), "z-brains", "regional", "allSubjects_unthr_sctx-z.csv"),
                   mvFull_s_unthr, delimiter=",")
        mvFull_s[(mvFull_s > -thr) & (mvFull_s < thr)] = 0
        np.savetxt(os.path.join(os.path.dirname(out), "z-brains", "regional", "allSubjects_sctx-z.csv"),
                   mvFull_s, delimiter=",")

        # Save individual features (unthresholded and thresholded)
        for (_, mvName) in enumerate(mvFullDic_s):
            mv = mvFullDic_s[mvName]
            np.savetxt(os.path.join(os.path.dirname(out), "z-brains", "regional", "allSubjects_unthr_sctx-{}.csv".
                                    format(str(mvName))), mv, delimiter=",")
            mv[(mv > -thr) & (mv < thr)] = 0
            np.savetxt(os.path.join(os.path.dirname(out), "z-brains", "regional", "allSubjects_sctx-{}.csv".
                                    format(str(mvName))), mv, delimiter=",")
    
    
    
    #=========================================================================================
    # Load, concatenate, and save HIPPOCAMPAL features (T2-FLAIR, qT1, AD, thickness)
    #=========================================================================================
    mvFull_h = []
    mvFullDic_h = {}
    if "flair" in featList_hipp:
        t2zFull_h = matrix("hipp", '_space-hipp_hemi-{}_midthickness_desc-flair_2mm.func.gii', TBL)
        mvFull_h.append(t2zFull_h.values)
        mvFullDic_h.update({'flair': t2zFull_h})
    if "qt1" in featList_hipp:
        qt1zFull_h = matrix("hipp", '_space-hipp_hemi-{}_midthickness_desc-qt1_2mm.func.gii', TBL)
        mvFull_h.append(qt1zFull_h.values)
        mvFullDic_h.update({'qt1': qt1zFull_h})
    if "adc" in featList_hipp:
        adczFull_h = matrix("hipp", '_space-hipp_hemi-{}_midthickness_desc-ADC_2mm.func.gii', TBL)
        mvFull_h.append(adczFull_h.values)
        mvFullDic_h.update({'adc': adczFull_h})
    if "thickness" in featList_hipp:
        ctzFull_h = matrix("hipp", '_space-hipp_hemi-{}_desc-thickness_2mm.func.gii', TBL)
        ctzFull_h = ctzFull_h * -1
        mvFull_h.append(ctzFull_h.values)
        mvFullDic_h.update({'thickness': ctzFull_h})

    if mvFull_h:
        mvFull_h = np.nanmean(mvFull_h, axis=0)
        mvFull_h_unthr = mvFull_h

        # Save unthrehsolded map, then filter out noise or non-significant findings and save
        np.savetxt(os.path.join(os.path.dirname(out), "z-brains", "regional", "allSubjects_unthr_hipp-z.csv"),
                   mvFull_h_unthr, delimiter=",")
        mvFull_h[(mvFull_h > -thr) & (mvFull_h < thr)] = 0
        np.savetxt(os.path.join(os.path.dirname(out), "z-brains", "regional", "allSubjects_hipp-z.csv"),
                   mvFull_h, delimiter=",")

        # Save individual features (unthresholded and thresholded)
        for (_, mvName) in enumerate(mvFullDic_h):
            mv = mvFullDic_h[mvName]
            np.savetxt(os.path.join(os.path.dirname(out), "z-brains", "regional", "allSubjects_unthr_hipp-{}.csv".
                                    format(str(mvName))), mv, delimiter=",")
            mv[(mv > -thr) & (mv < thr)] = 0
            np.savetxt(os.path.join(os.path.dirname(out), "z-brains", "regional", "allSubjects_hipp-{}.csv".
                                    format(str(mvName))), mv, delimiter=",")
    
    
    
    # ========================================================================
    # Plot multivariate cortical, subcortical, and hippocampal features (full)
    # ========================================================================
    # Get row number of subject
    rn = (TBL.loc[TBL['ID'] == subject].index).tolist()

    # Save and plot multivariate z-score | cortical
    mv_c_full = mvFull_c_unthr[rn, :]
    fname = os.path.join(os.path.dirname(out), "z-brains", "regional", subject, session, subject + "_ctx-mz-unthr_{}.mgh")
    nib.freesurfer.mghformat.MGHImage(np.float32(mv_c_full[:, :mv_c_full.shape[1]//2].flatten()),
                                      getaffine('conte69', 'lh')).to_filename(fname.format('lh'))
    nib.freesurfer.mghformat.MGHImage(np.float32(mv_c_full[:, mv_c_full.shape[1]//2:].flatten()),
                                      getaffine('conte69', 'rh')).to_filename(fname.format('rh'))
    
    mv_c_full = mvFull_c[rn, :]
    fname = os.path.join(os.path.dirname(out), "z-brains", "regional", subject, session, subject + "_ctx-mz_{}.mgh")
    nib.freesurfer.mghformat.MGHImage(np.float32(mv_c_full[:, :mv_c_full.shape[1]//2].flatten()),
                                      getaffine('conte69', 'lh')).to_filename(fname.format('lh'))
    nib.freesurfer.mghformat.MGHImage(np.float32(mv_c_full[:, mv_c_full.shape[1]//2:].flatten()),
                                      getaffine('conte69', 'rh')).to_filename(fname.format('rh'))

    fname = os.path.join(os.path.dirname(out), "z-brains", "regional", subject, session, subject + "_ctx-mz.png")
    plot_cortical(array_name=mv_c_full, surface_name='conte69', size=(800, 180), zoom=1.18, cmap='RdBu_r',
                  color_bar=True, color_range=(-3, 3), screenshot=True, view=['lateral', 'medial', 'medial', 'lateral'],
                  filename=fname)

    # Save and plot multivariate z-score | subcortical
    mv_tmp = np.loadtxt(os.path.join(os.path.dirname(out), "z-brains", "regional", "allSubjects_unthr_sctx-z.csv"),
                            delimiter=",")[rn, :].flatten()
    np.savetxt(os.path.join(os.path.dirname(out), "z-brains", "regional", subject, session, subject + "_sctx-mz-unthr.txt"), 
               mv_tmp, delimiter=",")
    
    mv_s_full = mvFull_s[rn, :].flatten()
    np.savetxt(os.path.join(os.path.dirname(out), "z-brains", "regional", subject, session, subject + "_sctx-mz.txt"), 
               mv_s_full, delimiter=",")
    
    fname = os.path.join(os.path.dirname(out), "z-brains", "regional", subject, session, subject + "_sctx-mz.png")
    plot_subcortical(array_name=mv_s_full, ventricles=False, size=(800, 180), zoom=1.18, cmap='RdBu_r',
                  color_bar=True, color_range=(-3, 3), screenshot=True, view=['lateral', 'medial', 'medial', 'lateral'],
                  filename=fname)

    # Save and plot multivariate z-score | hippocampal
    mv_h_full = np.loadtxt(os.path.join(os.path.dirname(out), "z-brains", "regional", "allSubjects_unthr_hipp-z.csv"),
                            delimiter=",")[rn, :].flatten()
    mv_h_lh = nib.gifti.gifti.GiftiImage()
    mv_h_lh.add_gifti_data_array(nib.gifti.gifti.GiftiDataArray(data=mv_h_full[:len(mv_h_full)//2]))
    mv_h_rh = nib.gifti.gifti.GiftiImage()
    mv_h_rh.add_gifti_data_array(nib.gifti.gifti.GiftiDataArray(data=mv_h_full[len(mv_h_full)//2:]))
    fname = os.path.join(os.path.dirname(out), "z-brains", "regional", subject, session, subject + "_hipp-mz-unthr_{}.gii")
    nib.save(mv_h_lh, fname.format('lh'))
    nib.save(mv_h_rh, fname.format('rh'))

    mv_h_full = mvFull_h[rn, :].flatten()
    mv_h_lh = nib.gifti.gifti.GiftiImage()
    mv_h_lh.add_gifti_data_array(nib.gifti.gifti.GiftiDataArray(data=mv_h_full[:len(mv_h_full)//2]))
    mv_h_rh = nib.gifti.gifti.GiftiImage()
    mv_h_rh.add_gifti_data_array(nib.gifti.gifti.GiftiDataArray(data=mv_h_full[len(mv_h_full)//2:]))
    fname = os.path.join(os.path.dirname(out), "z-brains", "regional", subject, session, subject + "_hipp-mz_{}.gii")
    nib.save(mv_h_lh, fname.format('lh'))
    nib.save(mv_h_rh, fname.format('rh'))

    fname = os.path.join(os.path.dirname(out), "z-brains", "regional", subject, session, subject + "_hipp-mz.png")
    plot_hippocampal(array_name=mv_h_full, size=(800, 180), zoom=1.18, cmap='RdBu_r',
                     color_bar=True, color_range=(-3, 3), screenshot=True,
                     view=['lateral', 'ventral', 'dorsal', 'medial'], filename=fname)
    
    
    
    # ======================================================================
    # Plot univariate cortical, subcortical, and hippocampal features (full)
    # ======================================================================
    # Plot univariate z-score | cortical
    for (_, mvName) in enumerate(mvFullDic_c):
        mv_tmp = np.loadtxt(os.path.join(os.path.dirname(out), "z-brains", "regional", "allSubjects_ctx-{}.csv".
                                    format(str(mvName))), delimiter=",")[rn, :]
        plot_cortical(array_name=mv_tmp, surface_name='conte69', size=(800, 180), zoom=1.18, cmap='RdBu_r',
                  color_bar=True, color_range=(-3, 3), screenshot=True, view=['lateral', 'medial', 'medial', 'lateral'],
                  filename=os.path.join(os.path.dirname(out), "z-brains", "regional", subject, session, subject + "_ctx-{}.png".
                  format(str(mvName))))
        
        mv_tmp = np.loadtxt(os.path.join(os.path.dirname(out), "z-brains", "regional", "allSubjects_unthr_ctx-{}.csv".
                                    format(str(mvName))), delimiter=",")[rn, :]
        fname = os.path.join(os.path.dirname(out), "z-brains", "regional", subject, session, subject + "_ctx-{}-unthr_{}.mgh")
        nib.freesurfer.mghformat.MGHImage(np.float32(mv_tmp[:, :mv_tmp.shape[1]//2].flatten()),
                                          getaffine('conte69', 'lh')).to_filename(fname.format(str(mvName), 'lh'))
        nib.freesurfer.mghformat.MGHImage(np.float32(mv_tmp[:, mv_tmp.shape[1]//2:].flatten()),
                                          getaffine('conte69', 'rh')).to_filename(fname.format(str(mvName), 'rh'))


    # Plot univariate z-score | subcortical
    for (_, mvName) in enumerate(mvFullDic_s):
        mv_tmp = np.loadtxt(os.path.join(os.path.dirname(out), "z-brains", "regional", "allSubjects_sctx-{}.csv".
                                    format(str(mvName))), delimiter=",")[rn, :].flatten()
        plot_subcortical(array_name=mv_tmp, ventricles=False, size=(800, 180), zoom=1.18, cmap='RdBu_r',
                  color_bar=True, color_range=(-3, 3), screenshot=True, view=['lateral', 'medial', 'medial', 'lateral'],
                  filename=os.path.join(os.path.dirname(out), "z-brains", "regional", subject,
                                        session, subject + "_sctx-{}.png".format(str(mvName))))
        
        mv_tmp = np.loadtxt(os.path.join(os.path.dirname(out), "z-brains", "regional", "allSubjects_unthr_sctx-{}.csv").
                                    format(str(mvName)), delimiter=",")[rn, :].flatten()
        np.savetxt(os.path.join(os.path.dirname(out), "z-brains", "regional", subject, session, subject + "_sctx-{}-unthr.txt").
                                    format(str(mvName)), mv_tmp, delimiter=",")
        

    # Plot univariate z-score | hippocampal
    for (_, mvName) in enumerate(mvFullDic_h):
        mv_tmp = np.loadtxt(os.path.join(os.path.dirname(out), "z-brains", "regional", "allSubjects_hipp-{}.csv".
                                    format(str(mvName))), delimiter=",")[rn, :].flatten()
        plot_hippocampal(array_name=mv_tmp, size=(800, 180), zoom=1.18, cmap='RdBu_r',
                     color_bar=True, color_range=(-3, 3), screenshot=True,
                     view=['lateral', 'ventral', 'dorsal', 'medial'],
                     filename=os.path.join(os.path.dirname(out), "z-brains", "regional", subject,
                                           session, subject + "_hipp-{}.png".format(str(mvName))))
        
        mv_tmp = np.loadtxt(os.path.join(os.path.dirname(out), "z-brains", "regional", "allSubjects_unthr_hipp-{}.csv".
                    format(str(mvName))), delimiter=",")[rn, :].flatten()
        mv_h_lh = nib.gifti.gifti.GiftiImage()
        mv_h_lh.add_gifti_data_array(nib.gifti.gifti.GiftiDataArray(data=mv_tmp[:len(mv_tmp)//2]))
        mv_h_rh = nib.gifti.gifti.GiftiImage()
        mv_h_rh.add_gifti_data_array(nib.gifti.gifti.GiftiDataArray(data=mv_tmp[len(mv_tmp)//2:]))
        fname = os.path.join(os.path.dirname(out), "z-brains", "regional", subject, session, subject + "_hipp-{}-unthr_{}.gii")
        nib.save(mv_h_lh, fname.format(str(mvName), 'lh'))
        nib.save(mv_h_rh, fname.format(str(mvName), 'rh'))
