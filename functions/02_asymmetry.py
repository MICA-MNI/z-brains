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

def new_gifti_image(data, intent=0, datatype=16, metadata=None):
  """NiBabel wrapper to generate a gifti image with data array and metadata.
  Parameters
  ----------
  data : ndarray
    1-D ndarray containing one hemisphere surface data.
  intent : int
    Intent code for Gifti File. Defaults to 0 (Intent = NONE).
    Available intent codes:
      NIFTI_INTENT_NONE - 0
      NIFTI_INTENT_CORREL - 2
      NIFTI_INTENT_TTEST - 3
      NIFTI_INTENT_ZSCORE - 5
      NIFTI_INTENT_PVAL - 22
      NIFTI_INTENT_LOGPVAL - 23
      NIFTI_INTENT_LOG10PVAL - 24
      NIFTI_INTENT_LABEL - 1002
      NIFTI_INTENT_POINTSET - 1008
      NIFTI_INTENT_TRIANGLE - 1009
      NIFTI_INTENT_TIME_SERIES - 2001
      NIFTI_INTENT_NODE_INDEX - 2002
      NIFTI_INTENT_SHAPE - 2005
    More intent codes can be found at: https://nifti.nimh.nih.gov/nifti-1/documentation/nifti1fields/nifti1fields_pages/group__NIFTI1__INTENT__CODES.html
  datatype : int
    Datatype for gifti image. Defaults to 16 (dtype = float32)
    Available datatypes:
      UINT8 - 2
      INT32 - 8
      FLOAT32 - 16
  metadata : nibabel.gifti.gifti.GiftiMetaData
    Metadata for gifti image.
  Returns
  -------
  nibabel.gifti.gifti.GiftiImage
    Gifti image with specified metadata and data array.
  """
  dtypes = {2: np.uint8, 8: np.int32, 16: np.float32}
  data = data.astype(dtypes[datatype])
  if metadata:
    metadata = nib.gifti.GiftiMetaData(metadata)
  gifti_data = nib.gifti.GiftiDataArray(data=data, intent=intent, datatype=datatype)
  gifti_img = nib.gifti.GiftiImage(meta=metadata, darrays=[gifti_data])
  return gifti_img


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
                dpath = f"{out}/../z-brains/scene-nativepro/{sub[i]}_{session}/"
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
                dpath = f"{out}/../z-brains/scene-nativepro/{sub[i]}_{session}/"
                mtrx[i] = np.loadtxt(dpath + sub[i] + "_" + session + filename, 
                                    delimiter=",", skiprows=1, usecols=range(1,15))
            except:
                print(sub[i] + " : NO " + filename)
                pass

    # Hippocampal feature matrix
    elif area == "hipp":
        for i in range(nSub):
            try:
                dpath = f"{out}/../z-brains/scene-nativepro/{sub[i]}_{session}/"
                d = []
                for _, h in enumerate(['L', 'R']):
                    d = np.append(d, nib.load(dpath + sub[i] +
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
    mvDic_c = {}
    if "flair" in featList_ctx:
        t2z_c = matrix("ctx", "_space-conte69_hemi-{}_midthickness_desc-flair_10mm.func.gii", TBL)
        mv_c.append(t2z_c.values)
        mvDic_c.update({'flair': t2z_c})
    if "qt1" in featList_ctx:
        qt1z_c = matrix("ctx", "_space-conte69_hemi-{}_midthickness_desc-qt1_10mm.func.gii", TBL)
        mv_c.append(qt1z_c.values)
        mvDic_c.update({'qt1': qt1z_c})
    if "adc" in featList_ctx:
        adcz_c = matrix("ctx", "_space-conte69_hemi-{}_midthickness_desc-ADC_10mm.func.gii", TBL)
        mv_c.append(adcz_c.values)
        mvDic_c.update({'adc': adcz_c})
    if "thickness" in featList_ctx:
        ctz_c = matrix("ctx", "_space-conte69_hemi-{}_desc-thickness_10mm.func.gii", TBL)
        ctz_c = ctz_c * -1
        mv_c.append(ctz_c.values)
        mvDic_c.update({'thickness': ctz_c})

    if mv_c:
        mv_c = np.nanmean(mv_c, axis=0)

        # Save unthrehsolded map, then filter out noise or non-significant findings and save
        m = load_mask(surface_name="conte69", join=True).astype(int)
        mv_c_orig = np.tile(mv_c, (1, 2)) * m
        mv_c_unthr = mv_c_orig
        np.savetxt(os.path.join(os.path.dirname(out), "z-brains", "asymmetry", "allSubjects_unthr_ctx-z.csv"),
                   mv_c_orig, delimiter=",")
                   
        mv_c[(mv_c > -thr) & (mv_c < thr)] = 0
        mv_c = np.tile(mv_c, (1, 2)) * m
        np.savetxt(os.path.join(os.path.dirname(out), "z-brains", "asymmetry", "allSubjects_ctx-z.csv"),
                   mv_c, delimiter=",")

        # Save individual features (unthresholded and thresholded)
        for (_, mvName) in enumerate(mvDic_c):
            mv = mvDic_c[mvName]
            np.savetxt(os.path.join(os.path.dirname(out), "z-brains", "asymmetry", "allSubjects_unthr_ctx-{}.csv".
                                    format(str(mvName))), mv, delimiter=",")
            mv[(mv > -thr) & (mv < thr)] = 0
            mv = np.tile(mv, (1, 2)) * m
            np.savetxt(os.path.join(os.path.dirname(out), "z-brains", "asymmetry", "allSubjects_ctx-{}.csv".
                                    format(str(mvName))), mv, delimiter=",")
    
    
    
    #================================================================================
    # Load, concatenate, and save SUBCORTICAL features (T2-FLAIR, qT1, AD, thickness)
    #================================================================================
    mv_s = []
    mvDic_s = {}
    if "flair" in featList_sctx:
        t2z_s = matrix("sctx", '_subcortical-flair.csv', TBL)
        mv_s.append(t2z_s.values)
        mvDic_s.update({'flair': t2z_s})
    if "qt1" in featList_sctx:
        qt1z_s = matrix("sctx", '_subcortical-qt1.csv', TBL)
        mv_s.append(qt1z_s.values)
        mvDic_s.update({'qt1': qt1z_s})
    if "adc" in featList_sctx:
        adcz_s = matrix("sctx", '_subcortical-ADC.csv', TBL)
        mv_s.append(adcz_s.values)
        mvDic_s.update({'adc': adcz_s})
    if "thickness" in featList_sctx:
        ctz_s = matrix("sctx", '_subcortical-thickness.csv', TBL)
        ctz_s = ctz_s * -1
        mv_s.append(ctz_s.values)
        mvDic_s.update({'thickness': ctz_s})

    if mv_s:
        mv_s = np.nanmean(mv_s, axis=0)
        mv_s_unthr = np.tile(mv_s, (1, 2))

        # Save unthrehsolded map, then filter out noise or non-significant findings and save
        np.savetxt(os.path.join(os.path.dirname(out), "z-brains", "asymmetry", "allSubjects_unthr_sctx-z.csv"),
                   mv_s, delimiter=",")
        mv_s[(mv_s > -thr) & (mv_s < thr)] = 0
        mv_s = np.tile(mv_s, (1, 2))
        np.savetxt(os.path.join(os.path.dirname(out), "z-brains", "asymmetry", "allSubjects_sctx-z.csv"),
                   mv_s, delimiter=",")

        # Save individual features (unthresholded and thresholded)
        for (_, mvName) in enumerate(mvDic_s):
            mv = mvDic_s[mvName]
            np.savetxt(os.path.join(os.path.dirname(out), "z-brains", "asymmetry", "allSubjects_unthr_sctx-{}.csv".
                                    format(str(mvName))), mv, delimiter=",")
            mv[(mv > -thr) & (mv < thr)] = 0
            mv = np.tile(mv, (1, 2))
            np.savetxt(os.path.join(os.path.dirname(out), "z-brains", "asymmetry", "allSubjects_sctx-{}.csv".
                                    format(str(mvName))), mv, delimiter=",")
    
    
    
    #=========================================================================================
    # Load, concatenate, and save HIPPOCAMPAL features (T2-FLAIR, qT1, AD, thickness)
    #=========================================================================================
    mv_h = []
    mvDic_h = {}
    if "flair" in featList_hipp:
        t2z_h = matrix("hipp", '_space-hipp_hemi-{}_midthickness_desc-flair_2mm.func.gii', TBL)
        mv_h.append(t2z_h.values)
        mvDic_h.update({'flair': t2z_h})
    if "qt1" in featList_hipp:
        qt1z_h = matrix("hipp", '_space-hipp_hemi-{}_midthickness_desc-qt1_2mm.func.gii', TBL)
        mv_h.append(qt1z_h.values)
        mvDic_h.update({'qt1': qt1z_h})
    if "adc" in featList_hipp:
        adcz_h = matrix("hipp", '_space-hipp_hemi-{}_midthickness_desc-ADC_2mm.func.gii', TBL)
        mv_h.append(adcz_h.values)
        mvDic_h.update({'adc': adcz_h})
    if "thickness" in featList_hipp:
        ctz_h = matrix("hipp", '_space-hipp_hemi-{}_desc-thickness_2mm.func.gii', TBL)
        ctz_h = ctz_h * -1
        mv_h.append(ctz_h.values)
        mvDic_h.update({'thickness': ctz_h})


    if "flair" in featList_hipp:
        t2z_h = matrix("hipp", '/anat/surfaces/flair/', '_hemi-{}_space-flair_desc-flair_N4_den-0p5mm_label-hipp_midthickness_10mm.func.gii', TBL)
        mv_h.append(t2z_h.values)
        mvDic_h.update({'flair': t2z_h})
    if "qt1" in featList_hipp:
        qt1z_h = matrix("hipp", '/anat/surfaces/qt1/', '_hemi-{}_space-qt1_desc-qt1_den-0p5mm_label-hipp_midthickness_10mm.func.gii', TBL)
        mv_h.append(qt1z_h.values)
        mvDic_h.update({'qt1': qt1z_h})
    if "adc" in featList_hipp:
        adcz_h = matrix("hipp", '/dwi/surfaces/', '_hemi-{}_space-dwi_desc-dwi-ADC_den-0p5mm_label-hipp_midthickness_10mm.func.gii', TBL)
        mv_h.append(adcz_h.values)
        mvDic_h.update({'adc': adcz_h})
    if "thickness" in featList_hipp:
        ctz_h = matrix("hipp", '/surf/', '_hemi-{}_space-T1w_den-0p5mm_label-hipp_thickness.shape.gii', TBL)
        ctz_h = ctz_h * -1
        mv_h.append(ctz_h.values)
        mvDic_h.update({'thickness': ctz_h})

    if mv_h:
        mv_h = np.nanmean(mv_h, axis=0)
        mv_h_unthr = np.tile(mv_h, (1, 2))

        # Save unthrehsolded map, then filter out noise or non-significant findings and save
        np.savetxt(os.path.join(os.path.dirname(out), "z-brains", "asymmetry", "allSubjects_unthr_hipp-z.csv"),
                   mv_h, delimiter=",")
        mv_h[(mv_h > -thr) & (mv_h < thr)] = 0
        mv_h = np.tile(mv_h, (1, 2))
        np.savetxt(os.path.join(os.path.dirname(out), "z-brains", "asymmetry", "allSubjects_hipp-z.csv"),
                   mv_h, delimiter=",")

        # Save individual features (unthresholded and thresholded)
        for (_, mvName) in enumerate(mvDic_h):
            mv = mvDic_h[mvName]
            np.savetxt(os.path.join(os.path.dirname(out), "z-brains", "asymmetry", "allSubjects_unthr_hipp-{}.csv".
                                    format(str(mvName))), mv, delimiter=",")
            mv[(mv > -thr) & (mv < thr)] = 0
            mv = np.tile(mv, (1, 2))
            np.savetxt(os.path.join(os.path.dirname(out), "z-brains", "asymmetry", "allSubjects_hipp-{}.csv".
                                    format(str(mvName))), mv, delimiter=",")
    
    
    
    # ==================================================================
    # Plot and save multivariate cortical, subcortical, and hippocampal asymmetry
    # ==================================================================
    # Get row number of subject
    rn = (TBL.loc[TBL['ID'] == subject].index).tolist()

    # Save and plot multivariate z-score | cortical
    fname = os.path.join(os.path.dirname(out), "z-brains", "asymmetry", subject, session, subject + "_ctx-mz-unthr_{}.gii")
    mv_c_half = mv_c_unthr[rn, :]
    mv_h_lh = new_gifti_image(mv_c_half[:len(mv_c_half)//2])
    mv_h_rh = new_gifti_image(mv_c_half[len(mv_c_half)//2:])
    nib.save(mv_h_lh, fname.format('lh'))
    nib.save(mv_h_rh, fname.format('rh'))
    
    fname = os.path.join(os.path.dirname(out), "z-brains", "asymmetry", subject, session, subject + "_ctx-mz_{}.gii")
    mv_c_half = mv_c[rn, :]
    mv_h_lh = new_gifti_image(mv_c_half[:len(mv_c_half)//2])
    mv_h_rh = new_gifti_image(mv_c_half[len(mv_c_half)//2:])
    nib.save(mv_h_lh, fname.format('lh'))
    nib.save(mv_h_rh, fname.format('rh'))

    fname = os.path.join(os.path.dirname(out), "z-brains", "asymmetry", subject, session, subject + "_ctx-mz.png")
    mv_c_half[:, :len(mv_c_half.flatten())//2] = np.nan
    plot_cortical(array_name=mv_c_half, surface_name='conte69', size=(800, 180), zoom=1.18, cmap='RdBu_r',
                  color_bar=True, color_range=(-3, 3), screenshot=True, view=['lateral', 'medial', 'medial', 'lateral'],
                  filename=fname)

    # Save and plot multivariate z-score | subcortical
    fname = os.path.join(os.path.dirname(out), "z-brains", "asymmetry", subject, session, subject + "_sctx-mz-unthr_{}.gii")
    mv_s_half = mv_s_unthr[rn, :].flatten()
    mv_h_lh = new_gifti_image(mv_s_half[:len(mv_s_half)//2])
    mv_h_rh = new_gifti_image(mv_s_half[len(mv_s_half)//2:])
    nib.save(mv_h_lh, fname.format('lh'))
    nib.save(mv_h_rh, fname.format('rh'))
    
    fname = os.path.join(os.path.dirname(out), "z-brains", "asymmetry", subject, session, subject + "_sctx-mz_{}.gii")
    mv_s_half = mv_s[rn, :].flatten()
    mv_h_lh = new_gifti_image(mv_s_half[:len(mv_s_half)//2])
    mv_h_rh = new_gifti_image(mv_s_half[len(mv_s_half)//2:])
    nib.save(mv_h_lh, fname.format('lh'))
    nib.save(mv_h_rh, fname.format('rh'))
    
    mv_s_half[:len(mv_s_half)//2] = np.nan
    fname = os.path.join(os.path.dirname(out), "z-brains", "asymmetry", subject, session, subject + "_sctx-mz.png")
    plot_subcortical(array_name=mv_s_half, ventricles=False, size=(800, 180), zoom=1.18, cmap='RdBu_r',
                  color_bar=True, color_range=(-3, 3), screenshot=True, view=['lateral', 'medial', 'medial', 'lateral'],
                  filename=fname)

    # Save and plot multivariate z-score | hippocampal unthresholded
    fname = os.path.join(os.path.dirname(out), "z-brains", "asymmetry", subject, session, subject + "_hipp-mz-unthr_{}.gii")
    mv_h_half = mv_h_unthr[rn, :].flatten()
    mv_h_lh = new_gifti_image(mv_h_half[:len(mv_h_half)//2])
    mv_h_rh = new_gifti_image(mv_h_half[len(mv_h_half)//2:])
    nib.save(mv_h_lh, fname.format('lh'))
    nib.save(mv_h_rh, fname.format('rh'))

    # Save and plot multivariate z-score | hippocampal thresholded
    fname = os.path.join(os.path.dirname(out), "z-brains", "asymmetry", subject, session, subject + "_hipp-mz_{}.gii")
    mv_h_half = mv_h[rn, :].flatten()
    mv_h_lh = new_gifti_image(mv_h_half[:len(mv_h_half)//2])
    mv_h_rh = new_gifti_image(mv_h_half[len(mv_h_half)//2:])
    nib.save(mv_h_lh, fname.format('lh'))
    nib.save(mv_h_rh, fname.format('rh'))

    mv_h_half[:len(mv_h_half) // 2] = np.nan
    fname = os.path.join(os.path.dirname(out), "z-brains", "asymmetry", subject, session, subject + "_hipp-mz.png")
    plot_hippocampal(array_name=mv_h_half, size=(800, 180), zoom=1.18, cmap='RdBu_r',
                     color_bar=True, color_range=(-3, 3), screenshot=True,
                     view=['lateral', 'medial', 'ventral', 'dorsal'], filename=fname)
    
    
    
    # ================================================================
    # Plot and save univariate cortical, subcortical, and hippocampal asymmetry
    # ================================================================
    # Plot univariate z-score | cortical
    for (_, mvName) in enumerate(mvDic_c):
        fname = os.path.join(os.path.dirname(out), "z-brains", "asymmetry", subject, session, subject + "_ctx-{}-unthr_{}.gii")
        mv_tmp = np.loadtxt(os.path.join(os.path.dirname(out), "z-brains", "asymmetry", "allSubjects_ctx-{}.csv".
                                    format(str(mvName))), delimiter=",")[rn, :]
       
        mv_h_lh = new_gifti_image(mv_tmp[:len(mv_tmp)//2])
        mv_h_rh = new_gifti_image(mv_tmp[len(mv_tmp)//2:])
        nib.save(mv_h_lh, fname.format(str(mvName), 'lh'))
        nib.save(mv_h_rh, fname.format(str(mvName), 'rh'))

        mv_tmp[:, :len(mv_c_half.flatten())//2] = np.nan
        plot_cortical(array_name=mv_tmp, surface_name='conte69', size=(800, 180), zoom=1.18, cmap='RdBu_r',
                  color_bar=True, color_range=(-3, 3), screenshot=True, view=['lateral', 'medial', 'medial', 'lateral'],
                  filename=os.path.join(os.path.dirname(out), "z-brains", "asymmetry", subject, session, subject + "_ctx-{}.png".
                  format(str(mvName))))      

    # Plot univariate z-score | subcortical
    for (_, mvName) in enumerate(mvDic_s):
        fname = os.path.join(os.path.dirname(out), "z-brains", "asymmetry", subject, session, subject + "_sctx-{}-unthr.gii")
        mv_tmp = np.loadtxt(os.path.join(os.path.dirname(out), "z-brains", "asymmetry", "allSubjects_sctx-{}.csv".
                                    format(str(mvName))), delimiter=",")[rn, :].flatten()
        
        mv_h_lh = new_gifti_image(mv_tmp[:len(mv_tmp)//2])
        mv_h_rh = new_gifti_image(mv_tmp[len(mv_tmp)//2:])
        nib.save(mv_h_lh, fname.format(str(mvName), 'lh'))
        nib.save(mv_h_rh, fname.format(str(mvName), 'rh'))

        mv_tmp[:len(mv_tmp)//2] = np.nan
        plot_subcortical(array_name=mv_tmp, ventricles=False, size=(800, 180), zoom=1.18, cmap='RdBu_r',
                  color_bar=True, color_range=(-3, 3), screenshot=True, view=['lateral', 'medial', 'medial', 'lateral'],
                  filename=os.path.join(os.path.dirname(out), "z-brains", "asymmetry", subject,
                                        session, subject + "_sctx-{}.png".format(str(mvName))))

    # Plot univariate z-score | hippocampal
    for (_, mvName) in enumerate(mvDic_h):
        fname = os.path.join(os.path.dirname(out), "z-brains", "asymmetry", subject, session, subject + "_hipp-{}-unthr_{}.gii")
        mv_tmp = np.loadtxt(os.path.join(os.path.dirname(out), "z-brains", "asymmetry", "allSubjects_hipp-{}.csv".
                                    format(str(mvName))), delimiter=",")[rn, :].flatten()
        
        mv_h_lh = new_gifti_image(mv_tmp[:len(mv_tmp)//2])
        mv_h_rh = new_gifti_image(mv_tmp[len(mv_tmp)//2:])
        nib.save(mv_h_lh, fname.format(str(mvName), 'lh'))
        nib.save(mv_h_rh, fname.format(str(mvName), 'rh'))

        mv_tmp[:len(mv_tmp) // 2] = np.nan
        plot_hippocampal(array_name=mv_tmp, size=(800, 180), zoom=1.18, cmap='RdBu_r',
                     color_bar=True, color_range=(-3, 3), screenshot=True,
                     view=['lateral', 'medial', 'ventral', 'dorsal'],
                     filename=os.path.join(os.path.dirname(out), "z-brains", "asymmetry", subject,
                                           session, subject + "_hipp-{}.png".format(str(mvName))))
        
