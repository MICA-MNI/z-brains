import numpy as np
import nibabel as nib
import glob
import sys
import pandas as pd
import warnings

controls_list = sys.argv[1]
patient = sys.argv[2]
featureStr = sys.argv[3]
smoothCtx = sys.argv[4]
smoothHipp = sys.argv[5]
outDir = sys.argv[6]
ses = sys.argv[7]


if featureStr == 'all':
    featureList = ['ADC','FA','T1map','thickness']
else:
    featureList = featureStr.split(',')

subList = pd.read_csv("/data_/mica1/03_projects/jessica/hackathon2023/lists/control.csv")


################# cortex #################


### get control distribution ###
# load example for shaping
example = nib.load(f'{outDir}/sub-{patient}/{ses}/maps/cortex/sub-{patient}_{ses}_hemi-L_feature-thickness_smooth-{smoothCtx}.func.gii').darrays[0].data
# initialize
controlData = np.ones((len(example),len(featureList),2,len(subList)))*np.nan
# load all control data
for f,feature in enumerate(featureList):
    for s,subject in enumerate(subList['ID']):
        for h,hemi in enumerate(['L','R']):
            session = f'ses-0' + str(subList["SES"][s])
            try:
                controlData[:,f,h,s] = nib.load(f'{outDir}/{subject}/{session}/maps/cortex/{subject}_{session}_hemi-{hemi}_feature-{feature}_smooth-{smoothCtx}.func.gii').darrays[0].data
            except:
                warnings.warn(f' data not found: {outDir}/{subject}/{session}/maps/cortex/{subject}_{session}_hemi-{hemi}_feature-{feature}_smooth-{smoothCtx}.func.gii')
controlMean = np.nanmean(controlData,axis=3)
controlStd = np.nanstd(controlData,axis=3)

### compare current subject ###
patientData = np.ones((len(example),len(featureList),2))*np.nan
for f,feature in enumerate(featureList):
    for h,hemi in enumerate(['L','R']):
        patientData[:,f,h] = nib.load(f'{outDir}/sub-{patient}/{ses}/maps/cortex/sub-{patient}_{ses}_hemi-{hemi}_feature-{feature}_smooth-{smoothCtx}.func.gii').darrays[0].data

zscores = patientData - controlMean
zscores = zscores / controlStd

### save ###
for f,feature in enumerate(featureList):
    for h,hemi in enumerate(['L','R']):
        data_array = nib.gifti.GiftiDataArray(data=zscores[:,f,h])
        image = nib.gifti.GiftiImage()
        image.add_gifti_data_array(data_array)
        nib.save(image, f'{outDir}/sub-{patient}/{ses}/norm-z/cortex/sub-{patient}_{ses}_hemi-{hemi}_feature-{feature}_smooth-{smoothCtx}.func.gii')


################# hippocampus #################


### get control distribution ###
# load example for shaping
example = nib.load(f'{outDir}/sub-{patient}/{ses}/maps/hippocampus/sub-{patient}_{ses}_hemi-L_feature-thickness_smooth-{smoothHipp}.func.gii').darrays[0].data
# initialize
controlData = np.ones((len(example),len(featureList),2,len(subList)))*np.nan
# load all control data
for f,feature in enumerate(featureList):
    for s,subject in enumerate(subList['ID']):
        for h,hemi in enumerate(['L','R']):
            session = f'ses-0' + str(subList["SES"][s])
            try:
                controlData[:,f,h,s] = nib.load(f'{outDir}/{subject}/{session}/maps/hippocampus/{subject}_{session}_hemi-{hemi}_feature-{feature}_smooth-{smoothHipp}.func.gii').darrays[0].data
            except:
                warnings.warn(f' data not found: {outDir}/{subject}/{session}/maps/hippocampus/{subject}_{session}_hemi-{hemi}_feature-{feature}_smooth-{smoothHipp}.func.gii')
controlMean = np.nanmean(controlData,axis=3)
controlStd = np.nanstd(controlData,axis=3)

### compare current subject ###
patientData = np.ones((len(example),len(featureList),2))*np.nan
for f,feature in enumerate(featureList):
    for h,hemi in enumerate(['L','R']):
        patientData[:,f,h] = nib.load(f'{outDir}/sub-{patient}/{ses}/maps/hippocampus/sub-{patient}_{ses}_hemi-{hemi}_feature-{feature}_smooth-{smoothHipp}.func.gii').darrays[0].data

zscores = patientData - controlMean
zscores = zscores / controlStd

### save ###
for f,feature in enumerate(featureList):
    for h,hemi in enumerate(['L','R']):
        data_array = nib.gifti.GiftiDataArray(data=zscores[:,f,h])
        image = nib.gifti.GiftiImage()
        image.add_gifti_data_array(data_array)
        nib.save(image, f'{outDir}/sub-{patient}/{ses}/norm-z/hippocampus/sub-{patient}_{ses}_hemi-{hemi}_feature-{feature}_smooth-{smoothHipp}.func.gii')


################# subcortex #################


featureList = list(map(lambda x: x.replace('thickness', 'volume'), featureList))

### get control distribution ###
# load example for shaping
example = pd.read_csv(f'{outDir}/sub-{patient}/{ses}/maps/subcortex/sub-{patient}_{ses}_feature-volume.csv')
print(example)
# initialize
controlData = np.ones((len(example.columns[1:]),len(featureList),len(subList)))*np.nan
# load all control data
for f,feature in enumerate(featureList):
    for s,subject in enumerate(subList['ID']):
            session = f'ses-0' + str(subList["SES"][s])
            try:
                controlData[:,f,s] = pd.read_csv(f'{outDir}/{subject}/{session}/maps/subcortex/{subject}_{session}_feature-{feature}.csv',header=[0],index_col=0).to_numpy()
            except:
                warnings.warn(f' data not found: {outDir}/{subject}/{session}/maps/subcortex/{subject}_{session}_feature-{feature}.csv')
controlMean = np.nanmean(controlData,axis=2)
controlStd = np.nanstd(controlData,axis=2)

### compare current subject ###
patientData = np.ones((len(example.columns[1:]),len(featureList)))*np.nan
for f,feature in enumerate(featureList):
    patientData[:,f] = pd.read_csv(f'{outDir}/sub-{patient}/{ses}/maps/subcortex/sub-{patient}_{ses}_feature-{feature}.csv',header=[0],index_col=0).to_numpy()

zscores = patientData - controlMean
zscores = zscores / controlStd

### save ###
for f,feature in enumerate(featureList):
    for i in example.columns[1:]:
        example[i] = zscores[[i,f,h]]
    example.to_csv(f'{outDir}/sub-{patient}/{ses}/norm-z/subcortex/sub-{patient}_{ses}_feature-{feature}.csv')


