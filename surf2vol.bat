wb_command -metric-resample^
 "C:\Users\Ian\Documents\GitHub\z-brains-IanTesting\src\data\sub-PX001_ses-02_hemi-L_surf-fsLR-32k_label-midthickness_feature-FA_smooth-10mm_analysis-asymmetry.func.gii"^
 "C:\Users\Ian\Documents\GitHub\z-brains-IanTesting\src\data\templates\fsLR-32k.L.sphere.reg.surf.gii"^
 "E:\data\derivatives\micapipe\sub-PX001\ses-02\surf\sub-PX001_ses-02_hemi-L_surf-fsnative_label-sphere.surf.gii"^
 "BARYCENTRIC"^
 "C:\Users\Ian\Documents\GitHub\z-brains-IanTesting\src\data\fsnative_surf_L.func.gii"

wb_command -metric-to-volume-mapping^
 "C:\Users\Ian\Documents\GitHub\z-brains-IanTesting\src\data\fsnative_surf_L.func.gii"^
 "E:\data\derivatives\micapipe\sub-PX001\ses-02\surf\sub-PX001_ses-02_hemi-L_space-nativepro_surf-fsnative_label-midthickness.surf.gii"^
 "C:\Users\Ian\Documents\GitHub\z-brains-IanTesting\src\data\patient_surfs\sub-PX001_ses-02_space-nativepro_map-T1map.nii.gz"^
 "C:\Users\Ian\Documents\GitHub\z-brains-IanTesting\src\data\Nativepro_L.nii.gz"^
 -ribbon-constrained^
 "E:\data\derivatives\micapipe\sub-PX001\ses-02\surf\sub-PX001_ses-02_hemi-L_space-nativepro_surf-fsnative_label-white.surf.gii"^
 "E:\data\derivatives\micapipe\sub-PX001\ses-02\surf\sub-PX001_ses-02_hemi-L_space-nativepro_surf-fsnative_label-pial.surf.gii"

@REM wb_command -volume-resample^
@REM  "C:\Users\Ian\Documents\GitHub\z-brains-IanTesting\src\data\Nativepro_L.nii.gz"^
@REM  "C:\Users\Ian\Documents\GitHub\z-brains-IanTesting\src\data\templates\MNI152_T1_0.8mm_brain.nii.gz"^
@REM  ENCLOSING_VOXEL^
@REM  "C:\Users\Ian\Documents\GitHub\z-brains-IanTesting\src\data\warpthenaffine.nii.gz"^
@REM  -warp^
@REM  "E:\data\derivatives\micapipe\sub-PX001\ses-02\xfm\sub-PX001_ses-02_from-nativepro_brain_to-MNI152_0.8mm_mode-image_desc-SyN_1Warp.nii.gz"^
@REM  -affine^
@REM  "C:\Users\Ian\Documents\GitHub\z-brains-IanTesting\src\data\transforms\canonical_to_fsaverage_noninv.txt"

