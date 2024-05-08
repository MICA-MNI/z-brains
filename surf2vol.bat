@REM wb_command -metric-to-volume-mapping^
@REM  "C:\Users\Ian\Documents\GitHub\z-brains-IanTesting\src\data\sub-PX001_ses-02_hemi-L_surf-fsLR-32k_label-midthickness_feature-ADC_smooth-10mm_analysis-asymmetry.func.gii"^
@REM  "E:\data\derivatives\micapipe\sub-PX001\ses-02\surf\sub-PX001_ses-02_hemi-L_space-nativepro_surf-fsLR-32k_label-midthickness.surf.gii"^
@REM  "C:\Users\Ian\Documents\GitHub\z-brains-IanTesting\src\data\mni152.nii.gz"^
@REM  "C:\Users\Ian\Documents\GitHub\z-brains-IanTesting\src\data\output_weirdbrain.nii.gz"^
@REM  -ribbon-constrained^
@REM  "E:\data\derivatives\micapipe\sub-PX001\ses-02\surf\sub-PX001_ses-02_hemi-L_space-nativepro_surf-fsLR-32k_label-white.surf.gii"^
@REM  "E:\data\derivatives\micapipe\sub-PX001\ses-02\surf\sub-PX001_ses-02_hemi-L_space-nativepro_surf-fsLR-32k_label-pial.surf.gii"

@REM pause


@REM wb_command -surface-apply-affine^
@REM  "E:\data\derivatives\micapipe\sub-PX001\ses-02\surf\sub-PX001_ses-02_hemi-L_space-nativepro_surf-fsLR-32k_label-midthickness.surf.gii"^
@REM  "E:\data\derivatives\micapipe\sub-PX001\ses-02\xfm\sub-PX001_ses-02_from-nativepro_brain_to-MNI152_0.8mm_mode-image_desc-SyN_0GenericAffine.mat"^
@REM  "C:\Users\Ian\Documents\GitHub\z-brains-IanTesting\src\data\mni152_corrected-midthickness.nii.gz"

@REM wb_command -surface-apply-affine^
@REM  "E:\data\derivatives\micapipe\sub-PX001\ses-02\surf\sub-PX001_ses-02_hemi-L_space-nativepro_surf-fsLR-32k_label-white.surf.gii"^
@REM  "E:\data\derivatives\micapipe\sub-PX001\ses-02\xfm\sub-PX001_ses-02_from-nativepro_brain_to-MNI152_0.8mm_mode-image_desc-SyN_0GenericAffine.mat"^
@REM  "C:\Users\Ian\Documents\GitHub\z-brains-IanTesting\src\data\mni152_corrected-white.nii.gz"

@REM wb_command -surface-apply-affine^
@REM  "E:\data\derivatives\micapipe\sub-PX001\ses-02\surf\sub-PX001_ses-02_hemi-L_space-nativepro_surf-fsLR-32k_label-pial.surf.gii"^
@REM  "E:\data\derivatives\micapipe\sub-PX001\ses-02\xfm\sub-PX001_ses-02_from-nativepro_brain_to-MNI152_0.8mm_mode-image_desc-SyN_0GenericAffine.mat"^
@REM  "C:\Users\Ian\Documents\GitHub\z-brains-IanTesting\src\data\mni152_corrected-pial.nii.gz"



wb_command -metric-to-volume-mapping^
 "C:\Users\Ian\Documents\GitHub\z-brains-IanTesting\src\data\sub-PX001_ses-02_hemi-L_surf-fsLR-32k_label-midthickness_feature-ADC_smooth-10mm_analysis-asymmetry.func.gii"^
 "E:\data\derivatives\micapipe\sub-PX001\ses-02\surf\sub-PX001_ses-02_hemi-L_space-nativepro_surf-fsLR-32k_label-midthickness.surf.gii"^
 "E:\data\derivatives\micapipe\sub-PX001\ses-02\maps\sub-PX001_ses-02_space-nativepro_model-DTI_map-ADC.nii.gz"^
 "C:\Users\Ian\Documents\GitHub\z-brains-IanTesting\src\data\output_okbrain.nii.gz"^
 -ribbon-constrained^
 "E:\data\derivatives\micapipe\sub-PX001\ses-02\surf\sub-PX001_ses-02_hemi-L_space-nativepro_surf-fsLR-32k_label-white.surf.gii"^
 "E:\data\derivatives\micapipe\sub-PX001\ses-02\surf\sub-PX001_ses-02_hemi-L_space-nativepro_surf-fsLR-32k_label-pial.surf.gii"

pause