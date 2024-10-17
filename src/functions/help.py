import colorama

pcolor = colorama.Fore.MAGENTA  # Purple
rcolor = colorama.Fore.RED  # Red
gcolor = colorama.Fore.GREEN  # Green
bcolor = colorama.Fore.BLUE  # Blue
gray = colorama.Fore.LIGHTBLACK_EX  # Gray
nc = colorama.Style.RESET_ALL  # No color

help = f"""
{pcolor}COMMAND:{nc}
zbrains.py


{pcolor}OPTIONS:{nc}
\t{rcolor}--sub{nc} ID                  : Subject ID. This is the target subject. Example: 'sub-PX001'. Also accepts "all".
\t{rcolor}--dataset{nc} path            : Path to the BIDS dataset containing the target subject's data.
                                    Example: '/path/to/BIDSDataset'.
\t{rcolor}--zbrains{nc} dir             : Name of the zbrains derivative folder in the target BIDS dataset. The
                                    folder will be created if it does not exist. Example: '--zbrains zbrainsdir' for
                                    '/path/to/BIDSDataset/derivatives/zbrainsdir'.
\t{gcolor}--run{nc} task                : Tasks to perform. Options:
                                    - {bcolor}proc{nc}          : perform post-processing of target subject (default).
                                    - analysis      : perform analysis (regional & asymmetry) and generate clinical
                                                        report for target subject. The analysis is performed wrt the
                                                        reference subjects (those provided with --demo_ref).
                                                        Post-processing of both target and reference subjects must be
                                                        performed beforehand.
                                    - all           : perform all tasks

\t{gcolor}--ses{nc} [ses]               : Identifier for a session in the target subject data. If omitted, all sessions will be used. Example: 'ses-001'.
\t{gcolor}--micapipe{nc} [dir]          : Name of the micapipe derivative folder in the target BIDS dataset. Required
                                    only for post-processing. Example: '--micapipe micapipedir' for
                                    '/path/to/BIDSDataset/derivatives/micapipedir'.
\t{gcolor}--hippunfold{nc} [dir]        : Name of the hippunfold derivative folder in the target BIDS dataset. Required
                                    only for post-processing. Example: '--hippunfold hipdir' for
                                    '/path/to/BIDSDataset/derivatives/hipdir'.
\t{gcolor}--plugin{nc} [dir]           : Name of a plugin derivative folder in the target BIDS dataset. zbrains can accept
                                data outside of micapipe and hippunfold as a 'plugin' folder. However, these data MUST
                                be formatted as BIDS-derivatives exactly as in micapipe and hippunfold. If hippocampal
                                surface data are present then they will be used but otherwise volumetric data will be
                                mapped to hippocampal and subcortical surfaces. 
                                '/path/to/BIDSDataset/derivatives/plugindir'.
\t{gcolor}--demo{nc} [path]             : CSV/TSV file with demographics for the target subject. Required only for
                                    analysis when provided by --normative or --deconfound. Additionally, the file is
                                    also used to extract the subject's age and sex, if available, to be included in the
                                    clinical report. The file must include one row with target subject ID and session.
                                    Expected column names:
                                    - participant_id: Subject ID. Required.
                                    - session_id    : Session ID. Use 'n/a' if no session available.
                                    - age           : Subject age. Required only when used by --normative or
                                                        --deconfound. If provided, it will also be included in the
                                                        clinical report.
                                    - sex           : Subject sex. Possible values: 'F' or 'M'. Required only when
                                                        used by --normative or --deconfound. If provided, it will
                                                        also be included in the clinical report.
                                    - site          : Acquisition site. Required only when used by --normative or
                                                        --deconfound.
                                    - Other         : Any other columns used by --normative or --deconfound.
                                    Use the --column_map option to indicate if these variables are under different
                                    column names in your file.

\t{gcolor}--dataset_ref{nc} [path ...]  : Paths to the BIDS datasets containing the reference subjects data. Required
                                    only for analysis. Each dataset must correspond to one file in --demo_ref.
\t{gcolor}--zbrains_ref{nc} [dir ...]   : Names of the zbrains derivative folder in each of the reference datasets.
                                    Required only for analysis. If only one folder name is provided but there are
                                    multiple datasets, we assume the same name in all datasets.
\t{gcolor}--demo_ref{nc} [path ...]     : CSV/TSV files with demographics for reference subjects. Required only for
                                    analysis. There must be one file for each reference dataset (--dataset_ref).
                                    Required only for analysis. See --demo for expected column names.

\t{gcolor}--struct{nc} [structure ...]  : Structures to use in processing and/or analysis. Options:
                                    - cortex        : cortical data
                                    - subcortex     : subcortical data
                                    - hippocampus   : hippocampal data
                                    - {bcolor}all{nc}           : all structures (default)
\t{gcolor}--feat{nc} [feature ...]      : Features to use in processing and/or analysis. Options:
                                    - ADC           : apparent diffusion coefficient
                                    - FA            : fractional anisotropy
                                    - flair         : FLAIR
                                    - qT1           : quantitative T1
                                    - thickness     : cortical thickness (for subcortex, volume is used)
                                    - {bcolor}all{nc}           : all features (default)
                                    - plugin-*      : when pulling data from a plugin, feature names must be given the 
                                                    'plugin-' prefix (but this is not needed in the actual file name)
\t{gcolor}--normative{nc} [cov ...]     : Normative modeling based on provided covariates. Covariates must match
                                    columns in --demo and --demo_ref files. Note that --normative expects some
                                    covariates to have specific names (see --column_map).
                                    Example: '--normative site age sex'.
\t{gcolor}--deconfound{nc} [[-]cov ...] : Deconfounding based on provided covariates. Covariates must match columns in
                                    --demo and --demo_ref CSV files. If the covariates include 'site', deconfounding is
                                    performed using ComBat. Otherwise, linear regression is used to regress out the
                                    effects of the covariates from the data. By default, ComBat preserves additional
                                    covariates. To remove the effects of a covariate, prepend with '-' (ignored when not
                                    using ComBat). Note that --deconfound expects some covariates to have specific names
                                    (see --column_map). Example: '--deconfound site -age -sex group' to harmonize data
                                    while preserving the effects of group and removing those of age and sex.
\t{gcolor}--resolution{nc} [res ...]    : Surface resolutions to use for cortex and hippocampus. Options:
                                    - low           : 5k cortical & 2mm hippocampal surfaces
                                    - high          : 32k cortical surfaces & 0p5mm hippocampal surfaces
                                    - {bcolor}all{nc}           : all resolutions (default)
\t{gcolor}--label_ctx{nc} [label]       : Cortical surfaces used in the volume to surface mapping. Options:
                                    - {bcolor}white{nc}         : WM surface (default)
                                    - midthickness  : Midthickness surface
                                    - pial          : Pial surface
                                    - swmD          : Superficial white matter, where D indicates the distance in
                                                        millimeters. Example: --label_ctx swm2
\t{gcolor}--label_hip{nc} [label]       : Hippocampal surface used in the volume to surface mapping. Options:
                                    - {bcolor}midthickness{nc}  : Midthickness surface (default)
\t{gcolor}--smooth_ctx{nc} [size]       : Size of gaussian smoothing kernel in mm used for cortical features.
                                    Default is {bcolor}5{nc}.
\t{gcolor}--smooth_hip{nc} [size]       : Size of gaussian smoothing kernel in mm used for hippocampal features.
                                    Default is {bcolor}2{nc}.
\t{gcolor}--threshold{nc} [th]          : Threshold for statistical maps used for clinical reports.
                                    Default is {bcolor}1.96{nc}.
\t{gcolor}--dicoms{nc} [bool]            : If set to 1, will convert NIfTI files to DICOMs during analysis. This can be time-intensive. 
                                    Default is {bcolor}1{nc}.
\t{gcolor}--volumetric{nc} [bool]            : If set to 1, will convert surface maps to volume space using ribbon-constrained mapping. 
                                    Default is {bcolor}1{nc}.
\t{gcolor}--column_map{nc} [VAR=col ...]: Map expected to actual column names in the CSV/TSV files:
                                    - participant_id: Subject ID is assumed to be provided by the 'participant_id'
                                                        column, unless indicated otherwise. For example, if subject ID
                                                        is under the column \u2018SubID\u2019, you can indicate this with
                                                        --column_map participant_id=SubID.
                                    - session_id    : Session ID is assumed to be provided by the \u2018session_id\u2019 column,
                                                        unless indicated otherwise (e.g., --column_map session_id=ses)
                                    - age           : Age is assumed to be provided by the \u2018age\u2019 column, unless
                                                        indicated otherwise (e.g., --column_map age=\"Subject age\")
                                    - sex           : Sex is assumed to be provided by the 'sex' column, unless
                                                        indicated otherwise (e.g., --column_map ses=\"Subject sex\")
                                    - site          : Acquisition site is assumed to be provided by the \u2018site\u2019 column,
                                                        unless indicated otherwise (e.g., --column_map site=center)                 
\t{gcolor}--n_jobs{nc} [number]         : Number of jobs to run in parallel. Default is {bcolor}1{nc}.
\t{gcolor}--wb_path{nc} [path]          : Path to the Connectome Workbench binaries. Default is {bcolor}/data/mica1/01_programs/workbench-1.4.2/bin_linux64{nc}.
\t{gcolor}--patient_prefix{nc} [prefix] : Prefix to use when determining patients versus controls. Default is {bcolor}PX{nc}.
\t{gcolor}--delete_temps{nc} [bool]     : If set to True, will delete any ragged temp files left from crashed analyses, then exit. Default is {bcolor}False{nc}.
\t{gcolor}--verbose{nc} [level]         : Verbosity level (default is {bcolor}-1{nc}). Levels:
                                    - 0             : Only errors
                                    - 1             : Warning messages and previous levels
                                    - 2             : Information messages and previous levels
                                    - 3             : Command logs and previous levels
                                    - >3 or <0      : All messages
\t{gcolor}--help{nc}                    : Print help
\t{gcolor}--version{nc}                 : Print software version
\t{gcolor}--pyinit{nc}                  : Specify a Python source, (e.g. a conda environment) to activate 


{pcolor}USAGE:{nc}
    {gray}# Post-processing{nc}
    {pcolor}zbrains.py{nc} {gcolor}--run{nc} proc
            {rcolor}--sub{nc} <participant_id>
            {gcolor}--ses{nc} <session_id>
            {rcolor}--dataset{nc} <path_bids_dataset>
            {rcolor}--zbrains{nc} <zbrains_directory>
            {gcolor}--micapipe{nc} <micapipe_directory>
            {gcolor}--hippunfold{nc} <hipunfold_directory>

    {gray}# Analysis{nc}
    {pcolor}zbrains.py{nc} {gcolor}--run{nc} analysis
            {rcolor}--sub{nc} <participant_id>
            {gcolor}--ses{nc} <session_id>
            {rcolor}--dataset{nc} <participant_dataset>
            {rcolor}--zbrains{nc} <participant_zbrains_dir>
            {gcolor}--demo_ref{nc} <reference_subjects1.csv> <reference_subjects2.csv>
            {gcolor}--dataset_ref{nc} <reference_dataset1> <reference_dataset2>
            {gcolor}--zbrains_ref{nc} <reference_zbrains_dir1> <reference_zbrains_dir2>



{pcolor}DEPENDENCIES:{nc}
    > workbench   1.4.2   (https://www.humanconnectome.org/software/workbench-command)
    > ANTs        2.3.4   (https://github.com/ANTsX/ANTs)
    > python      3.10    (https://www.python.org)

    To customize binary locations, use the following environment variables:
    - Set ANTSPATH for ANTs
    - Set WORKBENCH_PATH for Workbench

    Control the number of threads:
    - Set ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS for ANTs
    - Set OMP_NUM_THREADS for Workbench

    Example:
    {gray}# Set threads for ANTs{nc}
    $ export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=4

    {gray}# Set threads for Workbench{nc}
    $ export OMP_NUM_THREADS=4





McGill University, MNI, MICA lab, April 2023
https://github.com/MICA-MNI/micapipe
https://github.com/MICA-MNI/z-brains
http://mica-mni.github.io/
"""
