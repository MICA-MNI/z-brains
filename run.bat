setlocal enabledelayedexpansion
set "pth_dataset=E:\data"

set "micapipe_dir=micapipe"
set "hippunfold_dir=hippunfold"
set "zbrains_dir=zbrains_newtest"

set "demo_controls=E:\participants_mics_hc.csv"
set "demo_patients=E:\PX_participants.csv"

set "sids=sub-PX103"
set "sess=all"
call conda activate zbrains
call python -m src.zbrains --run "analysis" ^
                   --sub "%sids%" ^
                   --micapipe %micapipe_dir% ^
                   --hippunfold %hippunfold_dir% ^
                   --dataset %pth_dataset% ^
                   --zbrains %zbrains_dir% ^
                   --demo_ref %demo_controls% ^
                   --demo %demo_patients% ^
                   --dataset_ref %pth_dataset% ^
                   --zbrains_ref %zbrains_dir% ^
                   --smooth_ctx 10 ^
                   --smooth_hip 5 ^
                   --n_jobs 4 ^
                   --n_jobs_wb 4 ^
                   --dicoms 0 ^
		           --label_ctx "white" ^
                   --wb_path "C:/Users/Ian/Downloads/workbench-windows64-v1.5.0/workbench/bin_windows64" ^
                   --column_map participant_id=ID session_id=SES ^
                   --verbose 2
    )
pause