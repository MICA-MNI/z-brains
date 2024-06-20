@REM @echo off
setlocal enabledelayedexpansion
set "pth_dataset=E:\data"

set "micapipe_dir=micapipe"
set "hippunfold_dir=hippunfold"
set "zbrains_dir=zbrains_Ian_testing"

set "demo_controls=E:\HC_participants_full.csv"
set "demo_patients=E:\PX_participants.csv"

set "sids=all"
set "sess=all"
@REM set "sids=all"
@REM set "sess=all"
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
                   --wb_path "C:/Users/Ian/Downloads/workbench-windows64-v1.5.0/workbench/bin_windows64" ^
                   --column_map participant_id=ID session_id=SES ^
                   --verbose 2
    )
@REM for %%s in (%sids%) do (
@REM @REM    echo %%s
   
@REM    set "sess="
@REM    for /D %%d in ("%pth_dataset%\derivatives\%micapipe_dir%\%%s\*") do (
@REM     @REM echo %%d
@REM       set "part2=%%~nxd"
@REM       set "ses_part2=!part2:ses-=!"
@REM       set "sess=!sess! ses-!ses_part2!"
      
@REM    )
   
@REM    if not "!sess!"=="" (
@REM    for %%t in (!sess!) do (
@REM       call conda activate zbrains
@REM       @REM echo %%t
@REM       @REM echo !sess!
@REM    @REM     call python zbrains.py --run proc analysis ^
@REM    @REM                 --sub "%%s" ^
@REM    @REM                 --ses "%%t" ^
@REM    @REM                 --micapipe %micapipe_dir% ^
@REM    @REM                 --hippunfold %hippunfold_dir% ^
@REM    @REM                 --dataset %pth_dataset% ^
@REM    @REM                 --zbrains %zbrains_dir% ^
@REM    @REM                 --demo_ref %demo_controls% ^
@REM    @REM                 --dataset_ref %pth_dataset% ^
@REM    @REM                 --zbrains_ref %zbrains_dir% ^
@REM    @REM                 --n_jobs 2 ^
@REM    @REM                 --column_map participant_id=ID session_id=SES ^
@REM    @REM                 --verbose 2
@REM    @REM  )
@REM            call python zbrains.py --run proc ^
@REM                    --sub "all" ^
@REM                    --micapipe %micapipe_dir% ^
@REM                    --hippunfold %hippunfold_dir% ^
@REM                    --dataset %pth_dataset% ^
@REM                    --zbrains %zbrains_dir% ^
@REM                    --demo_ref %demo_controls% ^
@REM                    --dataset_ref %pth_dataset% ^
@REM                    --zbrains_ref %zbrains_dir% ^
@REM                    --n_jobs 1 ^
@REM                    --column_map participant_id=ID session_id=SES ^
@REM                    --verbose 2
@REM     )
@REM    @REM  call python zbrains_batch.py --run proc ^
@REM    @REM                 --sub "%%s" ^
@REM    @REM                 --ses "%%t" ^
@REM    @REM                 --micapipe %micapipe_dir% ^
@REM    @REM                 --hippunfold %hippunfold_dir% ^
@REM    @REM                 --dataset %pth_dataset% ^
@REM    @REM                 --zbrains %zbrains_dir% ^
@REM    @REM                 --demo %demo_controls% ^
@REM    @REM                 --demo_ref %demo_controls% ^
@REM    @REM                 --dataset_ref %pth_dataset% ^
@REM    @REM                 --zbrains_ref %zbrains_dir% ^
@REM    @REM                 --column_map participant_id=participant_id session_id=SES ^
@REM    @REM                 --verbose 2
@REM    @REM  )
@REM )
@REM )
pause