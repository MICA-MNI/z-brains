@REM @echo off
setlocal enabledelayedexpansion
set "pth_dataset=E:\data"

set "micapipe_dir=micapipe"
set "hippunfold_dir=hippunfold"
set "zbrains_dir=zbrains_out"

set "demo_controls=E:\data\part.csv"

set "sids="
for /D %%d in ("%pth_dataset%\derivatives\%micapipe_dir%\*") do (
   set "part2=%%~nxd"
   set "sid_part2=!part2:sub-=!"
   set "sids=!sids! sub-!sid_part2!"
)

@REM echo "%sid%"
@REM echo "%sids%"

set "sids=sub-PX070 sub-PX003 sub-PX002 sub-PX001"
echo "here"
for %%s in (%sids%) do (
@REM    echo %%s
   
   set "sess="
   for /D %%d in ("%pth_dataset%\derivatives\%micapipe_dir%\%%s\*") do (
    @REM echo %%d
      set "part2=%%~nxd"
      set "ses_part2=!part2:ses-=!"
      set "sess=!sess! ses-!ses_part2!"
      
   )
   
   if not "!sess!"=="" (
   for %%t in (!sess!) do (
      call conda activate zbrains
      @REM echo %%t
      @REM echo !sess!
   @REM     call python zbrains.py --run proc analysis ^
   @REM                 --sub "%%s" ^
   @REM                 --ses "%%t" ^
   @REM                 --micapipe %micapipe_dir% ^
   @REM                 --hippunfold %hippunfold_dir% ^
   @REM                 --dataset %pth_dataset% ^
   @REM                 --zbrains %zbrains_dir% ^
   @REM                 --demo_ref %demo_controls% ^
   @REM                 --dataset_ref %pth_dataset% ^
   @REM                 --zbrains_ref %zbrains_dir% ^
   @REM                 --n_jobs 2 ^
   @REM                 --column_map participant_id=ID session_id=SES ^
   @REM                 --verbose 2
   @REM  )
           call python zbrains.py --run proc ^
                   --sub "all" ^
                   --micapipe %micapipe_dir% ^
                   --hippunfold %hippunfold_dir% ^
                   --dataset %pth_dataset% ^
                   --zbrains %zbrains_dir% ^
                   --demo_ref %demo_controls% ^
                   --dataset_ref %pth_dataset% ^
                   --zbrains_ref %zbrains_dir% ^
                   --n_jobs 1 ^
                   --column_map participant_id=ID session_id=SES ^
                   --verbose 2
    )
   @REM  call python zbrains_batch.py --run proc ^
   @REM                 --sub "%%s" ^
   @REM                 --ses "%%t" ^
   @REM                 --micapipe %micapipe_dir% ^
   @REM                 --hippunfold %hippunfold_dir% ^
   @REM                 --dataset %pth_dataset% ^
   @REM                 --zbrains %zbrains_dir% ^
   @REM                 --demo %demo_controls% ^
   @REM                 --demo_ref %demo_controls% ^
   @REM                 --dataset_ref %pth_dataset% ^
   @REM                 --zbrains_ref %zbrains_dir% ^
   @REM                 --column_map participant_id=participant_id session_id=SES ^
   @REM                 --verbose 2
   @REM  )
)
)