Python Conversion Changelog
- 
### **Functional Changes**
- "**--sub**" argument now accepts the keyword "**all**", which will attempt to run zbrains on all available subjects while processing, and while performing analysis, it will run on all subjects matching the pattern given in the argument (**--patient_prefix**). 
- "**--ses**" argument now also accepts the keyword "**all**", as well as now being optional. If the "**--ses**" argument is not given, it is assumed to be "**all**".
- "**--column_map**" argument has been given a default value of "**participant_id=ID session_id=SES**"
- "**--init**" argument has been removed, in favour of two arguments: "**--wb_path**" and "**--pyinit**"
- "**--wb_path**" specified the path to the workbench bin folder. It is set to "**/data/mica1/01_programs/workbench-1.4.2/bin_linux64**" by default
- "**--pyinit**" specifies the root folder of your desired Python environment. By default it is set to "**/data/mica1/01_programs/miniforge3**". It can be set to any conda-like environment manager that has an environment called "**zbrains**". Setting "**--pyinit=false**" will use the currently activated python environment.
- "**--patient_prefix**" specifies the prefix given to patient cases. Defaults to "**PX**"
- "**--delete_temps**" is a utility function which, if enabled, will not run any analyses, and will instead clean up any loose temp folders in the zbrains dataset derivates folder which may have been leftover from a system crash.

### **How to speed up your analyses with new zbrains**
- zbrain and zbrains_bash are now identical. Batch processing is handled by Python, instead of bash.
- the new "**--n_jobs**" argument will declare how many separate jobs you would like to run in parallel. Speed will increase with higher n_jobs, with larger gains seen when running zbrains on many subjects. Currently, speed increases in the analysis section are limited to **n_jobs**=5. Larger numbers of jobs will have no effect on speed.

### **Backend Changes**
- There is now a testing framework for zbrains. Automated tests will determine if the whole pipeline runs following backend changes. This backend is written with pytest.
- All of the backend functions written in bash have been converted to python. 
- Zbrains is now importable, and importing src.zbrains will enable the programmatic use of zbrains.
