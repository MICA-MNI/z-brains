import os 
def setenv(wbpath):
    # wbpath = "C:/Users/Ian/Downloads/workbench-windows64-v1.5.0/workbench/bin_windows64"
    os.environ["WORKBENCH_PATH"] = wbpath
    return wbpath