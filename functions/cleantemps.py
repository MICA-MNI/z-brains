import fnmatch
import os
import shutil
def delete_temp_folders(directory):
    for root, dirs, files in os.walk(directory, topdown=False):
        for name in dirs:
            if fnmatch.fnmatch(name, 'z_brains_temp.????????'):
                print(f"Deleting temp folder {name}")
                # shutil.rmtree(os.path.join(root, name))