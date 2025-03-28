
from subprocess import check_output
from datetime import datetime
import glob
from os.path import dirname, basename, abspath, isfile

def _last_snap(path:str, prefix = "output_", verbose = False) -> int:
    folder = "%s/*"%path
    
    names  = glob.glob(folder)
    fnames = [] 
    for name in names: 
        if prefix in name: 
            fnames += [basename(name).removeprefix(prefix)]
    for idx, formatted_number in enumerate(fnames.copy()):
        fnames[idx] = int(formatted_number)
    
    return max(fnames) if len(fnames)>0 else 0 

 
def _last_backup(path:str, prefix = "backup_", verbose = False) -> int:
    folder = "%s/"%path
    searchpath = folder + prefix + "*/info*"
    names = glob.glob(searchpath)
    
    commanf = ['grep', '-ri','time' ]
    commanf.extend(names)
    out = check_output(commanf, text=True)
    lines = out.split("\n")
    lines.pop()
    lines = sorted(lines)

    backtime = [0.0]*len(lines)
    for idx,line in enumerate(lines): 
        backtime[idx] = float(line.split("=")[1])
    most_recent_idx = max(enumerate(backtime), key = lambda nv: nv[1])[0]
    most_recent     = most_recent_idx +1
    
    return  most_recent


def last_snap(*args,backup = False, **kwargs):
    return _last_backup(*args, **kwargs) if backup else _last_snap(*args, **kwargs) 


def last_in_realtime_backup(path:str, prefix = "backup_", verbose = True) -> int:
    folder = "%s/*"%path
    
    out = check_output([
        "ls -lrt  %s/namelist.txt | awk '{print $6, $7, $8, $9}'"%folder
        ], text=True, shell=True)
    
    if verbose : print(out)
    
    lines = out.split("\n")
    lines.pop()
    lines = sorted(lines, key= lambda k: k.split(prefix)[1])
    backtime = [0.0]*len(lines)
    for idx,line in enumerate(lines): 
        #backtime[idx] = float(line.split("=")[1])
        backtime[idx] = datetime.strptime((line.split(prefix)[0]).strip(), "%b %d %H:%M").replace(year=datetime.now().year)
       
    #print(backtime)
    most_recent_idx = max(enumerate(backtime), key=lambda nv: nv[1])[0]
    most_recent     = most_recent_idx +1
    return most_recent

#grep -ri --include="info*" "time" **/        
#-ri --include="info*" "time" **/        

#grep -ri --include="info*" "time" **/        
#-ri --include="info*" "time" **/        