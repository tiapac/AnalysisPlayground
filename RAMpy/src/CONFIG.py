from os.path import dirname, basename, abspath, isfile
import glob

import RAMpy.templates as tmpls
class config: 
    TEMPLATES_PATH = dirname(tmpls.__file__)

    _DEFAULT_PREFIX = "#SBATCH"
    _MODULE_CMD = "module"
        
    PARALLEL_RUN = "mpirun"
    DEF_EXEC = "/fast/pacicco/ramses_inkTurb/ramses_ink/bin/ramsesMHD3d "
    DEF_NML  = "./sb2.nml"
    NTASKS_FLAG = "-np"

    _JOB_SUFFIX = ".job"
    _NML_SUFFIX = ".nml"
    
    def __init__(self):
        pass
    
    def get_templates(self):
        files = glob.glob(self.TEMPLATES_PATH+"/*.*")
        return self.TEMPLATES_PATH, [basename(f) for f in files if (isfile(f) and "__.py" not in f)]