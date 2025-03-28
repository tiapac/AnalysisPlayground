
from RAMpy.src.logger import logger, oklog
from RAMpy.src.CONFIG import config
from os.path import exists, dirname
from os import makedirs
import subprocess#, os

from typing import Union
from numpy import ndarray
listlike = Union[list, tuple, ndarray]






class sbatcher: 
    _OVERWRITE = False
    
    ###############################à
    class sbatchline: 
        
        def __init__(self,line, prefix = None, var = None, value = None):
            if prefix is None and var is None and value is None:
                prefix, var, value = self.getparts(line)
            self.prefix = prefix
            dashes = ((prefix.lower()).strip("sbatch")).strip()
            self._COMMENTED = False
            if len(dashes)>1: self.set_active(False)
            self.var    = var   
            self.value  = value 
            pass
        
        @staticmethod
        def getparts(line):
            splitted = line.split("--")
            prefix = splitted[0].strip()
            cv = ("".join(splitted[1:])).strip("\n").split("=")
            var = cv[0]
            value = None if len(cv)<2 else cv[1]
            return prefix, var, value
        
        def set_active(self, value : bool):
            self._COMMENTED = not value 
            return 
        
        def show(self):
            return f"{self.var:s} is {self.value}\n"
           
        @property
        def active(self):
            return not self._COMMENTED
        def change(self, value):
            self.value = value
            return
        def __str__(self):
            return f"{self.prefix:s} --{self.var:s}={self.value}\n"
           
        def __repr__(self):
            return f"prefix: {self.prefix:s}, var: {self.var:s}, value: {self.value}"

    #####################################
    
           
    def __init__(self, filename ):
        self.filename = filename
        self.shebang = None  
        self.sbatchlines = []       # mm
        self.commands = []          # pre-run 
        self._run_command = None    # to run the simulation 
        self.post_run_commands = [] # run after simulations
        self.modules = []
        pass
    
    def set_overwriting(self, value: bool = True)-> None: 
        self._OVERWRITE = value
        if oklog: logger.warning("Overwriting of SBATCH file allowed!")
        return
    def _add_sbatchline(self, line = None, **kwargs):
        
        if line is None and len(kwargs)<3:
            raise Exception("Not enough arguments. Either supply a line to parse, or prefix, var and value.") 
        if line is not None:
            if not line.strip(): return
        
        tmp = self.sbatchline(line = line, **kwargs)
        if config._DEFAULT_PREFIX not in tmp.prefix: return
        self.sbatchlines += [tmp]
        self.__setattr__(tmp.var,tmp )
        return tmp
    
    def read(self, filename:str):
        with open(file = filename, mode = "r") as f: 
            self.shebang = self.sbatchline(f.readline())
            for line in f:
                self._add_sbatchline(line)
                
    def new_var(self, var : str, value : str, prefix = config. _DEFAULT_PREFIX):
        self._add_sbatchline(prefix = prefix, var = var, value=value)
        return
    
    def new_command(self, command ):
        self.commands += [command]
        return
    def new_post_command(self, command ):
        self.post_run_commands += [command]
        return
    def set_runner(self, command):
        self._run_command = command
        return
    
    @property
    def run_command(self):
        return self._run_command
        
    def write(self, filename:str):
        if not filename.endswith(config._JOB_SUFFIX): filename+=config._JOB_SUFFIX
        foldername = dirname(filename)
        
        if foldername.strip()!="" and not exists(foldername): makedirs(foldername)
        
        if self.run_command is None: raise Exception("Runner was not set. Set a runner befor writing to file.")
        mode = "x" if not self._OVERWRITE else "w"
        with open(file = filename, mode = mode) as f: 
            var = self.shebang
            f.write(f"{var.prefix} --{var.var}\n")
            for var in self.sbatchlines:
                f.write(str(var))#f"{var.prefix} {var.var} {var.value}\n") 
            for module_command in self.modules:
                f.write(f"{config._MODULE_CMD} {module_command}\n")  
            for command in self.commands:
                f.write(f"{command}\n")
                
            f.write(f"{self.run_command}")
            
            for command in self.post_run_commands:
                f.write(f"{command}\n")
            
            
            return
        
    def __getitem__(self, key):
        return  self.__getattribute__(key)
    
    def __repr__(self):
        msg = ""
        for line in self.sbatchlines:
            msg += f"{line}"
        return msg
    
    def run(self, filename):
        return subprocess.run(["sbatch", f"{filename}"], check=True)

    def umodule(self, names):
        if not isinstance(names,listlike): names=[names]
        for name in names: self.modules += [f"unload {name}"]
        return
    
    def lmodule(self, names):
        if not isinstance(names,listlike): names=[names]
        for name in names: self.modules += [f"load {name}"]
        return
    
    def pmodules(self):
        self.modules += [f"purge"]
        return
    
    
    
    
    
    
    
    
    