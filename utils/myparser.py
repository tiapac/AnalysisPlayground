import os
import argparse
from .logger import setup_logger
from .utils import checkstack
from utils.ErrorHandling import Errors

class myparser(argparse.ArgumentParser):
    def __init__(self, overwrite = True,**kwargs, ):
    

        conflict_handler = "resolve" if overwrite else "error"
        super().__init__(**kwargs, conflict_handler = conflict_handler)
        self.logger = setup_logger(f"MyParser.{self.prog}")
        self.logger.setLevel(self.logger.loglevels.WARNING)

        #if overwrite: self.logger.warning("Parser arguments can be overwritten. Careful.")
    def add_default(self):
        self.add_argument("-n1", 	help="enter at least one output number", type=int)
        self.add_argument("-n2"  , help = "Make plots up to output number n2. If None, sets n2 = n1.", type=int, default = None)
        self.add_argument("-nl","--nolog",  help ="plot linear variable."  , action  = "store_true", default=False)
        self.add_argument("-ax","--axis",	  help ="Set the axis to slice/project on.", default="x",type = str)
        self.add_argument("--path"  ,help ="path to the directory containing outputs.", default="./outputs",type = str)
        self.add_argument("--fields",help ="Show field and exit.", default=False,action  = "store_true")
        self.add_argument("-dir", help = "Name of the director where to store the profiles.", default = None)
        self.add_argument("-cpus", 	help="step in the data. Use only if n2 is used.", type=int, default=None)
        pass
    def parse_args(self,*args, fixdir = False, fix_n1n2 = False, **kwargs):
        self.args = super().parse_args(*args, **kwargs)
        
        if fixdir:
            if self.args.dir is None : self.args.dir = self.args.path+"./" 
        if fix_n1n2: self.args.n1, self.args.n2 = self._fix_n1n2(self.args)
        return self.args
    
    def parse_known_args(self,*args,fixdir = False, fix_n1n2 = False,  **kwargs):

        self.args, subres = super().parse_known_args(*args, **kwargs)
        if fixdir: 
            if self.args.dir is None: self.args.dir = self.args.path+"./" 
        if fix_n1n2: self.args.n1, self.args.n2 = self._fix_n1n2(self.args)
        return self.args, subres
    
    @staticmethod
    def _fix_n1n2(parser):
        n1 = parser.n1
        n2 = parser.n2
        if n2 is None: n2=int(n1)+1
        if  n1 > n2: 
            n2tmp	= int(n2)
            n2 = int(n1)
            n1 = n2tmp
        return [n1,n2]

    def set_outputdir(self, path):
        self.args.dir = path
        if not os.path.exists(self.args.dir): os.makedirs(self.args.dir)
        return
