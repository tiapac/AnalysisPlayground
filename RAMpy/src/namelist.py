
from RAMpy.src.logger import logger, oklog
from RAMpy.src.CONFIG import config
from os.path import exists, dirname
from os import makedirs
from typing import Union
from numpy import ndarray
listlike = Union[list, tuple, ndarray]


class block: 
    
    """
        Fortran namelists blocks.
    """
    check_values = True
    
    def __init__(self, 
                 name : str
                 ):
        """
        name: str, name of the namelist's block.
        """
        self.block_name = name
        self.var_list   = []
        pass
    
    @property
    def variables(self):
        """
            list of the name variables
        """
        return self.var_list
    @property
    def values(self):
        """
            list of the name variables
        """
        return [self[name] for name in self.var_list]
    
    @staticmethod
    def value_check(value:str) -> bool:
        #check boolean 
        if "true" in value.lower() or "false".lower() in value:
            
            if value.startswith(".") or value.endswith("."):
                if not value.startswith(".") and value.endswith("."): raise Exception("Boolean value must start and end with a dot '.' !")
    
    def change_var(self,
            name  : str,
            value : str|listlike
            ):
        # checking missing fields 
        _ = self[name]
        self.__setattr__(name, value)
        if oklog: logger.warning(f"Changing in block '{self.block_name:s}' name variable '{name:s}' with value " +  str(value))

    def new_var(self,
            name  : str,
            value : str|listlike
            ):
        """
            add a new variable.
        """
        
        self.__setattr__(name, value)
        
        self.var_list += [name]
        if oklog: logger.trace(f"Adding to block '{self.block_name:s}' name variable '{name:s}' with default value " +  str(value))
        return 
    
    def write_block(self, filename):
        
        #with open(file = filename, mode = "a") as f:
        if isinstance(filename,str):
            f = open(filename, mode="a")
        else:
            f = filename
        f.write(f"&{self.block_name.upper()}\n")
        for name, value in zip(self.variables,self.values):
            if isinstance(value, listlike):
                f.write(f"{name:s} = ")
                for ii, val in enumerate(value): 
                    f.write(f"{val}")
                    if ii < len(value)-1: 
                        f.write(f", ")
                    else: 
                        f.write(f"\n")
                        
            else:  
                f.write(f"{name.strip():s} = {value.strip():s}\n") 
        f.write("/\n\n")
        return 
       
    def __getitem__(self, key):
        try:
            return  self.__getattribute__(key)
        except AttributeError:
            raise Exception(f"Parameter {key} is not in namelist block {self.block_name}")
 
    
    @staticmethod
    def __get_first(val):
        val = val.split("!")[0]
        val = val.split()[0]
        # check if fortran double, very arbitrary
        try: 
            val = float(val.replace("d","e"))
        except Exception as e: 
            pass
        return val 

    def iterable(self): yield {name:self.__get_first(self[name]) for name in self.variables}

    def show(self):
        return self.__str__()
    
    def __iter__(self):
        return  iter({name:self.__get_first(self[name]) for name in self.variables})
    def __str__(self):
        string =  f"{self.block_name}&\n" + "".join([f"{name} = {self[name]} \n" for name in self.variables])
        return string
    def __repr__(self): return self.__str__()

class namelist: 
    """
        Fortran namelist.
    """
    _OVERWRITE = False
    def __init__(self, name):
        self.name = name
        self.blocks_list = []
        pass
    
    def add_block(self, name : str ) -> block:
        
        name = (name.strip("\n")).strip()
        try:
            assert name not in self.blocks_list
        except:
            raise Exception(f"Block {name:s}  of the variable already exist in this namelist.")
        
        self.__setattr__(name, block(name = name))
        self.blocks_list += [self[name]]
        if oklog: logger.trace(f"Adding to namelist '{self.name:s}' block '{name:s}' ")
        return self[name]
    @property
    def blocks(self):
        return [b for b in self.blocks_list]
    
    @property
    def noutput(self):
        tout = (self["OUTPUT_PARAMS"]["tend"].split("!")[0]).strip()
        dt = (self["OUTPUT_PARAMS"]["delta_tout"].split("!")[0]).strip()
        noutput = int(float(tout)/float(dt)) + 1
        return noutput   
    
    def read(self, filename:str):    
        with open(file = filename, mode = "r") as f: 
            for line in f:                 
                if line.startswith("&"): 
                    tmpblock = self.add_block(name = line[1:])                    
                    continue
                
                if "=" in line: 
                    splitted = line.split(sep = "=")
                    varname = splitted[0].strip()
                    varvalue = ("".join(splitted[1:])).strip("\n")
                    
                    tmpblock.new_var(name = varname, value = varvalue  ) 
                
                elif line[0] == "/" or line=="" or line=="  " or line==" ":
                    continue
    def set_overwriting(self, value: bool = True)-> None: 
        self._OVERWRITE = value
        if oklog: logger.trace("Overwriting of namelists allowed!")
        return
    # def show(self):
    #     print(self.blocks_list[1])
    #     quit()#: print(type(block) )

    def write(self, filename):
        if not filename.endswith(config._NML_SUFFIX): filename+=config._NML_SUFFIX
        foldername = dirname(filename)
        if foldername.strip()!="" and not exists(foldername): makedirs(foldername)
        mode = "x" if not self._OVERWRITE else "w"
        with open(file = filename, mode = mode) as f:
            for iblock in self.blocks:
                iblock.write_block(f)
                
    def __getitem__(self, key):
        return  self.__getattribute__(key)

    def __repr__(self):
        #for b in self.blocks_list: print(str(b))
        return "".join([str(blck)+"\n" for blck in self.blocks_list])
    
        
