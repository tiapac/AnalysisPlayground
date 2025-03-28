from .CONFIG import config


class runner:
    
    def __init__(self, ntasks, EXEC , nml, prunner:str = config.PARALLEL_RUN, mpi_args = {}, restart = 0):
        self.RUNNER   = prunner
        self.NTASKS   = {config.NTASKS_FLAG:str(ntasks)}
        self.EXEC     = EXEC
        self.NAMELIST = nml
        self.MPI_ARGS = mpi_args
        self.RESTART = restart
        pass
    
    def __str__(self):
        return self.__repr__()
    def __repr__(self):
        cmd =  f"{self.RUNNER} {config.NTASKS_FLAG} {self.NTASKS[config.NTASKS_FLAG]} "
        for k in (self.MPI_ARGS).keys():
            cmd += f" -{k} {self.MPI_ARGS[k]} "
        
        cmd += f"{self.EXEC} {self.NAMELIST}  {self.RESTART}"        
        return cmd
    
    
