
import RAMpy
from utils.myparser import myparser 
from utils.logger import setup_logger
sbatcher = RAMpy.sbatcher
namelist = RAMpy.namelist
runner   = RAMpy.runner

DEF_EXEC = RAMpy.config.DEF_EXEC

CONFIGS = RAMpy.config()

def quoted(string): return f"'{string}'"

if __name__ == "__main__":
    logger = setup_logger("Run")
    def print(*args, **kwargs): logger.info(*args, **kwargs)
    #plt.rcParams.update({'font.size': 50})
    parser = myparser(prog="Pynamelist")
    parser.add_argument("-d","--dry", action="store_true", default=False)
    parser.add_argument("-ow","--overwrite",	help="Write files.",action="store_true", default=True)
    parser.add_argument("-wn","--write_namelist",	help="Write files.",action="store_true", default=True)
    parser.add_argument("-s","--sbatch",	help="Sbatch the sbatchers.",action="store_true", default=False)
    parser.add_argument("-path", help="Path to namelist etc.", default="./")
    parser.add_argument("outname", help="Path to namelist etc.")
    parser.add_argument("-names", help="job name.",nargs="+", default=["./"])

    parser.add_argument("-r","--restart",	help="Sbatch the sbatchers.",action="store_true", default=False)
    parser.add_argument("-rfb","--restart_from_backup",	help="Sbatch the sbatchers.",action="store_true", default=False)

    args = parser.parse_args()
    if args.restart_from_backup: args.restart = True
    logger = setup_logger("RamsesRunner")
    loglevels = logger.loglevels
    logger.setLevel(loglevels.INFO)

    def main(args):

        path = args.path 
        # templates_path, templates = CONFIGS.get_templates()      
        job_names = args.names
        sbatch = sbatcher("myjob")
        
        sbatch.read("launch.job")    
        sbatch.umodule("ompi")
        sbatch.lmodule("ompi/intel/2.1.1")	

        nml_name = "sb2.nml"
        mynamelist = namelist("readRAMSES")
        #mynamelist.read(f"{nml_name}.nml")
        mynamelist.read("sb2.nml")
        outdirname = args.outname 
        backup_dirname = f"backup_{outdirname}/" 
        searchpath = outdirname if not args.restart_from_backup else backup_dirname
        if args.overwrite:
            sbatch.set_overwriting(args.overwrite and not args.dry)    
            mynamelist.set_overwriting(args.overwrite and not args.dry)
        for i, jn in enumerate(job_names):
        
            tmp_nml_name = f"{path}{nml_name}"
            tmp_sbatch_name = f"{path}/launch.job"
            last_snapshot = RAMpy.last_in_realtime_backup(searchpath, verbose = True)
            mynamelist["GRACKLE_PARAMS"].change_var("grackle_data_file","'/work/pacicco/grackle/grackle_data_files/input/CloudyData_UVB=HM2012.h5'")
            #last_snapshot = RAMpy.last_snap(searchpath, backup=args.restart, verbose = True)#, backup=True)
            if args.restart: 
                restart = last_snapshot
            else: 
                restart = 0 

            print(f"Starting from restart {11}")
            

            mynamelist["RUN_PARAMS"].change_var("nrestart", "%i"%11)

            myrunner =  runner( sbatch.ntasks.value,
                                DEF_EXEC,
                                nml_name,
                                restart = restart)
            
            sbatch.set_runner(myrunner)

            sbatch["job-name"].change(jn)

            
            # print(f"Sbatch job is:\n {sbatch}")
            
            # print(f"Namelist is: \n {mynamelist}")
            
            # print(f"running with:\n {myrunner}")
            
            if args.overwrite: 
                logger.warning("Overwriting set to True.")
                if input("Write yes to continue: ")!= "yes": quit("Exiting.")
            if args.write_namelist:
                if not args.dry:
                        mynamelist.write(nml_name)
            if args.sbatch:
                if not args.dry:
                    
                    sbatch.write(tmp_sbatch_name)
                    sbatch.run(tmp_sbatch_name)

    main(args)