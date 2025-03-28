import RAMpy


sbatcher = RAMpy.sbatcher
namelist = RAMpy.namelist
runner   = RAMpy.runner

DEF_EXEC = RAMpy.config.DEF_EXEC

print(RAMpy.config.TEMPLATES_PATH)

quit()
if __name__ == "__main__":
            
    job_names = ["HighBeta_2e-26","Mbeta_2e-26","lowBeta_2e-26","HighBeta_HM","Mbeta_HM","lowBeta_HM"]
    bfields   = ["2.7378525496449054d-07", "2.7378525496449053d-06","2.7378525496449054d-05","2.7378525496449054d-07", "2.7378525496449053d-06","2.7378525496449054d-05"]
    dirs      = ["HighBeta","Mbeta","lowBeta_HM","HighBeta_HM","Mbeta_HM","lowBeta_HM"]
    photo0    = ["0","0","0","1","1","1"]
    photo1    = ["1","1","1","0","0","0"]
    photo2    = (["/work/pacicco/grackle/grackle_data_files/input/CloudyData_noUVB.h5"]*3)
    photo2.extend(["/work/pacicco/grackle/grackle_data_files/input/CloudyData_UVBHM2012.h5"]*3)

    def quoted(string): return f"'{string}'"

    sbatch = sbatcher("myjob")
    
    sbatch.read("./templates/template.job")    
    sbatch.umodule("ompi")
    sbatch.lmodule("ompi/intel/2.1.1")	

    mynamelist = namelist("readRAMSES")
    mynamelist.read("./templates/template.nml")
    mynamelist["REFINEMENT"].change_var("tmp_max_refine","7")
    #mynamelist["INIT_PARAMS"].change_var("initfile(1)","ic_part0")
    sbatch.set_overwriting(True)    
    mynamelist.set_overwriting(True)
    for jn, dir, b0, ph0, ph1,ph2 in zip(job_names, dirs, bfields, photo0,photo1,photo2):
        nml_name = f"nmls/test_nml_{dir}.nml"
        sbatch_name = f"sbatchs/test_sbatch_{dir}.job"
        
        mynamelist["IOdump"].change_var("out_path", f"'{dir}'")
        mynamelist["IOdump"].change_var("back_path", f"'backup_{dir}'")
        mynamelist["IC_OBJECTS"].change_var("disk%By0", f"{b0}")
        mynamelist["GRACKLE_PARAMS"].change_var("grackle_UVbackground", f"{ph0}")
        mynamelist["GRACKLE_PARAMS"].change_var("grackle_photoelectric_heating", f"{ph1}")
        mynamelist["GRACKLE_PARAMS"].change_var("grackle_data_file", f"{quoted(ph2)}")
        mynamelist.write(nml_name)
        sbatch.set_runner(runner(200, DEF_EXEC, nml_name))
        sbatch["job-name"].change(dir)
        sbatch.write(sbatch_name)
        #sbatch.run(sbatch_name)

    