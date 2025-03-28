import numpy as np 
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
#plt.style.use('dark_background')
fields_total = ["thermal_energy", "kinetic_energy", "magnetic_energy"]#, "time","i"]
fields_ave   = ["density", "temperature", "pressure","velocity_magnitude", "magnetic_field_magnitude", "magnetic_field_xz"]#,"time","i"]

#print(fields_type["density"])
#quit()
path       = "quikplots/data/"#"hello_anal/"  
def gather_data(path):
    regions    = ["bubble", "shell"]
    sep        = "_"
    data_types = ["max","min","averages","std","total"]

    fields_type = {v:fields_ave for v in data_types[:len(data_types)-1]}
    fields_type.update({v:fields_total for v in data_types[len(data_types)-1:]})

    quantities = []
    file_ext   = ".txt"
    nfloats_ave = 6
    nfloats_tot = 4
    timecol = nfloats_ave
    indexcol= timecol + 1
    nint    = 1

    nfloats_type = {v:nfloats_ave for v in data_types[:len(data_types)-1]}
    nfloats_type.update({v:nfloats_tot for v in data_types[len(data_types)-1:]})

    database = {}
    database.update({"time":None})
    database.update({"index":None})
    first = True
    for region in regions:
        database.update({region:{}})    
        
        for data_type in data_types:
            fields = fields_type[data_type]
            nfloats = nfloats_type[data_type]
            database[region].update({data_type:{}})
            filename = path + region + sep + data_type + file_ext
            if first:  
                with open(filename, mode = "r") as f: database["time"] = np.loadtxt(f, dtype=float, usecols=timecol )
                with open(filename, mode = "r") as f: database["index"] = np.loadtxt(f, dtype=int, usecols =indexcol)
                first = False
            with open(filename, mode = "r") as f:
                
                data = np.loadtxt(f, dtype = float, usecols = range(0,nfloats+1))  #' %.8e'*7+" %d"    
                #data = np.loadtxt(f, dtype = int  , usecols = nfloats+1)
                
                for i, field in enumerate(fields[:len(fields)]):
                    #print(field)
                    #print(field, data[:,i])
                    #print(region, data_type, field)
                    database[region][data_type].update({field:data[:,i]})
                    
                    
    time   = database["time"]
    bubble = database["bubble"]
    shell  = database["shell"]
    return bubble, shell, time

bubble, shell, time = gather_data(path)



style_region = {"shell":"solid","bubble":"dashed"}

quantity = "magnetic_field_magnitude"
data_type="averages"
lw= 5
lw_lab=1
B0 = 1


    #magnetic field
if 1:
    fig0=plt.figure(0)
    ax0 = fig0.add_subplot(211)
    ax01 = fig0.add_subplot(212, sharex=ax0)
    bubble_patch = Line2D([0],[0],color='black', linestyle=style_region["bubble"], label='Hot bubble', lw=lw_lab)
    shell_patch = Line2D([0],[0], color='black', linestyle=style_region["shell"], label='Shell', lw=lw_lab)


    ax0.plot(time, bubble[data_type][quantity]/B0, color="midnightblue",label=r"$B$", lw=lw, linestyle=style_region["bubble"])
    ax0.plot(time, shell[data_type ][quantity]/B0, color="midnightblue", lw=lw, linestyle=style_region["shell"])
    print("bubble B, %.3e"%bubble["max"][quantity].max())
    print("shell  B, %.3e"%shell["max"][quantity].max())

    quantity = "magnetic_field_xz"
    ax01.plot(time, bubble[data_type][quantity]/B0, color="seagreen",label=r"$B_{\rm xz}$", lw=lw, linestyle=style_region["bubble"])
    ax01.plot(time, shell[data_type ][quantity]/B0, color="seagreen", lw=lw, linestyle=style_region["shell"])
    print("bubble Bxz %.3e"%bubble["max"][quantity].max())
    print("shell  Bxz %.3e"%shell["max"][quantity].max())

    #handles, labels = plt.gca().get_legend_handles_labels()
    #handles.extend([bubble_patch, shell_patch])
    #ax0.legend(handles=handles)

    #ax.set_ylabel(r"$n\,[\rm cm^{-3}]$")
    ax0.set_ylabel(r"$B\,[\rm \mu G]$")
    ax01.set_ylabel(r"$B\,[\rm \mu G]$")
    ax01.set_xlabel(r"$t\,[\rm kyr]$")
    ax01.set_ylim(1e-6,5e-5)
    #ax0.set_t
    ax0.set_yscale("log")
    
    
if 1:

    fig1=plt.figure(2)
    ax1 = fig1.add_subplot(111)
    bubble_patch = Line2D([0],[0],color='black', linestyle=style_region["bubble"], label='Hot bubble', lw=lw_lab)
    shell_patch = Line2D([0],[0], color='black', linestyle=style_region["shell"], label='Shell', lw=lw_lab)

    quantity = "temperature"
    ax1.plot(time, bubble[data_type][quantity]/B0, color="midnightblue",label=r"$T$", lw=lw, linestyle=style_region["bubble"])
    ax1.plot(time, shell[data_type ][quantity]/B0, color="midnightblue", lw=lw, linestyle=style_region["shell"])
    print("bubble B, %.3e"%bubble["max"][quantity].max())
    print("shell  B, %.3e"%shell["max"][quantity].max())


    handles, labels = plt.gca().get_legend_handles_labels()
    handles.extend([bubble_patch, shell_patch])
    ax1.legend(handles=handles)

    #ax.set_ylabel(r"$n\,[\rm cm^{-3}]$")
    ax1.set_ylabel(r"$T\,[\rm K]$")
    ax1.set_xlabel(r"$t\,[\rm kyr]$")
    ax1.set_yscale("log")
    fig1.savefig("temp_MHB")

    ####  MASS
    sTH0 = 8.95429951e+47
    sK0  = 2.45551595e+47
    sMB0 = 3.98392043e+50
    sV0  = 4.00582554e+60
    sM0  = 1.93324755e+36
    Ei0S = {"thermal_energy":sTH0,"kinetic_energy":sK0,"magnetic_energy":sMB0}

    fig=plt.figure(1)
    ax = fig.add_subplot(111)
    quantity="mass"
    data_type="total"
    bubble_patch = Line2D([0],[0],color="black", linestyle=style_region["bubble"], label='Hot bubble', lw=lw_lab)
    shell_patch = Line2D([0],[0], color="black", linestyle=style_region["shell"], label='Shell', lw=lw_lab)

    colors = {"thermal_energy":"crimson","kinetic_energy":"green","magnetic_energy":"darkblue"}
    labels = {"thermal_energy":"thermal","kinetic_energy":"kinetic","magnetic_energy":"magnetic"}

    for quantity in ["thermal_energy", "kinetic_energy", "magnetic_energy"]:

        ax.plot(time, (bubble[data_type][quantity])/1e51, color=colors[quantity], lw=lw, linestyle=style_region["bubble"],label =labels[quantity])
        #if quantity!="magnetic_energy": 
        ax.plot(time, (shell[data_type ][quantity])/1e51, color=colors[quantity], lw=lw , linestyle=style_region["shell"])

    handles, labels = plt.gca().get_legend_handles_labels()
    handles.extend([bubble_patch, shell_patch])
    ax.legend(handles=handles)
    #ax.set_ylabel(r"$n\,[\rm cm^{-3}]$")
    ax.set_ylabel(r"$E\,[10^{51}\,{\rm erg}]$")
    ax.set_ylim(0,1.1)

    ax.set_xlabel(r"$t\,[\rm kyr]$")
    ax.set_yscale("linear")

    fig.savefig("energies_MHB")



    data_type="averages"
    quantity = "temperature"
    fig3=plt.figure(3)
    ax3 = fig3.add_subplot(111)
    ax3.set_ylabel(r"$T_{\rm TI}/T_{\rm uni}$")
    ax3.set_ylim(0,1.5)

    ax3.set_xlabel(r"$t\,[\rm kyr]$")

  