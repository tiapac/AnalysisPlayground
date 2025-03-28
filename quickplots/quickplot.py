from utils.myparser import myparser


def MakePlot(parser = myparser(prog="Quickplot", overwrite=True), subc_args = None ):



	import matplotlib.pyplot as plt
	from multiprocessing import Pool
	import os
	import sys
	import numpy as np
	import utils.utils as utils
	import paper_settings as pg
	import localYt.yt.yt as yt
	#yt.enable_parallelism()
	#yt.frontends.ramses.api.OUTPUT_DIR_EXP="hello"
	#from yt.data_objects.level_sets.api import *
	from matplotlib.colors import LogNorm
	#yt.funcs.matplotlib_style_context("dark_background")
	import time as TIME
	#initfield(MHD=False)####NEEDED TO INITIALIZE NEW FIELS
	from analysis_playground.constants import factG_in_cgs

	from matplotlib.ticker import ScalarFormatter, LogFormatter

	def without_keys(d, keys):
		return {x: d[x] for x in d if x not in keys}

	paper = pg.PaperSettings()
	def add_pscalar(i):
		def _scalar(field, data):
			return (
				data["gas", "density"]
				* data["ramses", "hydro_scalar_%02d"%i]
			)
		yt.add_field(
		name=("gas", "dscalar_%02d"%i),
		function=_scalar,
		sampling_type="local",
		#units="dyne/cm**2",
		)
	def _jeansLength(field, data):
		G = yt.YTArray(6.67430e-8, "cm**3/g/s**2") #G = factG_in_cgs*(yt.units.dyne*yt.units.cm/yt.units.s**2)
		cs =  data["gas", "sound_speed"].in_cgs()
		rho=  data["gas", "density"].in_cgs()
		return np.sqrt(np.pi * cs**2 / (G * rho))
		
	yt.add_field(
		name=("gas", "jeans_length"),
		function=_jeansLength,
		sampling_type="local",
		#units="cm",
		)

	
	#parser = argparse.ArgumentParser(prog="Quickplot")
	parser.add_default()
	parser.add_argument("n1",	help="enter at least one output number", type=int)
	parser.add_argument("-n2", 	help = "Make plots up to output number n2. If None, sets n2 = n1.", type=int, default = None)
	parser.add_argument("quantity",	help ="Quantity to analyse.", default = "density"   , type = str)
	parser.add_argument("-nl","--nolog",  help ="plot linear variable."  , action  = "store_true", default=False)
	parser.add_argument("-b","--box",	  help ="Box size to be plotted.", type = int)
	parser.add_argument("-ax","--axis",	  help ="Set the axis to slice/project on.", default="x",type = str)
	parser.add_argument("-p","--particles",help ="Add particles to the plot.", action  = "store_true", default=False)
	parser.add_argument("-sb","--stream_B",help ="Add magnetic field streamlines.", action  = "store_true", default=False)
	parser.add_argument("-sv","--stream_V",help ="Add velocity field streamlines.", action  = "store_true", default=False)
	parser.add_argument("--conv_B",help ="Add magnetic field line integral convolution.", action  = "store_true", default=False)
	parser.add_argument("--conv_V",help ="Add velocity field line integral convolution.", action  = "store_true", default=False)

	parser.add_argument("-pj","--projection",help ="Do a projection instead of a slice.", action  = "store_true", default=False)
	parser.add_argument("--weight",help ="Weight to be used in the projection, if weighted.", default=None,type = str)
	parser.add_argument("--path",help ="path to the directory containing outputs.", default="./outputs",type = str)
	parser.add_argument("--units",help ="Units of the colormap.", default=None,type = str)

	parser.add_argument("--grid",help ="Plot the grid.", action  = "store_true", default=False)
	parser.add_argument("-c","--center",help ="Center of the plot, three cooridinates between 0.0 <= x_i <= 1.0] are needed ", default =(0.5,0.5,0.5), nargs="+")
	parser.add_argument("-d","--depth", help ="Depth of the projection along the --axis ", default = None,type=float)
	parser.add_argument("--vlim",help ="Limits on the colormap.", default = None, type=float,      nargs=2)
	parser.add_argument("--cmap",help ="Colormap to use.", default=None,type = str)
	parser.add_argument("--fields",help ="Show field and exit.", default=False,action  = "store_true")
	parser.add_argument("-dir", help = "Name of the director where to store the profiles.", default = None)
	parser.add_argument("-tuo","--time_units_offset", help = "Time unit and offset in the same units.", default =["Myr",0.0] , nargs=2)
	parser.add_argument("-xb","--xray_band", help = "minimum and maximum ernergy for the xray emission.", default =[0.5,3.0] , type= float, nargs=2)
	parser.add_argument("--xray_source", help = "Model to compute xray spectrums.", default ="cloudy")#"spex" )
	parser.add_argument("--arepo", help = "Analyse arepo files.", action="store_true", default=False )
	parser.add_argument("-s","--step", 	help="step in the data. Use only if n2 is used.", type=int, default=1)
	parser.add_argument("-cpus", 	help="step in the data. Use only if n2 is used.", type=int, default=None)

	parser.add_argument("-bcode","--box_code_units", 	help="use code units instade of image coordinates.", action  = "store_true", default=False)
	parser.add_argument("-clumps", 	help="Run the clump finder.", action  = "store_true", default=False)
	parser.add_argument("-clump_setting", 	help="Run the clump finder.", nargs=1)
	parser.add_argument("-info", 	help="Print yt info.", action  = "store_true", default = False)
	parser.add_argument("-noscale", 	help="Annotate scale.", action  = "store_true", default=False)
	parser.add_argument("-notime", 	help="Annotate scale.", action  = "store_true", default=False)
	parser.add_argument("-noaxlabels", 	help="Annotate scale.", action  = "store_true", default=False)

	parser.add_argument("-arrow", 	help="Annotate arrow.", default=None, choices=["x","y","z"])
	parser.add_argument("-filter", 	help="Filter with a gaussian beam.", default=None, nargs=2, type=int)
	parser.add_argument("-mpl", 	help="Second layout.", action  = "store_true", default=False)

	parser.add_argument("-contour", 	help="Second layout.", type  = str, default=None)
	parser.add_argument("-total", 	help="Second layout.", action  = "store_true", default = False)


	#slc.annotate_arrow((0.5, 0.5, 0.5), length=0.06, color="blue")

	colormaps = paper.cmaps.colormaps
	textcolor = paper.cmaps.textcolor
	args, _ = parser.parse_known_args(subc_args,fix_n1n2 = True, fixdir=True)
	print(args.path, args.n1, args.n2)
	if args.dir is None: args.dir = args.path+"/quikplots/maps"


	if args.cmap is None: args.cmap = colormaps.get(args.quantity, "viridis")
	print("using colormap %s"%colormaps.get(args.quantity, "viridis"))

	myunits = paper.units if not args.projection else paper.unitsPj

	if args.cpus is not None:
		from multiprocessing import Pool
		from itertools import repeat
	if not args.info: yt.set_log_level(50)### DO NOT SHOW THE INFO IN THE STDOUT

	if not args.arepo:
		try:
			from mylibraries.addedfields import initfield
			initfield_exitsts = True
			initfield(MHD=True)
		except Exception as e :
			print(e, "File 'mylibraries.addedfields' not found.") 
			initfield_exitsts = False
	args.time_units_offset[1]=float(args.time_units_offset[1])
	args.center = np.array(args.center,dtype=float)
	args.depth	= (args.depth,"pc")
	# if args.n2 is None: args.n2=int(args.n1)
	# if  args.n1 > args.n2: 
	# 	n2tmp	= int(args.n2)
	# 	args.n2 = int(args.n1)
	# 	args.n1 = n2tmp
	quant_name  = args.quantity
	weight_name = args.weight

	if "dscalar"in quant_name:  
		for i in range(0, 2 + 1): 
			add_pscalar(i)
	
	if "xray" in args.quantity: 
		import pyxsim
		source_model = pyxsim.CIESourceModel(args.xray_source,
											emin   	 = 0.05,
											emax 	 = 11. , 
											#kT_min = 0.025, kT_max = 64,
											nbins  	 = 1000, 
											#max_density=max_density,
											Zmet     = 1.0 ,
											binscale = "log" )
		bandmin = args.xray_band[0]
		bandmax = args.xray_band[1]
	snaplist=range(args.n1, args.n2, args.step)                   
			
	def make_a_plot(i,args=args):
		
		quant_name  = args.quantity
		weight_name = args.weight
		print("Plotting %s of snapshot %i"%(quant_name, i))
		quant= ('gas', quant_name) if "hydro" not in quant_name  else ("ramses", quant_name)
		if args.weight is not None: 
			weight_field = ('gas', weight_name) if "hydro" not in weight_name  else ("ramses",weight_name)
		else: 
			weight_field = None
		if not args.projection:
			plottype = "slices"  
		else:
			if args.weight is None:
				plottype = "projections"
			else:
				plottype = "weighted_projections/%s"%args.weight
		if args.path is None: 
			outpath = os.getcwd() + "/"
		else:
			outpath = args.path+"/"
		outname     = outpath+"output_%05d"%i if not args.arepo else "snap_%03d.hdf5"%i
		#yt.frontends.ramses.definitions.OUTPUT_DIR_EXP='backup_(\\d{5})'
		ds = yt.load(outname, default_species_fields='ionized')
		t  = ds.current_time.to("Myr")#.value
		if "xray" in args.quantity: 
			xray_fields  = source_model.make_source_fields(   ds, bandmin, bandmax  )	
			print("Generated the X-ray fields:")
			print(xray_fields)
			quant = ("gas", args.quantity + "_%s_%s_keV"%( str(bandmin), str(bandmax) ))
			if args.total:
				ad = ds.all_data()
				print(ad.quantities.total_quantity([quant]))
		if args.fields: 
			for f in ds.derived_field_list:
				print(f)
			quit()
		
		
		if args.box is None:
			width = float(ds.domain_width[0].to("pc").value)
		else: 
			width = args.box
		boxlen 	  = float(ds.domain_width[0].to("pc").value)

		if args.filter is not None:
			try: 
				os.makedirs(args.dir+"frb"+"/"+plottype+"/"+quant_name+"/"+args.axis)
			except:
				pass
			slc = ds.slice(args.axis, 0.5) if not args.projection else ds.proj(quant, args.axis)
			frb = slc.to_frb((width,"pc"), 1024)
			frb.apply_gauss_beam(nbeam = args.filter[0], sigma =args.filter[1])
			fig= plt.figure(1)
			ax = fig.add_subplot()
			ax.imshow(frb[quant].d, cmap = args.cmap, norm=None if args.nolog else LogNorm(), )
			ax.set_xticks([])
			ax.set_yticks([])
			fig.tight_layout()
			fig.savefig(args.dir+"frb"+"/"+plottype+"/"+quant_name+"/"+args.axis+"/%05d" %i, dpi=1024,transparent=True, pad_inches=0.0)

		else:

			center = (args.center).copy()
			if args.box_code_units: center = abs(np.array(args.center)+boxlen*0.5)/boxlen

			if not args.projection:
				s  = yt.SlicePlot(ds, args.axis, quant, width=(width,"pc"), center=center) #, center=center) #, no_ghost=False)### do the slice plto
				if ("velocity" in quant_name or "speed" in quant_name) and not "div" in quant_name and not  args.projection :s.set_unit(quant, "km/s")		
			else:
				s = yt.ProjectionPlot(ds, args.axis, quant, weight_field = weight_field, width=(width,"pc"))	
			if args.units is not None: s.set_unit(quant, args.units)
			s.set_cmap(field = "all", cmap = args.cmap)
			if  args.nolog: s.set_log(field = quant, log = False)
			if args.grid: s.annotate_grids()
			
			axses = utils.find_axes(args.axis)
	
			if args.axis == "y":
					s.swap_axes()
					axses.reverse()
		
			if args.conv_B:
				s.annotate_line_integral_convolution(field_x = ("gas", "magnetic_field_" + axses[0]),
													field_y = ("gas", "magnetic_field_" + axses[1]),
													texture = None, alpha = 0.5)

			if args.conv_V:
				s.annotate_line_integral_convolution(("gas","velocity_"+axses[0]),("gas","velocity_"+axses[1]))

			if args.stream_B:	
				s.annotate_streamlines( ("gas","magnetic_field_"+axses[0]),
										("gas","magnetic_field_"+axses[1]),
										# density	  = 5,
										# color	  = "midnightblue",
										# arrowsize = 0.2,
										# linewidth = 0.4)
										broken_streamlines =True,
										density	  = 2,
										color	  = "black",
										arrowsize = 0.0,
										linewidth = 1)
			if args.stream_V:
				s.annotate_streamlines(("gas","velocity_"+axses[0]),
									("gas","velocity_"+axses[1]),
										density	  = 5	   	 	   ,
										color	  = "crimson"	   ,
										arrowsize = 0.5  	 	   ,
										#broken_streamlines=False
										)
			if args.particles		: s.annotate_particles(width = width, marker = "*", p_size = 100,)
			if args.vlim is not None: s.set_zlim(quant, args.vlim[0],args.vlim[1])
			if args.clumps			: print("Not implemented in ramses.")
					
			#set the buffer size so that the dpi is 512
			dpi = paper.figure.dpi
			figuresize = (paper.onecol,paper.onecol) 
			buff_size  =  np.rint(np.ceil(max(figuresize) * dpi)) 
			
			s.set_figure_size((paper.onecol,paper.onecol))		
			s.set_buff_size(buff_size)
			
			
			if args.mpl:
				plt.rcParams['axes.facecolor'] = 'black'
				plt.rcParams['figure.facecolor'] = 'black'

				plt.rcParams.update({"text.usetex": True, **paper.mplfigfontstyle})
				plt.rcParams['axes.facecolor'] = 'black'


				s.hide_colorbar()
					# hide the axes, while still keeping the background color correct:
				s.hide_axes()
	
				s._setup_plots()
				fp     = s._font_properties
				plot   = s.plots[args.quantity if "xray" not in args.quantity else args.quantity + "_%s_%s_keV"%( str(bandmin), str(bandmax) )]
				fig    = plot.figure
				ax  = plot.axes
				img = ax.images[0]
				img.set_interpolation("bicubic")
				
				new_ax = fig.add_axes((0.1,0.9,0.8,0.035))
				ax.set_facecolor("black")
				fig.canvas.draw()
				fig.canvas.flush_events()
				#ax.draw(renderer=True)
				# img.cmap.set_under('black')
				cb     = plot.figure.colorbar(plot.image,
										new_ax,
										extend	    = 'both',
										orientation = 'horizontal'
											)
				
				cb.set_label(myunits.symbols[args.quantity], color=textcolor[args.quantity], **paper.figfontstyle)
				cb.ax.xaxis.set_label_position('top')
				cbxmin, cbxmax = cb.ax.get_xlim()
				
				# if (cbxmax/cbxmin < 10.) or force_calar:
				# # 	cb.ax.xaxis.set_major_formatter(ScalarFormatter())
				# # 	cb.ax.xaxis.set_minor_formatter(ScalarFormatter())
				# 	cb.ax.xaxis.set_major_formatter(LogFormatter())
				# 	cb.ax.xaxis.set_minor_formatter(LogFormatter())
					
				ax.tick_params(reset=True,
					axis       = "both" ,
					which      = "both",     # both major and minor ticks are affected
					bottom     = False ,      # ticks along the bottom edge are off
					top	       = False ,      # ticks along the top edge are off
					labelbottom= False ,      #
					left       = False ,      #
					right      = False )		# labels along the bottom edge are off
				cb.ax.xaxis.set_tick_params(colors = textcolor[args.quantity], which = "both")
				cb.ax.yaxis.set_tick_params(colors = textcolor[args.quantity], which = "both")
			
				plt.setp(plt.getp(cb.ax.axes, 'yticklabels'), color = textcolor[args.quantity])
				plt.setp(plt.getp(cb.ax.axes, 'xticklabels'), color = textcolor[args.quantity])
	
				# set colorbar edgecolor 
				cb.outline.set_edgecolor(textcolor[args.quantity])

				axes = ["x","y","z"]
				axes.remove(args.axis)
				axesl1, axsesl2 = axes
				if args.axis=="y": axesl1, axsesl2 = "z","x"
				
				# the x 
				if not args.noaxlabels:
					s.annotate_text((0.5, 0.05),axesl1, coord_system="axis", text_args=paper.figfontstyle)
					# the y

					s.annotate_text((0.05, 0.5),axsesl2, coord_system="axis", text_args=paper.figfontstyle)
					
			paper.figfontstyle.update({"color" : textcolor[args.quantity]})
			if not args.notime: s.annotate_timestamp(time_offset=(args.time_units_offset[1]  		   ,
											args.time_units_offset[0]) 		   ,
											time_unit = args.time_units_offset[0],
											text_args = paper.figfontstyle	   )
			if not args.noscale: s.annotate_scale(coeff = 0.1 * boxlen, size_bar_args={"color":textcolor[args.quantity]}, text_args    = without_keys(paper.figfontstyle, ["color"]))
			if args.contour is not None: 
				
				s.annotate_contour(("gas",args.contour),
						plot_args={"linewidths":0.2,
									#"colors":textcolor[args.quantity],
									#"linestyles":"dashed",
									"cmap":"jet"
								},
						#text_args={"colors":textcolor[args.quantity],
					#			"fontsize":4,
				#				"inline":1},
						#take_log=True,#args.nolog,
						label  = False,
						#levels=3,
						levels = 7,
						
						)
			#print(dpi)
			print("Saving %s"%args.dir+"/"+plottype+"/"+quant_name+"/"+args.axis+"/%05d" %i)
			s.save(args.dir+"/"+plottype+"/"+quant_name+"/"+args.axis+"/%05d" %i, mpl_kwargs=dict(dpi=dpi))
			#s.save(args.dir+"/"+plottype+"/"+quant_name+"/"+args.axis+"/%05d" %i)
			
	# quickplot 91 xray_emissivity -pj --path=./  -contour magnetic_field_magnitude --vlim 3e-6 4e-5 -xb 0.1 9.0 -mpl
	
	
	
	
	
			if args.filter is not None:
				try:
					os.makedirs(args.dir+"frb"+"/"+plottype+"/"+quant_name+"/"+args.axis)
				except:
					pass
				if 0:
					frb = s.frb
					img = frb.get_image((quant))  
					figure = plt.imshow(img.T,interpolation="gaussian", cmap="hot")
				else:
			
					slc = ds.slice("z", 0.) if not args.projection else ds.proj("z", 0.)
					frb = slc.to_frb((width,"pc"), 512)
					frb.apply_gauss_beam(nbeam=30, sigma=2.0)
					plt.imshow(frb["gas", "density"].d)
					plt.savefig("frb_filters.png")
				
				figure.savefig(args.dir+"frb"+"/"+plottype+"/"+quant_name+"/"+args.axis+"/%05d" %i)
		print(i, "Time=%.4f"%t)#,"Tmax=%.4f"%temp.min())

	if args.cpus is not None:
		with Pool(processes=args.cpus)as p:
			results=p.starmap_async(make_a_plot, zip(snaplist, repeat(args)) )
			results.get()###NEEDED TO OBTAIN THE RESULTS IN ASYNC
			p.close()
	else:
		for i in snaplist:
			make_a_plot(i,args)



	#s.annotate_cell_edges()
	print("Done.")
if __name__=="__main__":
	MakePlot()