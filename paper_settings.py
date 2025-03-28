from collections import defaultdict
from unities import cm, g, s, gauss, dyne, kelvin, erg, Pmes, adimensional

import unities as uu



aliases= {
            "dens":"density",
            "ndens":"ndensity",
            "temp":"temperature",
            "press":"pressure",
            "pmag":"magnetic_pressure",
            "|v|":"velocity_magnitude",
            "B":"magnetic_field_magnitude",
            "plasmaB":"plasma_beta"}


points_per_inch = 72.27

# Convert points to inches
def pt_to_in(pt): return pt / points_per_inch
class cmaps: 
    colormaps = {
            "density":"cividis",
            "ndensity":"cividis",
			"temperature":"hot",
			"magnetic_field_magnitude":"jet",
			"velocity_magnitude":"jet",
			"pressure":"plasma",
            "plasma_beta":"RdBu_r",
            "plasma_beta_inv":"RdBu_r",
            "magnetic_pressure":"plasma",
            "xray_emissivity":"hot",
            "hydro_scalar_00":"Reds",   
            }
    #colormaps.setdefault('gg', "")
    # colormaps = defaultdict(lambda:"jet",colormaps)
    textcolor = {"density":"white",
                "ndensity":"white", 
                "temperature":"white",
                "magnetic_field_magnitude":"black",
                "velocity_magnitude":"black",
                "pressure":"white",
                "magnetic_pressure":"white",
                "plasma_beta":"black",
                "plasma_beta_inv":"black",
                "xray_emissivity":"white",
                "hydro_scalar_00":"black",   

                }
    # textcolor = defaultdict(lambda:"black",textcolor)

    def __init__(self):
        
        for key in aliases.keys():
            self.colormaps.update({key:self.colormaps[aliases[key]]})
            self.textcolor.update({key:self.textcolor[aliases[key]]})

        pass

class mysym:
    def __init__(self, name, symbol, units):
        self.name  = name
        self.symbol= symbol
        self.units =  units
        pass


mysymbols = [
    mysym("density", r"$\rho$", g * cm**-3),
    mysym("ndensity", r"$n$", cm**-3),
    mysym("temperature", r"$T$", kelvin),
    mysym("magnetic_field_magnitude", r"$B$", gauss),
    mysym("pressure", r"$P_\mathrm{th}$", Pmes),
    mysym("magnetic_pressure", r"$P_\mathrm{B}$", Pmes),
    mysym("mass", r"$M$", g),
    mysym("velocity_magnitude", r"$v$", cm * s**-1),
    mysym("plasma_beta", r"$\beta$", adimensional),
    mysym("plasma_beta_inv", r"$\beta^{-1}$", adimensional),
    mysym("xray_emissivity", r"$\epsilon_X$", erg * cm**-3 * s**-1),
    mysym("xray_luminosity", r"$L_X$", erg * s**-1),
    mysym("radius", r"$r$", uu.symbols("pc")),
    mysym("theta", r"$\theta$", adimensional),
    mysym("theta2", r"$\theta$", adimensional),
    mysym("sintheta2", r"$\sin{\theta}$", adimensional),
    mysym("hydro_scalar_00", r"$S_0$", adimensional),
]
class Units:
    def __init__(self):
        
        for key in aliases.keys():
            self.symbols.update({key:self.symbols[aliases[key]]})
            self.units.update({key:self.units[aliases[key]]})

        pass

    symbols = {sym.name: sym.symbol for sym in mysymbols}
    
    # symbols = defaultdict(lambda: r"Sh", symbols)
   
    units = {sym.name: sym.units for sym in mysymbols}
    # units = defaultdict(lambda:r"Sh",units)

class UnitsP(Units):
    def __init__(self):
        super().__init__()
        self.units = super().units.copy()
        for key in self.units.keys():
            self.units[key] = self.units[key] * cm
        self.units = defaultdict(lambda: r"Sh", self.units)


class FigSettings:
    dpi = 512
    timeStampLoc = (0.6,0.9) # % of the x and y axes
    labelLoc     = (0.6,0.03)
    barlength    = (10, "pc")
    barLoc       = (0.1,0.1)
    linewidth    = 2
    markersize    = 20
    def __init__(self):
        pass
    
    
class Text:
    
    def __init__(self, name = "text", fontsize = "8", font = "Computer Modern", color = "black"):
        self.name     = "text"
        self.fontsize = 8
        self.font     = 'Times New Roman' #{'fontname':'Times New Roman'}#"Computer Modern"
        self.color    = "black"
        pass

class FontSettings:
    body = Text()
    
    def __init__(self):
        pass


class PaperSettings:
    onecol  = 3.5  #inches
    twocol  = 7.16 #inches
    pagelength = pt_to_in(705.0)#inches
    pagewidth  = pt_to_in(523.5307)
    psize = [
        (onecol, onecol),
        (onecol, pagewidth),
        (pagewidth, onecol*1.07),
        (onecol, pagelength*0.8),
        (onecol, onecol*0.5),
        
        ]
    figure  = FigSettings()
    font    = FontSettings()
    units   = Units()
    units2   = Units()
    unitsPj = UnitsP()
    cmaps   = cmaps()
    figfontstyle    = {'family': 'sans-serif', 'style': 'normal', 'weight': 'normal', 'size': 8}
    mplfigfontstyle = {'font.family': 'sans-serif', 'font.style': 'normal', 'font.size': 8}
    safecolor  = ['#00429d', '#315ba8', '#4f74b2', '#6b8eb8', '#8ba8bb', '#ffcab9', '#fd9291', '#e75d6f', '#c52a52', '#93003a']
    #colors     = ( "#352bf2", "#027295", "#04867b", "#1e9912", "#a08b01",  "#e27601", "#ff547b", "#eb08fb", "#7a7fff", "#0298cc", "#039b90", "#369a0f")
    niceColors = ("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")#["crimson","black","midnightblue", "magenta"]
  
    
    def __init__(self):
        pass 




