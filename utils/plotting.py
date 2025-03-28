import matplotlib.pyplot as plt
import string
import paper_settings as pp 

paper    = pp.PaperSettings()
alphabet = list(string.ascii_lowercase)


def set_panel( ax, ind, color = "black"): 
    
    txt = ind if isinstance(ind, str) else f"{alphabet[ind]})" 
        
    return ax.text(0.9, 0.1, txt, fontsize = paper.font.body.fontsize, color = color, transform = ax.transAxes)


