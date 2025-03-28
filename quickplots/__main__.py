
from utils.myparser import myparser

parser = myparser(prog = "Fastplosts") #argparse.ArgumentParser(description='Runs and anayses output from TurbGen.')
parser.add_argument("--makemovie", default=False) #argparse.ArgumentParser(description='Runs and anayses output from TurbGen.')
subparsers = parser.add_subparsers(title='sub-commands', dest='subcommand', help='see additional help with -h')

    # sub-command 'generate'
parser_phase    = subparsers.add_parser('phase', help="Call phase diagram maker.")
parser_plot     = subparsers.add_parser('plot', help="Call phase diagram maker.")
parser_plot.conflict_handler  = "resolve"
parser_hist     = subparsers.add_parser('hist', help="Call phase diagram maker.")
parser_profile  = subparsers.add_parser('profile', help="Call phase diagram maker.")
parser_fastplot = subparsers.add_parser('fastplot', help="Call phase diagram maker.")
parser_fastprofile = subparsers.add_parser('fastprofile', help="Call phase diagram maker.")
parser_movie = subparsers.add_parser('movie', help="Call phase diagram maker.")


#args = parser.parse_args()
args, subc_args = parser.parse_known_args(fixdir = False, fix_n1n2 = False,)
# print("Calling '{}'".format(args.subcommand))

# print("args '{}'".format(args))
# print("subargs '{}'".format(subc_args))


def moviecall(moviepath = None ):
        import os            
        parser_movie.add_argument("-path", type=str, default="./")            
        parser_movie.add_argument("-o", "--output", type=str, default="output.mp4")
        parser_movie.add_argument("-cpus", type=int, default=None)
        parser_movie.add_argument("-duration", type=str, default="0.06666")
        argss = parser_movie.parse_args(subc_args, fixdir=False )
        if moviepath is None: moviepath = argss.path 
        command = f"for i in {moviepath}/*.png; do echo file \\'$i\\'; echo duration {argss.duration}; done > {moviepath}/list.txt"
        
        done = os.system(command)

        command = f"yes | ffmpeg -f concat  -safe 0 -i {moviepath}/list.txt -quality best -c:v libx264  {argss.output}"
        if argss.cpus: command += f" -threads {argss.cpus}"
        done = os.system(command)
        if not done: quit()
        #if not done: raise Exception("Problem with movie making.")


match  args.subcommand:
    case 'phase':
        from phase import MakePhase

        MakePhase(parser_phase, subc_args)
    case 'plot':
        MakePlot(parser_plot, subc_args)
    case 'hist':
        from quickhist import MakeHist
        MakeHist(parser_hist, subc_args)

    case 'profile':
        from quickplot import MakePlot
        MakePlot(parser_profile, subc_args)
    case 'fastplot':
        
        from fastplots_scrpt import MakeFastPlot

        MakeFastPlot(parser_fastplot, subc_args)
    case 'fastprofile':
        from fast_profiles import MakeFastProfile
        
        MakeFastProfile(parser_fastprofile, subc_args)
    case 'movie':
        moviecall()
    case _:
        print("No subcommand given.")

        msg = "".join(["\n\t\t%s"%sc for sc in subparsers.choices.keys()])
        print(f"Use one of the following: {msg}")
        
        quit()
