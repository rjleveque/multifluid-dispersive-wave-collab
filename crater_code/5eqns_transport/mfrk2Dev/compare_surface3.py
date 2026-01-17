from pylab import *
from clawpack.pyclaw import Solution
from clawpack.visclaw import gridtools
import glob

# Need to adjust:
root_dir = '/Users/rjl/git/Forks/multifluid-dispersive-wave-collab/crater_code/5eqns_transport/'

try:
    # should be in recent Clawpack >= v5.13.0:
    from clawpack.clawutil.util import fullpath_import
except:
    print('*** Did not find fullpath_import in clawpack, defining here...')
    def fullpath_import(fullpath):
        """
        Return a module imported from a full path name.
        To reload the module, call this again (rather than using importlib.reload).
        """

        import os, sys, importlib
        fname = os.path.split(fullpath)[1]
        modname = os.path.splitext(fname)[0]
        spec = importlib.util.spec_from_file_location(modname, fullpath)
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        sys.modules[modname] = module
        print('loaded module from file: ',module.__file__)
        return module

mfclaw_tools = fullpath_import(f'{root_dir}/mfrk2Dev/mfclaw_tools.py')

# ==================
# Things to adjust:

interactive = True
if interactive:
    ion()

save_figs = not interactive  # or set to True to save figs when interactive

# define radial values r on which to compute surface:
rmax = 2e3

# y-axis limits for plots:
zmin = -4000
zmax = 1000


# rvals and zvals are now set in load_surface based on solution:

# outdirs to compare:

outdir1 = 'm4to6k300mcrater_250by125_3lev44/_output'
outdir2 = 'm4to6k300mcrater_250by125_3lev64/_output'
outdir3 = 'm4to6k300mcrater_500by250_3lev44/_output'
outdir4 = 'm4to6k300mcrater_250by125_3lev66/_output'

# list of outdirs to iterate over in each frame plot:
outdirs = [outdir1, outdir2,outdir3,outdir4]

# RJL testing:
outdir1 = '_output_1level'
outdir2 = '_output_2levels'
outdir3 = '_output_3levels'
outdirs = [outdir1,outdir2]

#labels = ['6k atm','7k atm','8k atm','finer 6k']  # labels to appear in legend
labels = 'dr values'  # indicates the calculated values should be used

#drs = [6.6667,6.6667,6.6667, 2.2222]  # not used - dr now calculated
colors = ['b','r','g','m','k']  # color to plot for each outdir

# ==================


# load framenos and times from outdir1, assuming the others are the same:
frame_time, find_frame = mfclaw_tools.load_times_mfclaw(outdir1)

# Note the find_frame function returned can be used to find the frame
# closest to a given time, if they aren't aligned between outdir's

# first column of frame_time is frameno, second column is time,
# make dictionary of times for each frameno:
framenos = frame_time[:,0]
t_outdir1 = {}
for k,frameno in enumerate(framenos):
    t_outdir1[frameno] = frame_time[k,1]


ans = ''
frameno = 0

if interactive:
    ans = input('Return for first frame, integer to jump, or q to quit... ')
    try:
        frameno = int(ans)
        print(f'Jumping to frameno = {frameno}')
    except:
        frameno = 0



while ans != 'q':
    if frameno not in framenos:
        print(f'*** frameno={frameno} is not in framenos found in {outdir1}')
        break
    figure(200, figsize=(10,6))
    clf()
    for k,outdir in enumerate(outdirs):
        try:
            rvals, eta1 = mfclaw_tools.load_surface(frameno, outdir)
        except:
            print(f"*** Could not load frame {frameno} from outdir={outdir}")
            raise

        dr = rvals[1] - rvals[0]
        if labels == 'dr values':
            plot(rvals, eta1, colors[k], label=f'dr = {dr:.2f}')
        else:
            plot(rvals, eta1, colors[k], label=labels[k])


    legend(loc='upper right', framealpha=1)
    grid(True)
    xlim(0,rmax)
    ylim(zmin,zmax)
    title(f'Frame {frameno} at time t = {t_outdir1[frameno]:.2f} seconds')
    draw()

    if save_figs:
        #fname = 'compHeightAtmos300m%s.png'%str(frameno).zfill(4)
        fname = 'compRes300m%s.png'%str(frameno).zfill(4)
        savefig(fname)
        print("created ",fname)

    if interactive:
        # interactive prompt:
        ans = input('Return for next frame, integer to jump, or q to quit... ')
        try:
            frameno = int(ans)
            print(f'Jumping to frameno = {frameno}')
        except:
            frameno += 1
    else:
        # always advance to next frame to make plots without prompting
        frameno += 1
