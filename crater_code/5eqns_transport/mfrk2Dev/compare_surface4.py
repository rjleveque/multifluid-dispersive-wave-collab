from pylab import *
from clawpack.pyclaw import Solution
from clawpack.visclaw import gridtools
import glob


from clawpack.clawutil.util import fullpath_import
mfclaw_tools = fullpath_import('./mfclaw_tools.py') # fix path if needed

# ==================
# Things to adjust:

interactive = True
if interactive:
    ion()

save_figs = not interactive  # or set to True to save figs when interactive

# x-axis limit for plots:
rmax = 2e3

# y-axis limits for plots:
zmin = -1000
zmax = 1000

colors = ['b','r','g','m','k']  # color to plot for each outdir

# labels distinguish different runs and are used as dictionary keys and
# for labels in plot legend

# For each label, set fname_nc[label] to the .nc file containing all the
# surfaces previously extracted using  extract_all_surfaces.py

fname_nc = {}

label1 = '1 level'
fname_nc[label1] = 'surface_1level.nc'

label2 = '2 levels'
fname_nc[label2] = 'surface_2levels.nc'

labels = [label1, label2]

# ==================

surface_data = {}
rvals = {}
tvals = {}
eta = {}

for label in labels:
    # we won't use tf or find_frame, so don't need to save these in dict's,
    # just the surface_data (which is an xarray.DataArray that contains
    # the coordinates r,t and a 2D array of values eta)
    tf, find_frame, surface_data[label] = \
            mfclaw_tools.load_surface_nc(fname_nc[label])


for label in labels:
    # convert DataArray info to ordinary numpy arrays:
    tvals[label] = surface_data[label].coords['t'].to_numpy()
    rvals[label] = surface_data[label].coords['r'].to_numpy()
    eta[label] = surface_data[label].to_numpy()

    # eta should have one column for each t value:
    assert eta[label].shape[1] == len(tvals[label]), '*** unexpected shape'
    
    print(f'{fname_nc[label]}: {len(tvals[label])} frames with ' \
      + f'{len(rvals[label])} radial values up to r = {rvals[label][-1]:.1f}m')

# Assume all data sets have the same number of frames at same times
framenos = range(len(tvals[labels[0]]))


ans = ''
frameno = 0

if interactive:
    ans = input('Return for first frame, integer to jump, or q to quit... ')
    try:
        frameno = int(ans)
        print(f'Jumping to frameno = {frameno}')
    except:
        frameno = 0

# Main loop on frames:

while ans != 'q':
    if frameno not in framenos:
        print(f'*** frameno={frameno} is too large')
        break
    figure(200, figsize=(10,6))
    clf()
    for k,label in enumerate(labels):
        plot(rvals[label], eta[label][:,frameno], colors[k], label=label)

    legend(loc='upper right', framealpha=1)
    grid(True)
    xlim(0,rmax)
    ylim(zmin,zmax)
    title(f'Frame {frameno} at time t = {tvals[label][frameno]:.2f} seconds')
    draw()

    if save_figs:
        fname = 'compareRes%s.png'%str(frameno).zfill(4)
        savefig(fname)
        print("created ",fname)

    if interactive:
        # interactive prompt:
        ans = input(f'Showing Frame {frameno}  ' \
                + 'Return for next frame, integer to jump, or q to quit... ')
        try:
            frameno = int(ans)
            print(f'Jumping to frameno = {frameno}')
        except:
            frameno += 1
    else:
        # always advance to next frame to make plots without prompting
        frameno += 1
