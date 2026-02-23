"""
Comparison of mfclaw results with Airy solution using data from Marsha

Using AMR results computed with mfclaw and Hankel function Airy solution.
"""

import sys
if 'matplotlib' not in sys.modules:
    # Use an image backend to insure animation has size specified by figsize
    import matplotlib
    matplotlib.use('Agg')

from pylab import *
from matplotlib import animation
from clawpack.visclaw import animation_tools
from scipy.signal import envelope
from clawpack.clawutil.util import fullpath_import
import xarray as xr
import os,sys


root_dir = '/Users/rjl/git/Forks/multifluid-dispersive-wave-collab/crater_code/5eqns_transport/'

AGT = os.environ['AGT']
LW = fullpath_import(os.path.join(AGT, 'linear_transforms','linear_waves.py'))
mfclaw_tools = fullpath_import(f'{root_dir}/mfrk2Dev/mfclaw_tools.py')


# use evolution data computed previously
RC = 600
fname_nc = 'airy_eta_u_RC600_t60-1800.nc'
print(f'Loading evolution data from {fname_nc}')
eta_u_rt = xr.open_dataset(fname_nc)

tvals = eta_u_rt.coords['t'].to_numpy()
rvals = eta_u_rt.coords['r'].to_numpy()
rkm = rvals / 1e3

rkm_max_plot = 100
ylimits = (-20,20)

# initial time for Airy:
t0airy = tvals[0]
eta0 = eta_u_rt['eta'].sel(t=t0airy).to_numpy()


# Create intial plot:

fig = figure(figsize=(12,7))

airy_plot, = plot(rkm, eta0, 'b',
         label='Airy, starting from mfclaw at t=%.0f' % t0airy)
grid(True)
xlabel('distance from crater (km)')
ylabel('surface elevation (m)')
xlim(0,rkm_max_plot)
ylim(ylimits)
grid(True)
title_text = title(f'Surface at t = {t0airy:6.1f}')

legend(loc='upper right', framealpha=1)


def update(t):
    """
    Update an existing plot by resetting the data.
    Use the frames from each simulation closest to the given time t.
    """

    if t in tvals:
        eta_airy = eta_u_rt['eta'].sel(t=t).to_numpy()
    else:
        print(f'*** Skipping time {t}, not in eta_u_rt')
        return

    airy_plot.set_data(rkm, eta_airy)
    title_text.set_text(f'Surface at t = {t:6.1f}')


if __name__ == '__main__':

    tmax = 1200
    times = [tv for tv in tvals if tv <= tmax]
    print(f'Making anim with {len(times)} frames, up to time {times[-1]}...')


    anim = animation.FuncAnimation(fig, update, frames=times,
                                   interval=200, blit=False)

    print(f'+++ done making anim = {anim}')

    # Output files:
    name = f'RC{RC:04d}_Airy_t{t0airy:.0f}-{times[-1]:.0f}'

    fname_mp4 = name + '.mp4'

    fname_html = None
    #fname_html = name + '.html'

    if fname_mp4:
        fps = 5
        print('Making mp4...')
        writer = animation.writers['ffmpeg'](fps=fps)
        anim.save(fname_mp4, writer=writer)
        print("Created %s" % fname_mp4)

    if fname_html:
        # html version:
        animation_tools.make_html(anim, file_name=fname_html, title=name)
