"""

"""

from pylab import *
from matplotlib import animation
from clawpack.clawutil.util import fullpath_import
import xarray as xr

import os,sys

root_dir = '/Users/rjl/git/Forks/multifluid-dispersive-wave-collab/crater_code/5eqns_transport/'

mfclaw_tools = fullpath_import(f'{root_dir}/mfrk2Dev/mfclaw_tools.py')

RC = 600  # crater radius
outdir = f'_output{RC}m44finer'
name = f'surface{RC}'

tf_mfclaw, find_frame_mfclaw = mfclaw_tools.load_times_mfclaw(outdir)

nframes = tf_mfclaw.shape[0]
eta = None
for k in range(nframes):
    frameno = tf_mfclaw[k,0]
    tframe = tf_mfclaw[k,1]
    r_k, eta_k = mfclaw_tools.load_surface(k, outdir)
    if eta is None:
        eta = zeros((len(r_k), nframes))
    eta[:,k] = eta_k

if 0:
    fname1 = f'{name}_tf.txt'
    savetxt(fname1, tf_mfclaw, fmt='%11.3f', header='frameno, time')
    print('Created ',fname1)

    r_eta = vstack((r_k, eta)).T

    fname = f'{name}_r_eta.npy'
    save(fname, r_eta)
    print('Created ',fname)

surface_data = xr.DataArray(eta, dims=('r','t'),
                        coords={'r':r_k, 't':tf_mfclaw[:nframes,1]})

fname_nc = f'{name}m44finer.nc'
surface_data.to_netcdf(fname_nc)
print('Created ',fname_nc)

# reload and plot at t=0 via:
#    surface_data = xr.open_dataarray(fname_nc)
#    plot(surface_data.coords['r'], surface_data.sel(t=0.))
# use array(...) to convert to numpy arrays
