from pylab import *
from matplotlib import animation
from clawpack.visclaw import animation_tools
from scipy.signal import envelope
from clawpack.clawutil.util import fullpath_import
import os,sys
from scipy.interpolate import interp1d

root_dir = '/Users/rjl/git/Forks/multifluid-dispersive-wave-collab/crater_code/5eqns_transport/'
AGT = os.environ['AGT']
LW = fullpath_import(os.path.join(AGT, 'linear_transforms','linear_waves.py'))
mfclaw_tools = fullpath_import(f'{root_dir}/mfrk2Dev/mfclaw_tools.py')

h0 = 4000.
RC = 600

kmax = 0.06  # min wavelength 2*pi/kmax ~ 104.
mk = 2000
k = linspace(1e-6,kmax,mk)

fname_nc = f'surface{RC}m44finer.nc'
print(f'Loading surface_data from {fname_nc}')
tf_mfclaw, find_frame_mfclaw, surface_data \
                        = mfclaw_tools.load_surface_nc(fname_nc)

# initial time for Airy:
t0airy =  60
frameno0, t0frame = find_frame_mfclaw(time=t0airy)
print(f'Using mfclaw frame {frameno0} at time {t0frame} for t0airy={t0airy}')

rvals = surface_data.coords['r'].to_numpy()
eta0vals = surface_data.sel(t=t0frame).to_numpy()

r0 = rvals # may be extended below

eta0fcn = interp1d(rvals, eta0vals, fill_value=0., bounds_error=False)
eta0 = eta0fcn(r0)  # sample on r grid, extending by 0 if necessary for large r

dr = rvals[-1] - rvals[-2]

if 1:
    # damp out near origin and/or for large r:

    r1damp = 1500.
    wdamp = 0.3*r1damp   #e-folding width
    eta0 = where(r0<r1damp, eta0*exp((r0-r1damp)/wdamp), eta0)

    eta2 = eta0vals[-1]
    r2 = rvals[-1]
    slope2 = (eta0vals[-1] - eta0vals[-5]) / (rvals[-1] - rvals[-5])
    dr = rvals[-1] - rvals[-2]
    rextend = []
    rint = r2 - eta2/slope2
    print(f'slope2 = {slope2}, intersects at {rint}')
    if rint > r2:
        if rint > 1.1*r2:
            print(f'*** rint = {rint:.0f} too large, restricing to 1.1*r2')
            rint = 1.1*r2
        rnew = arange(r2+dr, rint, dr)
        print(f'Adding {len(rnew)} additional points to r')
        r0 = hstack((rvals, rnew))
        eta0new = eta2 + (rnew-r2)*slope2
        eta0 = hstack((eta0, eta0new))
    else:
        # extrapolation doesn't work, so damp out
        r2damp = 0.9*rvals[-1] # damp for r > r2damp
        w2damp = 0.3*0.1*rvals[-1]
        eta0 = where(r0>r2damp, eta0*exp((r2damp-r0)/w2damp), eta0)


print(f'with dr = {dr}, r0 has {len(r0)} points up to {r0[-1]}')


print('Computing eta0hat transform...')
print(f'  max(abs(eta0)) = {abs(eta0).max()}')
eta0hat = LW.Htransform(r0,eta0,k)
print(f'  max(abs(eta0hat)) = {abs(eta0hat).max()}')
#omega = lambda k: LW.omega_airy(k,h0)


tfinal = 1800.
dt = 120.
t = arange(t0airy, tfinal+0.9*dt, dt)

# grid for evaluation:
#rupper = 2*RC + sqrt(g*h0)*t[-1]
dr = 10.
rupper = 100e3
r = arange(dr/2, rupper, dr)
print(f'with dr = {dr}, r has {len(r)} points up to {r[-1]}')

eta_u_rt = LW.evolve_airy_eta_u(t, r, t0airy, r0, eta0, k, h0, eta0hat=eta0hat,
                             direction='outgoing', fname_nc=None)

if 1:
    fname_nc = f'airy_eta_u_RC{RC}_t{t0airy:.0f}-{tfinal:.0f}.nc'
    eta_u_rt.to_netcdf(fname_nc)
    print('Created ',fname_nc)
