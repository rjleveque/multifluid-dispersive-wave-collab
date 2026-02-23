
from pylab import *

import os,sys
#import linear_waves as LW
from scipy.interpolate import interp1d
from clawpack.clawutil.util import fullpath_import

# hardwired for now...
root_dir = '/Users/rjl/git/Forks/multifluid-dispersive-wave-collab/crater_code/5eqns_transport/'
mfclaw_tools = fullpath_import(f'{root_dir}/mfrk2Dev/mfclaw_tools.py')

AGT = os.environ['AGT']
LW = fullpath_import(os.path.join(AGT, 'linear_transforms','linear_waves.py'))

grav = 9.81

h0 = 4000.
RC = 600.
t0airy = 180.
tAiry_start = t0airy   # initial eta(r,t) from mfluid solution at this time
tAiry_end = 7200.      # compute up to this time
rAiry_end = 100e3      # capture the solution eta(r,t) as fcn of t at this r



# directory for saving boundary conditions time series and plots:
savefile_dir = f'./BC_RC{int(RC):04d}_h{int(h0):04d}_' \
                + f'tstart{int(tAiry_start)}'

os.system('mkdir -p %s' % savefile_dir)

# file for saving time series of eta and u (for SWE and Bouss) at rAiry_end:
# (will be saved in directory savefile_dir)
savefile_name = f'eta_hu_bc2_RC{int(RC):04d}_h{int(h0):04d}_' \
                + f'tstart{int(tAiry_start)}_rbc{int(rAiry_end/1000)}km'


# Hankel transform wave numbers k:
mk = 9000  # Is this large enough?  How to choose?
#kmax = 1.5 * 2*pi/RC
kmax = 0.02
k = linspace(1e-6,kmax,mk)

# Create initial data eta0 for Airy solution based on mfclaw surface
# at t = t0airy, possibly damped out near origin and/or for large r...

# where to obtain mfclaw results:

if 0:
    # load directly from fort.b files and extract surface
    outdir = f'_output{RC}m44finer'
    print(f'Loading solutions and extracting surface from {outdir}')
    tf_mfclaw, find_frame_mfclaw = mfclaw_tools.load_times_mfclaw(outdir)
    surface_data = None
else:
    # use surface_data previously extracted from fort.b files
    fname_nc = f'surface{int(RC):d}m44finer.nc'
    print(f'Loading surface_data from {fname_nc}')
    tf_mfclaw, find_frame_mfclaw, surface_data = \
                                    mfclaw_tools.load_surface_nc(fname_nc)

# when to switch from mfclaw solution to Airy:
frameno0, t0frame = find_frame_mfclaw(time=t0airy)
print(f'Using mfclaw frame {frameno0} at time {t0frame} for t0airy={t0airy}')


if surface_data is not None:
    rvals = surface_data.coords['r'].data  # locations where surface was evaluated
    eta0vals = surface_data.sel(t=t0frame).data
else:
    # load full fort.b files and choose resolution of rvals based on amr data:
    rvals, eta0vals = mfclaw_tools.load_surface(frameno0, outdir)

# set r, values to use in computing Hankel transform eta0hat of eta0:
r = rvals  # may need to extend if eta0 isn't 0 near rvals[-1]?

eta0fcn = interp1d(rvals, eta0vals, fill_value=0., bounds_error=False)
eta0 = eta0fcn(r)  # sample on r grid, extending by 0 if necessary for large r


if 1:
    # damp out near origin and/or for large r:
    #eta0 = where(r>35e3, eta0*exp(-(r-35e3)/5e3), eta0)
    #r1damp = 0.2e3
    #wdamp = 0.5*r1damp   #e-folding width
    r1damp = 1500.  # damp for 0 <= r < r1damp
    w1damp = 0.15*r1damp   # e-folding width
    eta0 = where(r<r1damp, eta0*exp((r-r1damp)/w1damp), eta0)

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
        r = hstack((rvals, rnew))
        eta0new = eta2 + (rnew-r2)*slope2
        eta0 = hstack((eta0, eta0new))
    else:
        # extrapolation doesn't work, so damp out
        r2damp = 0.9*rvals[-1] # damp for r > r2damp
        w2damp = 0.3*0.1*rvals[-1]
        eta0 = where(r>r2damp, eta0*exp((r2damp-r)/w2damp), eta0)

figure(figsize=(9,5))
plot(rvals/1e3, eta0vals, 'r', label=f'eta0vals at t={t0airy:.0f} from mfclaw')
plot(r/1e3, eta0, 'b--',label=f'eta0 used for Airy')
xlim(0,max(rvals[-1],r[-1])/1e3)
grid(True)
legend(framealpha=1)
#fname = f'{savefile_dir}/eta0_t{t0airy:03d}.png'
title(f'eta0 for RC={RC:.0f}m, h0={h0:.0f}m')
fname = f'{savefile_dir}/{savefile_name.replace('eta_hu_bc', 'eta0')}.png'
savefig(fname)
print('Created ',fname)

print('Computing eta0hat transform...')
print(f'  max(abs(eta0)) = {abs(eta0).max()}')
eta0hat = LW.Htransform(r,eta0,k)
print(f'  max(abs(eta0hat)) = {abs(eta0hat).max()}')
omega = lambda k: LW.omega_airy(k,h0)
reval = r  # grid to evaluate Airy solution

fname = f'{savefile_dir}/{savefile_name}.txt'
t,eta,hu_swe,hu_sgn = mfclaw_tools.make_bc(rAiry_end, tAiry_end, tAiry_start,
                              eta0hat, k, h0, t=None, savefile=fname)

fname = f'{savefile_dir}/{savefile_name}.png'
mfclaw_tools.plot_bc(t,eta,hu_swe,hu_sgn,rAiry_end,tAiry_end,
                    tAiry_start,RC,h0,savefile=fname)
