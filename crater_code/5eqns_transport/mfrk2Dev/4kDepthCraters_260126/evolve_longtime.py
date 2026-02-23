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

import os,sys
#import linear_waves as LW
from scipy.interpolate import interp1d

global eta_bc, etamax

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

root_dir = '/Users/rjl/git/Forks/multifluid-dispersive-wave-collab/crater_code/5eqns_transport/'

AGT = os.environ['AGT']
LW = fullpath_import(os.path.join(AGT, 'linear_transforms','linear_waves.py'))
mfclaw_tools = fullpath_import(f'{root_dir}/mfrk2Dev/mfclaw_tools.py')


xlower = 0
#xupper = 10e3 # for RC=150
xupper = 40e3
h0 = 4000.

# Initial eta and u
L = xupper - xlower
mx_H1 = 2000
dx = L/mx_H1
r = linspace(dx/2, xupper-dx/2, mx_H1)
RC = 600  # crater radius
#kmax = 1.5 * 2*pi/RC
kmax = 0.2
mk = mx_H1 * 2  # better choice?
k = linspace(1e-6,kmax,mk)

if 0:
    # wave packet:
    width = 10e3; r0 = 40e3; wavelength = 5e3; ampl = 100.
    eta0 = ampl*exp(-((r-r0)/width)**2) * cos(r*2*pi/wavelength)

if 0:
    # load directly from fort.b files and extract surface
    outdir = f'_output{RC}m44finer'
    print(f'Loading solutions and extracting surface from {outdir}')
    tf_mfclaw, find_frame_mfclaw = mfclaw_tools.load_times_mfclaw(outdir)
    surface_data = None
else:
    # use surface_data previously extracted from fort.b files
    fname_nc = f'surface{RC}m44finer.nc'
    print(f'Loading surface_data from {fname_nc}')
    tf_mfclaw, find_frame_mfclaw, surface_data = mfclaw_tools.load_surface_nc(fname_nc)



# initial time for Airy:
t0airy =  30
frameno0, t0frame = find_frame_mfclaw(time=t0airy)
print(f'Using mfclaw frame {frameno0} at time {t0frame} for t0airy={t0airy}')
#rkm, eta0, t0 = C.load_surf_mfclaw(outdir, j)

if surface_data is not None:
    rvals = surface_data.coords['r'].data
    eta0vals = surface_data.sel(t=t0frame).data
else:
    rvals, eta0vals = mfclaw_tools.load_surface(frameno0, outdir)

r = rvals # may be extended below

eta0fcn = interp1d(rvals, eta0vals, fill_value=0., bounds_error=False)
eta0 = eta0fcn(r)  # sample on r grid, extending by 0 if necessary for large r


if 0:
    # damp out near origin and/or for large r:
    #eta0 = where(r>35e3, eta0*exp(-(r-35e3)/5e3), eta0)
    #r1damp = 0.2e3
    #wdamp = 0.5*r1damp   #e-folding width
    r1damp = 1500.
    wdamp = 0.3*r1damp   #e-folding width
    eta0 = where(r<r1damp, eta0*exp((r-r1damp)/wdamp), eta0)
    eta0 = where(r<r1damp, eta0*exp((r-r1damp)/wdamp), eta0)

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

rkm = r/1e3


figure(figsize=(9,5))
plot(rvals/1e3, eta0vals, 'r', label=f'eta0vals at t={t0airy:.0f} from mfclaw')
plot(r/1e3, eta0, 'b--',label=f'eta0 used for Airy')
#xlim(0,1.1*xupper/1e3)
xlim(0, 100.)
fname = f'eta0_t{t0airy:03d}.png'
savefig(fname)
print('Created ',fname)



print('Computing eta0hat transform...')
print(f'  max(abs(eta0)) = {abs(eta0).max()}')
eta0hat = LW.Htransform(r,eta0,k)
print(f'  max(abs(eta0hat)) = {abs(eta0hat).max()}')
omega = lambda k: LW.omega_airy(k,h0)
reval = r  # grid to evaluate Airy solution

if 1:
    figure(33, figsize=(9,5))
    clf()
    plot(k, eta0hat, 'k')
    xlabel('wavenumber k')
    ylabel('eta0hat')
    grid(True)
    title('Hankel transform of eta0')
    fname = 'Htransform.png'
    savefig(fname)
    print('Created ',fname)


# Create intial plot:

fig = figure(figsize=(12,7))


# mfclaw:
#frameno = find_frame_mfclaw(t)
#rkm_mfclaw, eta_mfclaw, t_mfclaw = C.load_surf_mfclaw(outdir,
#                                                    frameno_mfclaw)
#r_mfclaw, eta_mfclaw = mfclaw_tools.load_surface(frameno, outdir)


eta_mfclaw = eta0vals  # on grid rvals
eta_airy = eta0   # agrees with eta_mfclaw at t0airy (but on grid reval)
airy_plot, = plot(reval/1e3, eta_airy, 'b',
         label='Airy, starting from mfclaw at t=%.0f' % t0airy)
#mfclaw_plot, = plot(rvals/1e3, eta_mfclaw, 'r', label='mfclaw')
grid(True)
xlabel('distance from crater (km)')
ylabel('surface elevation (m)')
#rmax_plot = RC*20
rmax_plot = 100e3
xlim(0,rmax_plot/1e3)
#ylim(-RC/6,RC/6)
ylim(-20,20)
grid(True)
title_text = title(f'Surface at t = {t0airy:6.1f}')

if 0:
    print('Evaluating eta(r,t) with Airy...')
    t = 600.
    omega = lambda k: LW.omega_airy(k,h0)

    eta,u = LW.eta_u_radial(t-t0,reval,k,eta0hat,omega,h0)
    plot(reval/1e3, eta, 'b', label='Airy evolution to t=%.0f' % t)

    j = find_frame(time=t)
    r_mfclaw, eta_mfclaw = mfclaw_tools.load_surface(j, outdir)
    airy_plot, = plot(rkm, eta_mfclaw, 'b',
         label='Airy at t=%.0f starting from mfluid at t=%.0f' \
                % (t_mfclaw,t0airy))
    mfclaw_plot, = plot(rkm, eta_mfclaw, 'r', label='mfluid at t=%.0f' % t)

legend(loc='upper right', framealpha=1)

eta_bc = []
r_bc_km = None

etamax = zeros(eta_airy.shape)
r_etamax = reval.copy()
etamax_plot, = plot(r_etamax/1e3, etamax, 'r', label='envelope so far')

def update(t):
    """
    Update an existing plot by resetting the data.
    Use the frames from each simulation closest to the given time t.
    """
    global eta_bc, etamax
    #if t > t0airy:
    if 1:
        if 0:
            # mfluid:
            frameno_mfclaw, tframe_mfclaw = find_frame_mfclaw(t)
            if tframe_mfclaw != t:
                print('*** tframe_mfclaw does not agree with t')

            if surface_data is not None:
                r_mfclaw = surface_data.coords['r']
                eta_mfclaw = surface_data.sel(t=tframe_mfclaw)
            else:
                r_mfclaw, eta_mfclaw = mfclaw_tools.load_surface(frameno_mfclaw,
                                                        outdir, rmax=rmax_plot)

            mfclaw_plot.set_data(r_mfclaw/1e3, eta_mfclaw)

        title_text.set_text(f'Surface at t = {t:6.1f}')

        # Airy:
        print('Evaluating eta(r,t=%.0f) with Airy...' % t)
        omega = lambda k: LW.omega_airy(k,h0)

        #print(f'+++ t = {t}, eta0hat.max() = {abs(eta0hat).max()}')
        rmax = 20e3 + 80e3 * min(600, t) / 600
        dr = 20.
        reval = arange(dr,rmax,dr)
        print(f't = {t}, rmax = {rmax}, dr = {dr}, len(reval) = {len(reval)}')

        eta_airy,u_airy = LW.eta_u_radial(t-t0airy,reval,k,eta0hat,omega,h0,
                                          direction='outgoing')

        eta_airy = real(eta_airy)

        # computed integral of eta^2 * r*dr for PE in radial coordinates
        eta2r = eta_airy**2 * reval
        eta2integral = eta2r.sum()*dr
        print(f'integral of eta^2 * r*dr = {eta2integral:.2f}')

        if r_bc_km is not None:
            j = where(reval/1e3 < r_bc_km)[0].max()
            eta_bc.append([t,eta_airy[j]])

        eta_env = envelope(eta_airy)
        etamax[:len(reval)] = maximum(etamax[:len(reval)], eta_env[0,:])
        etamax_plot.set_data(r_etamax/1e3, etamax)

        airy_plot.set_data(reval/1e3, eta_airy)


if __name__ == '__main__':

    print('Making anim...')
    #times = tf_mfclaw[:,1]
    times = arange(t0airy,1801,300)
    #times = [t0airy, t0airy+20, t0airy+40]


    anim = animation.FuncAnimation(fig, update, frames=times,
                                   interval=200, blit=False)


    print(f'+++ done making anim = {anim}')

    # Output files:
    #name = 'Gaussian_AirySwitch_t%s' % str(t0airy).zfill(3)
    name = f'RC{RC:04d}_AirySwitch_t{t0airy:03d}_longtime'

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

    if r_bc_km is not None:
        eta_bc = array(eta_bc)
        figure()
        plot(eta_bc[:,0], eta_bc[:,1])
        title(f'eta_airy vs t at {r_bc_km}km')
        fname = f'eta_airy_bc_{r_bc_km}km.png'
        savefig(fname)
        print('Created ',fname)

    fname = name + '_env.txt'
    savetxt(fname, vstack((rvals, etamax)).T, fmt='%16.6f')
    print('Created ',fname)
