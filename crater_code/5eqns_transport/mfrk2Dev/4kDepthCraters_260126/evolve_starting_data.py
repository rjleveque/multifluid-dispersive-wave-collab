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

import os,sys
#import linear_waves as LW
from scipy.interpolate import interp1d

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
xupper = 10e3
h0 = 4000.

# Initial eta and u
L = xupper - xlower
mx_H1 = 2000
dx = L/mx_H1
r = linspace(dx/2, xupper-dx/2, mx_H1)
RC = 600  # crater radius
kmax = 1.5 * 2*pi/RC
mk = mx_H1 * 2  # better choice?
k = linspace(1e-6,kmax,mk)

if 0:
    # wave packet:
    width = 10e3; r0 = 40e3; wavelength = 5e3; ampl = 100.
    eta0 = ampl*exp(-((r-r0)/width)**2) * cos(r*2*pi/wavelength)

#starting_data = loadtxt('starting.data',skiprows=2)

outdir = f'_output{RC}m44finer'
#tf_mfclaw, find_frame_mfclaw = C.load_times_mfclaw(outdir)
tf_mfclaw, find_frame_mfclaw = mfclaw_tools.load_times_mfclaw(outdir)

# initial time for Airy:
t0airy = 60
frameno0, t0frame = find_frame_mfclaw(time=t0airy)
print(f'Using mfclaw frame {frameno0} at time {t0frame} for t0airy={t0airy}')
#rkm, eta0, t0 = C.load_surf_mfclaw(outdir, j)

rvals, eta0vals = mfclaw_tools.load_surface(frameno0, outdir)

#r = rvals
rkm = r/1e3

eta0fcn = interp1d(rvals, eta0vals, fill_value=0., bounds_error=False)
eta0 = eta0fcn(r)  # sample on r grid, extending by 0 if necessary for large r


if 1:
    # damp out near origin and/or for large r:
    #eta0 = where(r>35e3, eta0*exp(-(r-35e3)/5e3), eta0)
    r1damp = 0.2e3
    wdamp = 0.5*r1damp   #e-folding width
    eta0 = where(r<r1damp, eta0*exp((r-r1damp)/wdamp), eta0)

figure(figsize=(9,5))
plot(rvals/1e3, eta0vals, 'r', label=f'eta0vals at t={t0airy:.0f} from mfclaw')
plot(r/1e3, eta0, 'b--',label=f'eta0 used for Airy')
xlim(0,1.1*xupper/1e3)
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
mfclaw_plot, = plot(rvals/1e3, eta_mfclaw, 'r', label='mfclaw')
grid(True)
xlabel('distance from crater (km)')
ylabel('surface elevation (m)')
rmax_plot = RC*20
xlim(0,rmax_plot/1e3)
ylim(-RC/6,RC/6)
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


def update(t):
    """
    Update an existing plot by resetting the data.
    Use the frames from each simulation closest to the given time t.
    """

    if t > t0airy:
        # mfluid:
        frameno_mfclaw, tframe_mfclaw = find_frame_mfclaw(t)
        if tframe_mfclaw != t:
            print('*** tframe_mfclaw does not agree with t')

        r_mfclaw, eta_mfclaw = mfclaw_tools.load_surface(frameno_mfclaw, outdir,
                                                         rmax=rmax_plot)
        mfclaw_plot.set_data(r_mfclaw/1e3, eta_mfclaw)

        title_text.set_text(f'Surface at t = {t:6.1f}')

        # Airy:
        print('Evaluating eta(r,t=%.0f) with Airy...' % t)
        omega = lambda k: LW.omega_airy(k,h0)

        print(f'+++ t = {t}, eta0hat.max() = {abs(eta0hat).max()}')
        eta_airy,u_airy = LW.eta_u_radial(t-t0airy,reval,k,eta0hat,omega,h0,
                                          direction='outgoing')
        #import pdb; pdb.set_trace()
        airy_plot.set_data(reval/1e3, eta_airy)


if __name__ == '__main__':

    print('Making anim...')
    #times = tf_mfclaw[:,1]
    times = arange(t0airy,121,10)
    times = [60, 120]
    anim = animation.FuncAnimation(fig, update, frames=times,
                                   interval=200, blit=False)

    print(f'+++ done making anim = {anim}')

    # Output files:
    #name = 'Gaussian_AirySwitch_t%s' % str(t0airy).zfill(3)
    name = f'RC{RC:04d}_AirySwitch_t{t0airy:03d}'

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
