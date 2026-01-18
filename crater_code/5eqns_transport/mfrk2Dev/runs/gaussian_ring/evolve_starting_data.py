"""
Comparison of mfclaw results with Airy solution for initial data consisting
of a radial Gaussian ring.

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
xupper = 20e3
h0 = 1500.  # Note: also used to set outdir below

# Initial eta and u
L = xupper - xlower
mx_H1 = 2000
dx = L/mx_H1
r = linspace(dx/2, xupper-dx/2, mx_H1)
w_ring = 500
kmax = 1.5 * 2*pi/w_ring
mk = mx_H1 * 2  # better choice?
k = linspace(1e-6,kmax,mk)

if 0:
    # wave packet:
    width = 10e3; r0 = 40e3; wavelength = 5e3; ampl = 100.
    eta0 = ampl*exp(-((r-r0)/width)**2) * cos(r*2*pi/wavelength)

#starting_data = loadtxt('starting.data',skiprows=2)

outdir = f'mfclaw_outputs/_output_h0_{h0:04.0f}'
#tf_mfclaw, find_frame_mfclaw = C.load_times_mfclaw(outdir)
tf_mfclaw, find_frame_mfclaw = mfclaw_tools.load_times_mfclaw(outdir)

# initial time for Airy:
t0airy = 0
j = find_frame_mfclaw(time=t0airy)
#rkm, eta0, t0 = C.load_surf_mfclaw(outdir, j)
rvals, eta0 = mfclaw_tools.load_surface(j, outdir)

r = rvals
rkm = r/1e3

etafcn = interp1d(r, eta0, fill_value=0., bounds_error=False)
#plot(r/1e3,eta0,'k--',label='orig eta at t0=%.0f' % t0)

if 0:
    # damp out near origin and for large r:
    eta0 = etafcn(r)
    #eta0 = where(r>35e3, eta0*exp(-(r-35e3)/5e3), eta0)
    eta0 = where(r<4e3, eta0*exp((r-4e3)/1e3), eta0)

print('Computing eta0hat transform...')
eta0hat = LW.Htransform(r,eta0,k)
omega = lambda k: LW.omega_airy(k,h0)
reval = r

if 1:
    figure(33, figsize=(9,5))
    clf()
    plot(k, eta0hat, 'k')
    xlabel('wavenumber k')
    ylabel('eta0hat')
    grid(True)
    title('Hankel transform of eta0')


# Create intial plot:

fig = figure(figsize=(12,7))

t = 0.

# mfclaw:
frameno = find_frame_mfclaw(t)
#rkm_mfclaw, eta_mfclaw, t_mfclaw = C.load_surf_mfclaw(outdir,
#                                                    frameno_mfclaw)
r_mfclaw, eta_mfclaw = mfclaw_tools.load_surface(frameno, outdir)

rkm = r/1e3

eta_airy = nan*eta_mfclaw
airy_plot, = plot(rkm, eta_airy, 'b',
         label='Airy, starting from mfclaw at t=%.0f' % t0airy)
mfclaw_plot, = plot(rkm, eta_mfclaw, 'r', label='mfclaw')
grid(True)
xlabel('distance from crater (km)')
ylabel('surface elevation (m)')
xlim(0,10.)
#ylim(-1200,1200)
ylim(-150,150)
grid(True)
title_text = title('Gaussian,' \
            + '  t = %6.1f' % t)

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

    # mfluid:
    frameno_mfclaw = find_frame_mfclaw(t)

    r_mfclaw, eta_mfclaw = mfclaw_tools.load_surface(frameno_mfclaw, outdir)
    mfclaw_plot.set_data(rkm, eta_mfclaw)

    title_text.set_text('Gaussian,' \
            + '  t = %6.1f' % t)


    if t > t0airy:
        print('Evaluating eta(r,t=%.0f) with Airy...' % t)
        omega = lambda k: LW.omega_airy(k,h0)

        eta_airy,u_airy = LW.eta_u_radial(t-t0airy,reval,k,eta0hat,omega,h0,
                                          direction='both')
        airy_plot.set_data(reval/1e3, eta_airy)


if __name__ == '__main__':

    print('Making anim...')
    times = tf_mfclaw[:,1]
    #times = [0, 5,10]
    anim = animation.FuncAnimation(fig, update, frames=times,
                                   interval=200, blit=False)

    # Output files:
    #name = 'Gaussian_AirySwitch_t%s' % str(t0airy).zfill(3)
    name = f'Gaussian_Airy_depth{h0:04.0f}'

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
