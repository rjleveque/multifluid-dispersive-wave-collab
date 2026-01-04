"""
make anim from mfrk2Dev simulation
"""

import sys
if 'matplotlib' not in sys.modules:
    # Use an image backend to insure animation has size specified by figsize
    import matplotlib
    matplotlib.use('Agg')

from pylab import *
from matplotlib import animation
from scipy import io as sio
import os,sys,glob
import pathlib


root_dir = '/Users/rjl/git/Forks/multifluid-dispersive-wave-collab/crater_code/5eqns_transport/'

try:
    # should be in recent Clawpack >= v5.13.0:
    from clawpack.clawutil.util import fullpath_import
except:
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

# code to read Darrel's ALE3D sim's for comparison:
C = fullpath_import(f'{root_dir}/compareALE3D.py')

def load_surface(frameno, rvals, zvals, outdir='_output', file_format='binary'):
    """
    Read in the mfclaw AMR solution for this time frame, interpolating
    to the radial coordinate r values in rvals
    Requires zvals, vertical values used for integrating the indicator,
    function, chosen so zvals[0] and zvals[-1] always bracket the free
    surface, with dz small relative to the perturbation of the interface.
    """
    from clawpack.pyclaw import Solution
    from clawpack.visclaw import gridtools

    framesoln = Solution(frameno, path=outdir, file_format=file_format)

    print(f'Loaded solution at time {framesoln.t:.3f}, computing eta...')

    eta = zeros(len(rvals))  # initialize array

    for i,r in enumerate(rvals):
        ri = r*ones(len(zvals))
        # extract zfa = q[5,:,:] on vertical transect at r=ri:
        # (taking data from finest level that exists at each point in zvals)
        zfa = gridtools.grid_output_2d(framesoln, 5, ri, zvals, method='linear')
        have = (1. - zfa).sum() * dz
        eta[i] = have + zmin
        #import pdb; pdb.set_trace()
    return eta


include_ale = False

if include_ale:
    #datadir_ale = '/Users/rjl/D/Darrel_craters/surface_data_250905'
    datadir_ale = '/Users/rjl/D/Darrel_craters/surface_data_250910'
    fname_prefix_ale = 'hemi1000m_4km'
    tf_ale, find_frame_ale = C.load_times_ale(datadir_ale, fname_prefix_ale)

# ------------------------------------
# set parameters for this movie:

outdir_mfluid = f'{root_dir}/mfrk2Dev/_outputc'
title_base = 'Hemispherical crater with radius 300 m on 4km ocean'
movie_name = 'Test'
rmax = 20e3
dr = 100 #20e3 / (250*12)   # finest grid resolution from setrun.py
rvals = arange(dr/2,rmax,dr)

# vertical z values on which to evaluate 2D solution for computing surface:
zmin = -400. # any depth in domain that's always below surface
zmax = 1000.  # any elevation in domain that's always above surface
dz = 0.01  # needs to be pretty small
zvals = arange(zmin, zmax, dz)

# ------------------------------------

tf_mfluid, find_frame_mfluid = C.load_times_mfluid(outdir_mfluid)
#print(f'{len(tf_mfluid)} frames found in mfluid solution: ')
#print(f'    last frameno = {tf_mfluid[-1][0]} at time {tf_mfluid[-1][1]:.2f}')

rkm_mfluid = rvals / 1e3  # rvals in km

# Create intial plot:

fig = figure(figsize=(12,7))

t = 0.

# mfluid:
frameno_mfluid = find_frame_mfluid(t)

#eta_mfluid = load_surface(frameno_mfluid, rvals, zvals, outdir_mfluid)
# don't bother loading first frame properly since this plot is not used in anim
eta_mfluid = 200 + 0 * rkm_mfluid
mfluid_plot, = plot(rkm_mfluid, eta_mfluid, 'r', label='mfluid')

if include_ale:
    # ALE:
    frameno_ale = find_frame_ale(t)
    fname = '%s_surface_%s.mat' % (fname_prefix_ale, str(frameno_ale).zfill(4))
    rkm, eta, t = C.load_surf_ale(datadir_ale, fname)
    ale_plot, = plot(rkm, eta, 'b', label='ALE3D')
else:
    ale_plot = None


legend(loc='upper right', framealpha=1)
xlabel('radial distance (km)')
ylabel('surface elevation (m)')
xlim(0,rmax/1e3)
ylim(-500,500)
grid(True)
title_text = title(title_base)


def update(t):
    """
    Update an existing plot by resetting the data.
    Use the frames from each simulation closest to the given time t.
    """
    print(f'+++ update called at t = {t}')

    if include_ale:
        # ALE:
        frameno_ale = find_frame_ale(t)
        fname = '%s_surface_%s.mat' % (fname_prefix_ale, str(frameno_ale).zfill(4))
        rkm, eta, t_ale = C.load_surf_ale(datadir_ale, fname)
        ale_plot.set_data(rkm, eta)

    # mfluid:
    frameno_mfluid = find_frame_mfluid(t)
    eta = load_surface(frameno_mfluid, rvals, zvals, outdir_mfluid)
    print(f'+++ len(rvals)={len(rvals)}, len(eta)={len(eta)}')
    mfluid_plot.set_data(rkm_mfluid, eta)
    #import pdb; pdb.set_trace()

    if include_ale:
        if max(abs(t_ale-t),abs(t_mfluid-t)) > 0.1:
            # t not found exactly in one sim or the other:
            print('MISMATCH t = %.1f: frameno_ale = %i at t = %1.f, frameno_mfluid = %i at t = %1.f' \
                % (t,frameno_ale, t_ale, frameno_mfluid, t_mfluid))

        title_text.set_text(f'{title_base} at t = {t:6.1f}' \
                + f'  (ALE3D at t = {t_ale:6.1f}')
    else:
        title_text.set_text(f'{title_base} at t = {t:6.1f}')

if __name__ == '__main__':

    print('Making anim...')
    #times = tf_mfluid[:,1]  # all times found for mfluid solution
    times = [0, 10]  # subset of times for testing
    anim = animation.FuncAnimation(fig, update, frames=times,
                                   interval=200, blit=False)

    # Output files:

    fname_mp4 = movie_name + '.mp4'

    fname_html = None
    #fname_html = movie_name + '.html'

    if fname_mp4:
        fps = 5
        print('Making mp4...')
        writer = animation.writers['ffmpeg'](fps=fps)
        anim.save(fname_mp4, writer=writer)
        print("Created %s" % fname_mp4)

    if fname_html:
        # html version:
        animation_tools.make_html(anim, file_name=fname_html, title=name)
