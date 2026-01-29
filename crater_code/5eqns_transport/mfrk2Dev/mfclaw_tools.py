"""
tools for working with mfclaw AMR output (mfrk2Dev)
"""

from pylab import *
import os,sys
#from scipy.interpolate import interp1d
import glob
from clawpack.pyclaw import Solution
from clawpack.visclaw import gridtools
from clawpack.clawutil.data import ClawData


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


AGT = os.environ['AGT']
#LW = fullpath_import(os.path.join(AGT, 'linear_transforms', 'linear_waves.py'))


# vertical z values on which to evaluate 2D solution for computing surface:
zmin = -900. # any depth in domain that's always below surface
zmax = 1000.  # any elevation in domain that's always above surface
dz = 1  #0.01  # needs to be pretty small
zvals = arange(zmin, zmax, dz)

grav = 9.81


def load_times_mfclaw(outdir):
    """
    load the times from the fort.t files and create a function find_frame
    that finds the frame closest to a given time.

    :Outputs:
    - tf: array of times
    - fund_frame: function of t
    """
    #print('+++ outdir = ',outdir)
    fortt_files = glob.glob('%s/fort.t*' % outdir)
    fortt_files.sort()
    fortt_files = [f for f in fortt_files if 'tck' not in f]
    #print(fortt_files)
    if len(fortt_files) == 0:
        raise IOError('No fort.t files found in ',outdir)
    times = []
    for f in fortt_files:
        frameno = int(f[-4:])
        t = float(open(f).readline().split()[0])
        times.append((frameno,float(round(t,1))))
    tf = array(times)

    def find_frame(time):
        """
        find frameno for best match to time
        """
        if time < tf[0,1]:
            k = 0
        else:
            k = where(tf[:,1] < time+1e-6)[0].max()
            #print('+++ k = %i' % k, '   tf[k:k+2, :] = ',tf[k:k+2, :])
            if (k < tf.shape[0]-1):
                if (tf[k+1,1] - time) < (time - tf[k,1]):
                    k = k+1
        return int(tf[k,0]), tf[k,1]  # also return time of this frame

    print('Found %i mfclaw frames up to time %.1f seconds' \
                % (tf.shape[0], tf[:,1].max()))

    return tf, find_frame



def load_surface(frameno, outdir='_output', file_format='binary', rmax=inf):
    """
    Read in the mfclaw AMR solution for this time frame, interpolating
    to the radial coordinate r values determined by the finest grid resolution
    in this run (based on level 1 dx and refinement ratios).

    rmax is optional upper bound on r of interest,
         if inf, use clawdata.upper[0], the right boundary of the domain
         and in any case replace rmax by min(rmax, clawdata.upper[0])
    """

    clawdata = ClawData()
    clawdata.read(f'{outdir}/claw.data',force=True)
    amrdata = ClawData()
    amrdata.read(f'{outdir}/amr.data',force=True)

    h0 = -clawdata.lower[1]     # ocean depth
    zupper = clawdata.upper[1]  # top of atmosphere
    print(f'load_surface called with rmax = {rmax}')
    rmax = min(rmax, clawdata.upper[0])    # don't go past right boundary

    # determine finest AMR resolution dr_fine:
    num_cells_r_fine = clawdata.num_cells[0] * \
              prod(amrdata.refinement_ratios_x[:(amrdata.amr_levels_max-1)])
    dr_fine = clawdata.upper[0] / num_cells_r_fine
    rvals = arange(0.5*dr_fine, rmax, dr_fine)  # grid at this resolution dr_fine
    print(f'For frameno={frameno} from {outdir}:')
    print(f'    rmax = {rmax}, dr_fine = {dr_fine}, len(rvals) = {len(rvals)}')

    # load AMR solution:
    framesoln = Solution(frameno, path=outdir, file_format=file_format)

    print(f'Loaded solution at time {framesoln.t:.3f}, computing eta...')

    eta = zeros(len(rvals))  # initialize array
    eta_lowerbound = zeros(len(rvals))
    eta_upperbound = zeros(len(rvals))

    dz_coarse = 10.
    zvals_coarse = arange(-h0+dz_coarse/2, zupper, dz_coarse)

    for i,r in enumerate(rvals):
        ri_coarse = r*ones(len(zvals_coarse))

        # first find indicator zfa = q[5,:,:] on coarse grid in z:
        zfa_coarse = gridtools.grid_output_2d(framesoln, 5, ri_coarse,
                                              zvals_coarse, method='linear')

        # determine range of z over which zfa is varying between 0 and 1:
        tol = 1e-5
        k1 = where(zfa_coarse > tol)[0].min()
        k2 = where(zfa_coarse < 1-tol)[0].max()
        if k2 <= k1:
            k1 = k1 - 1
            k2 = k2 + 1

        #print(f'+++ k1={k1}, k2={k2}')

        # now evaluate zfa on fine grid in z between these limits:
        dz_fine = 0.01
        zvals_fine = arange(zvals_coarse[k1], zvals_coarse[k2], dz_fine)
        ri_fine = r*ones(len(zvals_fine))
        #print(f'zvals_fine has length {len(zvals_fine)} from z={zvals_fine[0]:.1f} to z = {zvals_fine[-1]:.1f}')

        zfa = gridtools.grid_output_2d(framesoln, 5, ri_fine,
                                       zvals_fine, method='linear')

        # integrate the indicator function that is 0 in water, 1 in air
        # to estimate location of interface:
        have = (1. - zfa).sum() * dz_fine
        eta[i] = have + zvals_fine[0]

        eta_lowerbound[i] = zvals_fine[0]
        eta_upperbound[i] = zvals_fine[-1]


    if 0:
        # for debugging return lower and upper limits of zvals_fine too:
        return rvals, eta, eta_lowerbound, eta_upperbound
    else:
        return rvals, eta


def load_surface_fromtop(frameno, outdir='_output', file_format='binary'):
    """
    Read in the mfclaw AMR solution for this time frame, interpolating
    to the radial coordinate r values determined by the finest grid resolution
    in this run (based on level 1 dx and refinement ratios).

    Version that sums from top down rather than bottom up, maybe different
    results if there are bubbles near the axis?
    """

    try:
        clawdata = ClawData()
        clawdata.read(f'{outdir}/claw.data',force=True)
        amrdata = ClawData()
        amrdata.read(f'{outdir}/amr.data',force=True)

        h0 = -clawdata.lower[1]     # ocean depth
        zupper = clawdata.upper[1]  # top of atmosphere
        rmax = clawdata.upper[0]    # right boundary

        # determine finest AMR resolution dr_fine:
        num_cells_r_fine = clawdata.num_cells[0] * \
                  prod(amrdata.refinement_ratios_x[:(amrdata.amr_levels_max-1)])
        dr_fine = rmax / num_cells_r_fine
    except:
        print('*** could not load from claw.data, amr.data')
        dr_fine = 50.
        h0 = -4000.
        zupper = 6000.
        rmax = 6000.

    rvals = arange(0.5*dr_fine, rmax, dr_fine)  # grid at this resolution dr_fine
    print(f'rmax = {rmax}, dr_fine = {dr_fine}, len(rvals) = {len(rvals)}')

    # load AMR solution:
    framesoln = Solution(frameno, path=outdir, file_format=file_format)

    print(f'Loaded solution at time {framesoln.t:.3f}, computing eta...')

    eta = zeros(len(rvals))  # initialize array
    eta_lowerbound = zeros(len(rvals))
    eta_upperbound = zeros(len(rvals))

    dz_coarse = 10.
    zvals_coarse = arange(-h0+dz_coarse/2, zupper, dz_coarse)

    for i,r in enumerate(rvals):
        ri_coarse = r*ones(len(zvals_coarse))

        # first find indicator zfa = q[5,:,:] on coarse grid in z:
        zfa_coarse = gridtools.grid_output_2d(framesoln, 5, ri_coarse,
                                              zvals_coarse, method='linear')

        # determine range of z over which zfa is varying between 0 and 1:
        tol = 1e-5
        k1 = where(zfa_coarse <= tol)[0].max()  # topmost pure water cell
        k2 = where(zfa_coarse < 1-tol)[0].max() # topmost pure air cell
        if k2 <= k1:
            k1 = k1 - 1
            k2 = k2 + 1

        #print(f'+++ k1={k1}, k2={k2}')

        # now evaluate zfa on fine grid in z between these limits:
        dz_fine = 0.01
        zvals_fine = arange(zvals_coarse[k1], zvals_coarse[k2], dz_fine)
        ri_fine = r*ones(len(zvals_fine))
        #print(f'zvals_fine has length {len(zvals_fine)} from z={zvals_fine[0]:.1f} to z = {zvals_fine[-1]:.1f}')

        zfa = gridtools.grid_output_2d(framesoln, 5, ri_fine,
                                       zvals_fine, method='linear')

        # integrate the indicator function that is 0 in water, 1 in air
        # to estimate location of interface:

        # from bottom:
        #have = (1. - zfa).sum() * dz_fine
        #eta[i] = have + zvals_fine[0]

        # from top:
        have = zfa.sum() * dz_fine
        eta[i] = zvals_fine[-1] - have

        eta_lowerbound[i] = zvals_fine[0]
        eta_upperbound[i] = zvals_fine[-1]

    if 0:
        # for debugging return lower and upper limits of zvals_fine too:
        return rvals, eta, eta_lowerbound, eta_upperbound
    else:
        return rvals, eta


def plotzfa(frameno, r, outdir='_output', file_format='binary'):
    """
    For debugging, plot zfa indicator as function of z at specified r
    """

    figure(2)
    clf()
    framesoln = Solution(frameno, path=outdir, file_format=file_format)

    clawdata = ClawData()
    clawdata.read(f'{outdir}/claw.data',force=True)
    amrdata = ClawData()
    amrdata.read(f'{outdir}/amr.data',force=True)

    h0 = -clawdata.lower[1]     # ocean depth
    zupper = clawdata.upper[1]  # top of atmosphere
    rmax = clawdata.upper[0]    # right boundary

    dz_coarse = 10.
    zvals_coarse = arange(-h0+dz_coarse/2, zupper, dz_coarse)

    ri_coarse = r*ones(len(zvals_coarse))

    # first find indicator zfa = q[5,:,:] on coarse grid in z:
    zfa_coarse = gridtools.grid_output_2d(framesoln, 5, ri_coarse,
                                          zvals_coarse, method='linear')

    have = (1. - zfa_coarse).sum() * dz_coarse
    eta_coarse = have + zvals_coarse[0]
    print(f'Coarse: eta = {eta_coarse:.3f}')
    #plot(zvals_coarse, zfa_coarse, 'r-')
    grid(True)
    #plot([eta_coarse,eta_coarse],[0,1],'r--')


    # determine range of z over which zfa is varying between 0 and 1:
    tol = 1e-5
    k1 = where(zfa_coarse > tol)[0].min()
    k2 = where(zfa_coarse < 1-tol)[0].max()
    if k2 <= k1:
        k1 = k1 - 1
        k2 = k2 + 1

    #print(f'+++ k1={k1}, k2={k2}')


    # now evaluate zfa on fine grid in z between these limits:
    dz_fine = 0.01
    zvals_fine = arange(zvals_coarse[k1], zvals_coarse[k2], dz_fine)
    ri_fine = r*ones(len(zvals_fine))
    #print(f'zvals_fine has length {len(zvals_fine)} from z={zvals_fine[0]:.1f} to z = {zvals_fine[-1]:.1f}')

    zfa_fine = gridtools.grid_output_2d(framesoln, 5, ri_fine,
                                   zvals_fine, method='linear')

    # integrate the indicator function that is 0 in water, 1 in air
    # to estimate location of interface:
    have = (1. - zfa_fine).sum() * dz_fine
    eta_fine = have + zvals_fine[0]

    print(f'Fine: eta = {eta_fine:.3f}')
    print(f'zvals_fine has length {len(zvals_fine)} from z={zvals_fine[0]:.1f} to z = {zvals_fine[-1]:.1f}')

    figure(2)
    plot(zvals_fine, 1-zfa_fine, 'b')

    plot([eta_fine,eta_fine],[0,1],'b--')
    ylim(-0.05, 1.05)

    return zvals_fine, zfa_fine
