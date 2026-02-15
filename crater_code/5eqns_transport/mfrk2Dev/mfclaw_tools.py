"""
tools for working with mfclaw AMR output (mfrk2Dev)
"""

from pylab import *
import os,sys
from scipy.interpolate import interp1d
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


# vertical z values on which to evaluate 2D solution for computing surface:
zmin = -900. # any depth in domain that's always below surface
zmax = 1000.  # any elevation in domain that's always above surface
dz = 1  #0.01  # needs to be pretty small
zvals = arange(zmin, zmax, dz)

grav = 9.81

def load_surface_nc(fname_nc):
    """
    Load the surface as a function of r at a sequence of times t from
    a netcdf file previously created by `extract_surface.py`
    The netcdf file should contain an xarray DataArray with coordinates r,t.
    """
    import xarray as xr

    surface_data = xr.open_dataarray(fname_nc)
    tf = vstack((array(range(len(surface_data.coords['t']))),
                 surface_data.coords['t'])).T

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

    return tf, find_frame, surface_data

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
        k2 = where(zfa_coarse < 1-tol)[0].max() + 1
        if k2 <= k1:
            k1 = k1 - 1
            k2 = k2 + 1
        k1 = max(k1,0)
        k2 = min(k2,len(zvals_coarse)-1)

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

def make_Airy_data(tAiry_start, r, k):
    """
    Load the mfluid solution at time tAiry_start
    (or from frame that is closest to specified time and adjust tAiry_start)
    and then:
      - smooth near the origin if needed (also for large r?)
      - resample at uniform times specified by the `r` array (linear interp)
      - compute etahat as Hankel transform of eta, using wavenumbers specified
        by `k` array
    Returns eta, etahat, u
        (where u was also read from mfluid solution, assuming
        depth-averaged values of u were also in fort.c* files, but is not
        currently being used in computing BCs)

    The etahat function returned is can then be passed into make_bc as
    etahat_start to compute BCs at rAiry_end based on this solution.

    Adapted from mfclaw_to_shelf.py
    """
    LW = fullpath_import(f'{AGT}/linear_transforms/linear_waves.py')

    # load mfluid solution at time tAiry_start:
    # (or frame that is closest to specified time)
    frameno_mfluid = find_frame_mfluid(tAiry_start)

    #rkm_mfluid, eta_mfluid, u_mfluid, t_mfluid = \
    #                C.load_mfluid(outdir_mfluid, frameno_mfluid, h0)

    r_mfluid, eta_mfluid = mfclaw_tools.load_surface(frameno_mfluid,
                                                     outdir_mfluid)


    if 0:
        if t_mfluid != tAiry_start:
            print('Resetting tAiry_start from %s to %s' % (tAiry_start,t_mfluid))
            tAiry_start = t_mfluid
        r_mfluid = rkm_mfluid * 1000.  # convert from km to m

    if 0:
        # damp out near origin:
        eta_mfluid = where(r_mfluid<3e3,
                            eta_mfluid*exp(-((r_mfluid-3e3)/1e3)**2), eta_mfluid)
        u_mfluid = where(r_mfluid<3e3,
                            u_mfluid*exp(-((r_mfluid-3e3)/1e3)**2), u_mfluid)


    # resample eta_mfluid at r, possibly extending by 0 to larger distance:
    print('resampling eta_mfluid from r_mfluid (out to %s m with dr=%s)' \
            % (r_mfluid[-1], r_mfluid[-1]-r_mfluid[-2]) \
            + '\n                               to r (out to %s m with dr=%s)' \
            % (r[-1],dx))
    etafcn = interp1d(r_mfluid, eta_mfluid, fill_value=0., bounds_error=False)
    eta = etafcn(r)
    #ufcn = interp1d(r_mfluid, u_mfluid, fill_value=0., bounds_error=False)
    #u = ufcn(r)  # not needed?

    # Also need to damp out for large r since mfluid soln may not fully decay?

    print('Computing Hankel transform...')
    etahat = LW.Htransform(r,eta,k)

    return eta, etahat

def make_bc(rAiry_end, tAiry_end, tAiry_start, etahat_start, k, h0,
            t=None, omega=None, savefile=None, return_hu_airy=False):

    """
    Create boundary data eta(rAiry_end,t) at fixed radius rAiry_end for
    a series of times from t=tAiry_start to tAiry_end
    (currently time increment dt is hard-wired).

    Also compute the corresponding hu(rAiry_end,t) for SWE or Bouss
    (by using each of those dispersion relations when estimating hu from eta)
    """
    from scipy.signal import find_peaks
    LW = fullpath_import(f'{AGT}/linear_transforms/linear_waves.py')

    if t is None:
        dt = 1.
        t = arange(tAiry_start, tAiry_end+dt/2, dt)

    print(f'Computing time series at time {tAiry_end:.1f}...')
    elapsed_time = t - tAiry_start

    if omega is None:
        omega = lambda k: LW.omega_airy(k,h0)

    eta = LW.eta_tseries_radial(elapsed_time,rAiry_end,k,etahat_start,omega,h0)

    if 1:
        # for use with SWE on shelf: (using right-going eigen-vector)
        cphase = sqrt(grav*h0) # independent of k
        hu_swe = sqrt(grav*h0) * eta  # for SWE

        # for use with SGN on shelf:

        # currently estimating period of wave at rAiry_end
        # and frequency omega around each time t in time series,
        # and from this estimating wave number k in space by inverting
        # dispersion relation omega(k) for the Bouss solver that will be
        # used on the shelf (assumed to be SGN with alpha = 1.153)
        # so that this can be used as BC in that 1D geoclaw simulation

        # Using the fact that hu = cphase * eta in general for any linear
        # dispersion relation, where cphase is the phase speed omega(k)/k.

        omega = lambda k: LW.omega_sgn(k,h0,alpha=1.153)
        jpeaks = find_peaks(eta)[0]
        jpeaks = jpeaks[4:]  # throw away first few points  - ADJUST?
        tpeaks = t[jpeaks]
        Tperiod = diff(t[jpeaks])
        tpeaks = hstack((0, tpeaks, tpeaks[-1]+1000))
        #Tperiod = hstack((240, Tperiod, Tperiod[-1], Tperiod[-1]))
        Tperiod = hstack((Tperiod[0], Tperiod, Tperiod[-1], Tperiod[-1]))
        Tfcn = interp1d(tpeaks, Tperiod) #, fill_value=0., bounds_error=False)
        Tt = Tfcn(t) # estimate of period at each time in t
        kk = linspace(0,0.01,1000)
        ww = omega(kk)
        kfcn = interp1d(ww,kk)  # k as a function of omega
        kt = kfcn(2*pi/Tt)  # estimate of k at each t
        cphase = omega(kt)/kt
        hu_sgn = cphase * eta

        if return_hu_airy:
            # also compute for Airy
            kk = linspace(0,0.015,1500)
            omega = lambda k: LW.omega_airy(k,h0)
            ww = omega(kk)
            kfcn = interp1d(ww,kk)  # k as a function of omega
            kt = kfcn(2*pi/Tt)  # estimate of k at each t
            cphase = omega(kt)/kt
            hu_airy = cphase * eta

    else:
        raise NotImplementedError('have not implemented cgroup approach')
        # approximate group velocity of waves arriving at time t:
        #cgroup = (rAiry_end - 1) / t
        # need to invert for k and then compute cphase = omega(k)/k
        #hu = cphase * eta

    # save the time series for eta and both momenta hu_swe and hu_sgn
    # so that this provides BCs at rAiry_end for either SWE or SGN on shelf:
    d = vstack((t,eta,hu_swe,hu_sgn)).T
    if savefile is not None:
        #fname = 'eta_hu_bc_%skm.txt' % int(rAiry_end/1e3)
        #fname = os.path.join(savefile_dir, savefile_name + '.txt')
        savetxt(savefile, d, header='%.0f\n%.1f\n%i' \
                                 % (rAiry_end,tAiry_start,len(eta)),
                comments='',fmt='%20.10e')
        print('Created ',savefile)

    if return_hu_airy:
        return t, eta, hu_swe, hu_sgn, hu_airy
    else:
        return t, eta, hu_swe, hu_sgn


def plot_bc(t,eta,hu_swe,hu_sgn,rAiry_end,tAiry_end,
            tAiry_start,RC,h0,savefile=None):
    """
    Plot the time series of eta at rAiry_end in upper subplot,
    and both hu_swe and hu_sgn in lower subplot,
    for the times in t
    and with xlim determined by tAiry_end.
    """
    figure(figsize=(12,8))
    subplot(211)
    plot(t,eta,'b')
    title(f'time series of eta, hu at rAiry_end = {int(rAiry_end/1e3)} km\n' \
            + f'RC={RC:.0f}, h0={h0:.0f}, tAiry_start={tAiry_start:.0f}',
          fontsize=15)
    grid(True)
    xlim(0,tAiry_end)
    ylabel('surface eta (m)', fontsize=13)
    xticks(fontsize=13)
    yticks(fontsize=13)

    subplot(212)
    plot(t,hu_swe,'b',label='hu for swe')
    plot(t,hu_sgn,'r',label='hu for sgn')
    #title('time series of hu(r2,t) at r2 = rAiry_end = %s km' \
    #        % int(rAiry_end/1e3), fontsize=15)
    grid(True)
    legend(loc='upper left', framealpha=1, fontsize=11)
    xlim(0,tAiry_end)
    xlabel('time (seconds)', fontsize=13)
    ylabel('momentum hu (m**2/s)', fontsize=13)
    xticks(fontsize=13)
    yticks(fontsize=13)

    tight_layout()

    if savefile is not None:
        savefig(savefile, bbox_inches='tight')
        print('Created ',savefile)
