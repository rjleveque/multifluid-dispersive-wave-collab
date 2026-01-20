"""
mfluid_to_shelf with mfrk2Dev AMR output
"""

from pylab import *
import os,sys,glob
from scipy.interpolate import interp1d

root_dir = '/Users/rjl/git/Forks/multifluid-dispersive-wave-collab/crater_code/5eqns_transport/'

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

mfclaw_tools = fullpath_import(f'{root_dir}/mfrk2Dev/mfclaw_tools.py')

C = fullpath_import(os.path.join(root_dir, 'compareALE3D.py'))

AGT = os.environ['AGT']
LW = fullpath_import(os.path.join(AGT, 'linear_transforms', 'linear_waves.py'))

# directory for mfluid simulation data used to initialize:
outdir_mfluid = 'RC300/_output'

h0 = 4000.
RC = 300.
tAiry_start = 60.   # initial eta(r,t) from mfluid solution at this time
tAiry_end = 3600.    # compute up to this time
rAiry_end = 100e3    # capture the solution eta(r,t) as fcn of t at this r

# domain for computing radial Hankel solution:
xlower = 0
xupper = 120e3  # at least as large as rAiry_end, how much larger??

# directory for saving boundary conditions time series:
savefile_dir = './BC'
os.system('mkdir -p %s' % savefile_dir)

# file for saving time series of eta and u (for SWE and Bouss) at rAiry_end:
# (will be saved in directory savefile_dir)
savefile_name = 'eta_hu_bc_RC%i_h%s_tstart%i_rbc%ikm' \
    % (int(RC), int(h0), int(tAiry_start), int(rAiry_end/1000))


grav = 9.81

print('outdir_mfluid = ',outdir_mfluid)
#tf_mfluid, find_frame_mfluid = C.load_times_mfluid(outdir_mfluid)


tf_mfluid, find_frame_mfluid = mfclaw_tools.load_times_mfclaw(outdir_mfluid)


def save_eta_u(t, r, eta, u, fname):
    """
    save frame of solution at one time t for all r
    """
    d = vstack((r,eta,u)).T
    savetxt(fname, d, header='%.1f\n%i' % (t,len(eta)), comments='')
    print('Created ',fname)


# Initial eta and u
L = xupper - xlower
mx = 5000  # Is this large enough?  How to choose?
mk = 4000  # Is this large enough?  How to choose?

dx = L/mx
r = linspace(dx/2, xupper-dx/2, mx)

kmax = 1.5 * 2*pi/RC
k = linspace(1e-6,kmax,mk)
#k = linspace(1e-6,0.005,mk)

def make_Airy_data(tAiry_start, r=r, k=k):
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
    """
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


def evolve_Airy(etahat,tAiry_start,tAiry_end,k=k,r=r):
    """
    Computes Airy solution eta(r,elapsed_time), given transform at initial time.
    Not used now, instead we make time series at rAiry_end for use as
    boundary condition.
    """
    omega = lambda k: LW.omega_airy(k,h0)
    elapsed_time = tAiry_end - tAiry_start
    eta, u = LW.eta_u_radial(elapsed_time,r,k,etahat,omega,h0)
    return eta, u

# ===============
# main function:

def make_bc(rAiry_end, tAiry_end, tAiry_start, etahat_start, omega=None,
            save_txt=False):

    """
    Create boundary data eta(rAiry_end,t) at fixed radius rAiry_end for
    a series of times from t=tAiry_start to tAiry_end
    (currently time increment dt is hard-wired).

    Also compute the corresponding hu(rAiry_end,t) for SWE or Bouss
    (by using each of those dispersion relations when estimating hu from eta)
    """
    from scipy.signal import find_peaks
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
    else:
        raise NotImplementedError('have not implemented cgroup approach')
        # approximate group velocity of waves arriving at time t:
        #cgroup = (rAiry_end - 1) / t
        # need to invert for k and then compute cphase = omega(k)/k
        #hu = cphase * eta

    # save the time series for eta and both momenta hu_swe and hu_sgn
    # so that this provides BCs at rAiry_end for either SWE or SGN on shelf:
    d = vstack((t,eta,hu_swe,hu_sgn)).T
    if save_txt:
        #fname = 'eta_hu_bc_%skm.txt' % int(rAiry_end/1e3)
        fname = os.path.join(savefile_dir, savefile_name + '.txt')
        savetxt(fname, d, header='%.0f\n%.1f\n%i' \
                                 % (rAiry_end,tAiry_start,len(eta)),
                comments='',fmt='%20.10e')
        print('Created ',fname)
    return t,eta,hu_swe,hu_sgn


def plot_bc(t,eta,hu_swe,hu_sgn,rAiry_end,tAiry_end,
            save_png=False):
    """
    Plot the time series of eta at rAiry_end in upper subplot,
    and both hu_swe and hu_sgn in lower subplot,
    for the times in t (converted from seconds to minutes)
    and with xlim determined by tAiry_end.
    """
    figure(figsize=(10,7))
    subplot(211)
    plot(t/60,eta,'b')
    title('time series of eta, hu at rAiry_end = %s km' \
            % int(rAiry_end/1e3), fontsize=15)
    grid(True)
    xlim(0,tAiry_end/60)
    #xlabel('time (minutes)', fontsize=13)
    ylabel('surface eta (m)', fontsize=13)
    xticks(fontsize=13)
    yticks(fontsize=13)

    subplot(212)
    plot(t/60,hu_swe,'b',label='hu for swe')
    plot(t/60,hu_sgn,'r',label='hu for sgn')
    #title('time series of hu(r2,t) at r2 = rAiry_end = %s km' \
    #        % int(rAiry_end/1e3), fontsize=15)
    grid(True)
    legend(loc='upper left', framealpha=1, fontsize=11)
    xlim(0,tAiry_end/60)
    xlabel('time (minutes)', fontsize=13)
    ylabel('momentum hu (m**2/s)', fontsize=13)
    xticks(fontsize=13)
    yticks(fontsize=13)

    tight_layout()

    if save_png:
        #fname = 'eta_hu_bc_%skm.png' % int(rAiry_end/1e3)
        fname = os.path.join(savefile_dir, savefile_name + '.png')
        savefig(fname, bbox_inches='tight')
        print('Created ',fname)

def test_bc(tAiry_start=250, etahat_start=None):

    if etahat_start is None:
        eta_start, etahat_start, u_start = make_Airy_data(tAiry_start=tAiry_start,
                                                          r=r, k=k)
    #rAiry_end = 200.e3
    #tAiry_end = 3.*3600.
    rAiry_end = 100.e3
    tAiry_end = 4000.

    t,eta,hu_swe,hu_sgn = make_bc(rAiry_end, tAiry_end, tAiry_start, etahat_start,
                       save_txt=True)
    plot_bc(t, eta,hu_swe,hu_sgn, rAiry_end, tAiry_end, save_png=True)

    return t,eta,hu_swe,hu_sgn

def compare_airy_sgn(etahat_start=None):
    """
    Make the BC's at rAiry_end using Airy as usual
    Also make them using the SGN dispersion relation instead, to see
    how much difference there would be if SGN were used to propagate from
    mfluid regime to shelf.
    """
    tAiry_start = 250
    if etahat_start is None:
        eta_start, etahat_start, u_start = make_Airy_data(tAiry_start=tAiry_start,
                                                          r=r, k=k)
    rAiry_end = 200.e3
    tAiry_end = 3*3600.
    print('tAiry_end = %.1f sec' % tAiry_end)

    # make_bc using the SGN dispersion relation:
    omega = lambda k: LW.omega_sgn(k,h0,1.153)
    t,eta,hu = make_bc(rAiry_end, tAiry_end, tAiry_start,
                       etahat_start,omega=omega, save_txt=False)
    plot_bc(t, eta, rAiry_end, tAiry_end, plot_eta=False)
    plot(t/60,eta,'r',label='SGN')

    # make_bc using the Airy dispersion relation:
    t,eta,hu = make_bc(rAiry_end, tAiry_end, tAiry_start, etahat_start,
                       save_txt=False)
    #plot_bc(t, eta, rAiry_end, tAiry_end, plot_eta=False)
    plot(t/60,eta,'b',label='Airy')
    legend(loc='upper right', framealpha=1, fontsize=12)

    return t,eta,hu

def plot_etahat(k, etahat, t_etahat, save_png=False):
    figure(11,figsize=(10,6))
    plot(k, etahat, 'r')
    title(f'Hankel transform etahat(k) at time {t_etahat/60:.2f} minutes')
    xlabel('wavenumber k')
    ylabel('etahat(k)')
    grid(True)
    if save_png:
        fname = f'Htransform_t{t_etahat:.0f}.png'
        fname = savefile_name.replace('eta_hu_bc','Htransform')
        fname = os.path.join(savefile_dir, fname)
        savefig(fname, bbox_inches='tight')
        print('Created ',fname)

if __name__ == '__main__':

    eta_start, etahat_start = make_Airy_data(tAiry_start=tAiry_start, r=r, k=k)

    plot_etahat(k, etahat_start, tAiry_start, save_png=True)

    t,eta,hu_swe,hu_sgn = make_bc(rAiry_end, tAiry_end, tAiry_start,
                                  etahat_start, save_txt=True)

    plot_bc(t, eta,hu_swe,hu_sgn, rAiry_end, tAiry_end, save_png=True)
