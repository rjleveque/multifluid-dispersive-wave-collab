from pylab import *
from clawpack.pyclaw import Solution
from clawpack.visclaw import gridtools
import glob
ion()


# Will need to adjust rmax, dr, zmin, zmax for a larger crater...
# (r and z refer to x,y from 2D mfluid simulation)

# define radial values r on which to compute surface:
rmax = 2e3
#rmax = 3e3
#rmax = 20e3
#dr = 10.
#rvals = arange(dr,rmax,dr)

# vertical z values on which to evaluate 2D solution for computing surface:
zmin = -900. # any depth in domain that's always below surface
zmax = 1000.  # any elevation in domain that's always above surface
dz = 0.01  # needs to be pretty small
zvals = arange(zmin, zmax, dz)

def load_surface2(frameno, rvals, outdir='_output', file_format='binary'):

    # read in the AMR solution for this time frame:
    framesoln = Solution(frameno, path=outdir, file_format=file_format)

    eta = zeros(len(rvals))  # initialize array

    for i,r in enumerate(rvals):
        ri = r*ones(len(zvals))
        # extract zfa = q[5,:,:] on vertical transect at r=ri:
        # (taking data from finest level that exists at each point in zvals)
        zfa = gridtools.grid_output_2d(framesoln, 5, ri, zvals, method='linear')
        have = (1. - zfa.sum()) * dz
        eta[i] = have - zmin
    return eta



def load_surface(frameno, rvals,  outdir='_output', file_format='binary'):
    """
    read in the mfclaw AMR solution for this time frame, interpolating
    to the radial coordinate r values in rvals
    """
    from clawpack.pyclaw import Solution
    from clawpack.visclaw import gridtools

    framesoln = Solution(frameno, path=outdir, file_format=file_format)

    #print(f'Loaded solution at time {framesoln.t:.3f}, computing eta...')

    eta = zeros(len(rvals))  # initialize array

    for i,r in enumerate(rvals):
        ri = r*ones(len(zvals))
        # extract zfa = q[5,:,:] on vertical transect at r=ri:
        # (taking data from finest level that exists at each point in zvals)
        zfa = gridtools.grid_output_2d(framesoln, 5, ri, zvals, method='linear')
        have = (1. - zfa).sum() * dz
        eta[i] = have + zmin
    return eta


if __name__=='__main__':

    # loop over frames as specified by framenos below

    outdir1 = 'm4to6k300mcrater_250by125_3lev44/_output'
    outdir2 = 'm4to6k300mcrater_250by125_3lev64/_output'
    outdir3 = 'm4to6k300mcrater_500by250_3lev44/_output'
    outdir4 = 'm4to6k300mcrater_250by125_3lev66/_output'

    #outdir1 = '_outputa'
    #outdir2 = '_outputb'
    #outdir3 = '_outputc'
    #outdir4 = 'm4to6k300mcrater_250by125_3lev66/_output'

    # list of outdirs to iterate over in each frame plot:
    outdirs = [outdir1, outdir2,outdir3,outdir4]
    #labels = ['6k atm','7k atm','8k atm','finer 6k']  # labels to appear in legend
    #drs = [6.6667,6.6667,6.6667, 2.2222]     
    labels = ['h=5.0','h=3.333','h=2.5','h=2.222']  # labels to appear in legend
    drs = [5.0,3.3333333,2.5,2.22222]
    colors = ['b','r','g','m','k']  # color to plot for each outdir

    # time t for each frame in outdir1  (Assume they are the same in outdir2)
    t_outdir1 = {}
    fortq_files = glob.glob(f'{outdir1}/fort.q*')
    fortq_files.sort()
    for frameno,qfile in enumerate(fortq_files):
        tfile = qfile.replace('q','t')
        firstline = open(tfile).readline()
        t_outdir1[frameno] = float(firstline.split()[0])

    # all frames found in outdir1
    framenos = array(range(len(fortq_files)))
    #framenos = array(range(4)) # test

    figure(200, figsize=(10,6))
    show()

    ans = ''
    frameno = 0
    while ans != 'q':
        clf()
        for outdir,label,color,dr in zip(outdirs,labels,colors,drs):
            rvals = arange(dr/2,rmax,dr)
            if 1: 
                eta1 = load_surface(frameno, rvals, outdir)
                plot(rvals, eta1, color, label=label)
            else:  #except:
               print(f"frame not found in outdir={outdir}")
               pass
        legend(loc='upper right', framealpha=1)
        grid(True)
        xlim(0,rmax)
        ylim(zmin,zmax)
        title(f'Frame {frameno} at time t = {t_outdir1[frameno]:.2f} seconds')
        draw()
        #fname = 'compHeightAtmos300m%s.png'%str(frameno).zfill(4)
        fname = 'compRes300m%s.png'%str(frameno).zfill(4)
        savefig(fname)
        print("created ",fname)

        #if 1:
        if 0:
            # interactive prompt:
            ans = input('Return for next frame, integer to jump, or q to quit... ')
            try:
                frameno = int(ans)
                print(f'Jumping to frameno = {frameno}')
            except:
                frameno += 1
        else:
            # always advance to next frame to make plots without prompting
            frameno += 1

        if frameno > framenos.max():
            print(f'max frame reached framenos.max() = {framenos.max()}')
            ans = 'q'

