"""
Code to perform multiple mfclaw runs, with a different parameters
for each, and then make plots for each run.

This version varies the ocean depth h0.

If run_code is False the code will not be run (e.g. if the runs have been
done and only plotting should be done).

If make_plots is False the plotting will not be done.

It is assumed all parameters set in setrun are the same for all realizations,
except for the depth h0 to use.  The setrun_case.py module contains a
function setrun that takes case as a parameter so that case['h0']
can be used in setrun.

For making plots, it is assumed a setplot_case.py module contains a
function setplot that takes case as a parameter so that the event name
can be extracted for adding to plots / filenames.

The function make_all_cases_h0 returns *caselist*, a list of
dictionaries.  Each dictionary should define whatever parameters
are needed for one case.

The function run_one_case_clawpack(case) takes a dictionary *case* as input
and does what's necessary to run mfclaw for that one case.

Before running, make sure the following are properly set:

    h0vals, the list of ocean depths

    runs_dir is set to the destination for _output and _plots directories
        (possibly on a scratch disk).

    nprocs is set to the number of mfclaw runs that can be run in parallel.
        (each run will also use OpenMP with OMP_NUM_THREADS threads)

Set dry_run = False before executing to actually run mfclaw.

Note that the call to multip_tools.run_many_cases_pool should be in __main__

This script can be executed from the command line as:
    python runclaw_makeplots_h0vals.py NPROCS
where NPROCS is the number of jobs to run in parallel using the
Python mulitprocessing tools.

"""

from numpy import *
import os,sys,glob

from clawpack.clawutil.util import fullpath_import
from clawpack.clawutil import multip_tools, clawmultip_tools
import make_cases


#dry_run = True  # If True, only print out settings, do not run mfclaw
dry_run = False  # If True, only print out settings, do not run mfclaw

# what to do:
run_code = True
make_plots = False


# location for big files for different computer environments:
this_dir = os.getcwd()
HOME = os.environ['HOME']

if 'rjl/git' in this_dir:
    computer = 'rjl-laptop'
    #scratch_dir = this_dir.replace('rjl/git', 'rjl/scratch')

    # not using separate scratch_dir:
    scratch_dir = this_dir

elif '/home1' in this_dir:
    computer = 'tacc'
    #scratch_dir = this_dir.replace('/home1', '/scratch')
    try:
        SCRATCH = os.environ['SCRATCH']
        scratch_dir = this_dir.replace(HOME, SCRATCH)
    except:
        scratch_dir = this_dir  # if $SCRATCH not set

else:
    computer = 'unknown'
    scratch_dir = this_dir

print('scratch_dir = ',scratch_dir)

# where to put output for all the runs:
# (in a subdirectory of runs_dir named mfclaw_outputs)
runs_dir = os.path.abspath(scratch_dir)

# create scratch directory for output, if it doesn't exist:
os.system('mkdir -p %s' % runs_dir)

# path to mfclaw executable:
# (should agree with how EXE is set in Makefile used to compile)
if run_code:
    xmfclaw_path = '/Users/rjl/git/Forks/multifluid-dispersive-wave-collab/crater_code/5eqns_transport/mfrk2Dev/runs/gaussian_ring/xamr'
    if computer == 'tacc':
        # not using TACC yet...
        xmfclaw_path = '/work2/04137/rjl/CHTshare/clawpack-share/tacc/xmfclaw'
else:
    xmfclaw_path = None  # do not run mfclaw code


# Specify the list of h0 values to loop over for mfclaw runs:

h0vals = [1000]



if __name__ == '__main__':

    import sys

    # see the doc string at the top of this file for more details on
    # how to run this script

    # parse command line arguments:

    # always expect at least one argument NPROCS, how many to run in parallel
    try:
        nprocs = int(sys.argv[1])
    except:
        raise Exception('*** Missing integer argument nprocs on command line')


    print('\n--------------------------')
    if dry_run:
        # just print out settings, no runs...
        print('DRY RUN - settings in runclaw_makeplots_h0vals.py')

    print('Will run mfclaw for %i h0 values' % len(h0vals))
    print('list of h0 values: ', h0vals)
    print('output will go in \n    %s/mfclaw_outputs/' % runs_dir)
    print('plots will go in \n    %s/mfclaw_plots/' % runs_dir)
    print('nprocs = %i jobs will run simultaneously' % nprocs)
    print('OMP_NUM_THREADS = ', os.environ['OMP_NUM_THREADS'])
    print('xmfclaw executable:\n    ',xmfclaw_path)

    if dry_run:
        print('Set dry_run=False and re-execute to run mfclaw')
        print('--------------------------')
    else:
        foutfile = f'{runs_dir}/mfclaw_outputs/fortran_output.txt'
        poutfile = f'{runs_dir}/mfclaw_outputs/python_output.txt'
        print(f'Fortran output will go to \n    {foutfile}')
        print(f'Python output will go to \n    {poutfile}')

        # make list of dictionaries with parameters for each case:
        caselist = make_cases.make_all_cases_h0(runs_dir, xmfclaw_path,
                                                make_plots, h0vals)

        # copy setprob_local.txt into each outdir just created:
        for case in caselist:
            os.system(f'cp setprob_local.txt {case['outdir']}/')

        # run all cases using nprocs processors:
        multip_tools.run_many_cases_pool(caselist, nprocs,
                                         clawmultip_tools.run_one_case_clawpack)
