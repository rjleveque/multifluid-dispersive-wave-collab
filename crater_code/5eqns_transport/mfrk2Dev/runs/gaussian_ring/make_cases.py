

def make_all_cases_h0(runs_dir, xmfclaw_path, make_plots, h0vals):
    """
    Output: *caselist*, a list of cases to be run.
    Each case should be dictionary of any parameters needed to set up an
    individual case.  These will be used by run_one_case.

    For this example, the ocean depth h0 is varied.

    If xmfclaw_path is not None, AMRClaw will be run for each of
    the h0 values, based on code in setrun_case.py.
    The setrun function in this directory should take an argument case
    so that parameters for each case can be passed in.

    If make_plots is True, plotting will also be done based on code in
    setplot_case.py.
    The setplot function in this directory should take an argument case
    so that parameters for each case can be passed in.

    A unique directory will be created for each run, with names
    based on h0, and residing within runs_dir.
    The _output and _plots directories will be within the run directory.

    xupper, yupper are hard-wired to fixed values below, but if a list
    of values is to be tested there could also be loops on these.
    """
    import os

    # Create a list of the cases to be run:
    caselist = []

    for h0 in h0vals:

        # Create all directories needed first, in case there's a problem:

        # name for each run e.g. h0_0900 if h=900. (depth as 0-padded integer)
        run_name = f'h0_{h0:04d}'

        runs_dir = os.path.abspath(runs_dir)
        outdir = os.path.join(runs_dir, 'mfclaw_outputs/_output_%s' \
                    % run_name)
        plotdir = os.path.join(runs_dir, 'mfclaw_plots/_plots_%s' \
                    % run_name)
        os.system('mkdir -p %s' % outdir)
        print('Created %s' % outdir)
        os.system('mkdir -p %s' % plotdir)
        print('Created %s' % plotdir)

        case = {}

        case['case_name'] = f'depth{h0:.0f}'
        case['outdir'] = outdir

        case['xclawcmd'] = xmfclaw_path   # executable created by 'make .exe'
                                          # or None, in which case code is
                                          # not run (plots may be made below)

        # setrun parameters:
        case['setrun_file'] = 'setrun_case.py'

        # parameter(s) that vary:
        case['h0'] = h0
        case['xupper'] = 20e3  # fixed for this set of cases
        case['yupper'] = 2e3   # fixed for this set of cases

        if make_plots:
            case['plotdir'] = plotdir
        else:
            case['plotdir'] = None  # if None, will not make plots

        case['setplot_file'] = 'setplot_case.py'

        #case['redirect_python'] = False  # for debugging send to stdout

        caselist.append(case)

    return caselist
