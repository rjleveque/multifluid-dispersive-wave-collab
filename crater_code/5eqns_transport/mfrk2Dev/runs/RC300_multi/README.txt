
The command
    python runclaw_makeplots_h0vals.py <nprocs>

will run mfclaw for a set of ocean depths h0 with all other parameters fixed
and create separate output directories and plot directories for each,
as subdirectories of ./mfclaw_outputs and ./mfclaw_plots.

<nprocs> should be an integer indicating how many jobs to run in parallel.

Some things to adjust in this script before running this:

    Set dry_run to True initially,
        or to False to actually do runs and/or make plots

    Set run_code and make_plots
        to indicate what should be done (when dry_run==False)

    Set xmfclaw_path to the path to the executable xamr (after compiling)

    Set h0vals to the desired set of h0 values (currently [4000, 3000])

See the comments in runclaw_makeplots_h0vals.py for more details.

Also note that setrun_case.py is set for a 1-level run so this runs quickly
for testing.

------------
Other files:

make_cases.py 
    contains make_all_cases_h0, 
    which makes the caselist (list of case dictionaries).
    Note that it loops over h0vals but hard-wires the other case parameters
    'xupper' and 'yupper', which are also used in setrun_case.py
    (A version could be created to loop over one of these instead of h0)

setrun_case.py 
    contains setrun with an argument `case` that will be
    passed in when running a particular case.

setplot_case.py 
    contains setplot with an argument `case` that will be
    passed in when running a particular case.


