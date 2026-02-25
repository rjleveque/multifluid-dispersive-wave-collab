
"""
Set up the plot figures, axes, and items to be done for each frame.

This module is imported by the plotting routines and then the
function setplot is called to set the plot parameters.

"""


import os

from clawpack.clawutil.util import fullpath_import
mfclaw_tools = fullpath_import('./mfclaw_tools.py') # fix path if needed

# set fname_nc to the file created by extract_all_surfaces.py
# or set to None if no such file and surface must be extracted from fort.b
#fname_nc = 'surface600m44finer.nc'
fname_nc = 'surface_RC300.nc'
#fname_nc = None

if fname_nc is not None:
    # load surfaces from a netCDF file that was already created using
    # extract_all_surfaces.py
    tf, find_frame, surface_data = mfclaw_tools.load_surface_nc(fname_nc)



#--------------------------
def setplot(plotdata=None):
#--------------------------

    """
    Specify what is to be plotted at each frame.
    Input:  plotdata, an instance of clawpack.visclaw.data.ClawPlotData.
    Output: a modified version of plotdata.

    """


    from clawpack.visclaw import colormaps

    if plotdata is None:
        from clawpack.visclaw.data import ClawPlotData
        plotdata = ClawPlotData()


    plotdata.clearfigures()  # clear any old figures,axes,items data

    plotdata.format = 'binary'
    #plotdata.format = 'ascii'

    # loop
    # --------------------------

    clims = [(0.00,4.0), (0,1000), (-2000,0), (-15e3, 15e3), (0,1e9), (0,1)]

    for m in range(6):
        plotfigure = plotdata.new_plotfigure(figno=101+m)
        #plotfigure.show = False

        # Set up for axes in this figure:
        plotaxes = plotfigure.new_plotaxes()
        plotaxes.xlimits = 'auto'
        plotaxes.ylimits = 'auto'
        #plotaxes.xlimits = [0,2000]
        #plotaxes.ylimits = [-2100,2100]
        plotaxes.title = f'q({m+1} ) '
        plotaxes.scaled = True      # so aspect ratio is 1

        # Set up for item on these axes:
        plotitem = plotaxes.new_plotitem(plot_type='2d_pcolor')
        plotitem.plot_var = m
        plotitem.pcolor_cmap = colormaps.blue_yellow_red
        plotitem.add_colorbar = True
        plotitem.show = True       # show on plot?
        plotitem.pcolor_cmin = clims[m][0]
        plotitem.pcolor_cmax = clims[m][1]
        #plotitem.amr_data_show = [1,0]  # which levels to plot data
        plotitem.amr_patchedges_show = [0,1,1]
        plotitem.amr_celledges_show = [0,0,0]

    # 2D pcolor plot if desired:
    # --------------------------

    plotfigure = plotdata.new_plotfigure(name='2D', figno=0)
    plotfigure.show = True

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.xlimits = 'auto'
    plotaxes.ylimits = 'auto'
    plotaxes.title = 'q(6) indicator'
    plotaxes.scaled = True      # so aspect ratio is 1

    # Set up for item on these axes:
    plotitem = plotaxes.new_plotitem(plot_type='2d_pcolor')
    plotitem.plot_var = 0
    plotitem.pcolor_cmap = colormaps.blue_yellow_red
    plotitem.add_colorbar = True
    plotitem.show = True       # show on plot?
    plotitem.pcolor_cmin = 0.
    plotitem.pcolor_cmax = 1.0
    #plotitem.amr_data_show = [1,0]  # which levels to plot data
    plotitem.amr_patchedges_show = [1,1,1]
    #plotitem.amr_celledges_show = [1,0,0]



    # Figure for surface plot
    # -----------------------

    plotfigure = plotdata.new_plotfigure(name='surface', figno=2)

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.xlimits = [0,6000]
    plotaxes.ylimits = [-550, 950]
    plotaxes.title = 'Surface plot'
    plotaxes.grid = True

    def plot_surface(current_data):
        from pylab import plot, legend
        if fname_nc is not None:
            r = surface_data.coords['r'].to_numpy()  # convert to nd array
            eta = surface_data.sel(t=current_data.t).to_numpy()
        else:
            # this takes a lot longer (integrates vertically):
            r,eta = mfclaw_tools.load_surface(current_data.frameno,
                                              plotdata.outdir)
        plot(r, eta, 'm')

    plotaxes.afteraxes = plot_surface

    #-----------------------------------------
    # Figures for gauges
    #-----------------------------------------
    plotfigure = plotdata.new_plotfigure(name='q', figno=300,
                    type='each_gauge')
    plotfigure.show = False
    plotfigure.clf_each_gauge = True

    plotaxes = plotfigure.new_plotaxes()
    plotaxes.xlimits = 'auto'
    plotaxes.ylimits = 'auto'
    plotaxes.title = 'Pressure'
    plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
    plotitem.plot_var = 0
    plotitem.plotstyle = 'b-'


    # Parameters used only when creating html and/or latex hardcopy
    # e.g., via clawpack.visclaw.frametools.printframes:

    plotdata.printfigs = True                # print figures
    plotdata.print_format = 'png'            # file format
    plotdata.print_framenos = 'all'          # list of frames to print
    plotdata.print_fignos = 'all'            # list of figures to print
    plotdata.html = True                     # create html files of plots?
    plotdata.html_homelink = '../README.html'   # pointer for top of index
    plotdata.latex = True                    # create latex file of plots?
    plotdata.latex_figsperline = 2           # layout of plots
    plotdata.latex_framesperline = 1         # layout of plots
    plotdata.latex_makepdf = False           # also run pdflatex?
    plotdata.parallel = True                 # make multiple frame png's at once

    return plotdata
