
"""
Set up the plot figures, axes, and items to be done for each frame.

This module is imported by the plotting routines and then the
function setplot is called to set the plot parameters.

This version has a case argument to make plots for different cases.
Set up assuming we vary only h0.
"""


import os


#-----------------------------------
def setplot(plotdata=None, case={}):
#-----------------------------------

    """
    Specify what is to be plotted at each frame.
    Input:  plotdata, an instance of clawpack.visclaw.data.ClawPlotData.
    Output: a modified version of plotdata.

    case is a dictionary that will be passed in, and it is assumed that
    this dictionary has at least the keys 'outdir', 'plotdir', 'h0'
    with the appropriate values for one particular case.
    """


    from clawpack.visclaw import colormaps

    if plotdata is None:
        from clawpack.visclaw.data import ClawPlotData
        plotdata = ClawPlotData()


    plotdata.clearfigures()  # clear any old figures,axes,items data

    # these parameters should be set properly when calling this setplot:

    plotdata.outdir = case['outdir']
    plotdata.plotdir = case['plotdir']
    h0 = case['h0'] # not actually used, but could be used for plotaxes.ylimits

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
    plotaxes.xlimits = [0,3000]
    plotaxes.ylimits = [-550, 950]
    plotaxes.title = 'Surface plot'
    plotaxes.grid = True

    def surface(current_data):
        # Return radius of each grid cell and surface elevation
        from numpy import nan, isnan, nanmin, nanmax
        h0 = 4000.  # ocean depth
        r = current_data.x[:,0]
        y = current_data.y
        level = current_data.level
        dy = y[0,1] - y[0,0]
        ylow = y.min() - dy/2.
        q = current_data.q
        zfa = q[5,:,:]  # q(6,:,:) in Fortran, called zfa in Keh-Ming's code
        # compute have as in Keh-Mings out1.f:
        #have = (1-zfa).sum(axis=1) * dy   # sum up in y
        have = nan*zfa[:,0]  # all nan to start
        for j in range(len(r)):
            #if (zfa[j,:].min() < 1e-6) and (zfa[j,:].max() > 0.9):
            if (zfa[j,:].min() < 0.5) and (zfa[j,:].max() > 0.5):
                # surface seems to exist in jth column of this patch
                have[j] = (1-zfa[j,:]).sum() * dy   # sum up in y on this patch
                assert have[j] > 0, '*** expected have[j]>0'
                have[j] = have[j] + (ylow+h0)  # add depth of water below patch
                have[j] = have[j] - h0
            else:
                have[j] = nan  # not needed now
        #import pdb; pdb.set_trace()
        #if not isnan(nanmax(have)):
            #print(f'r: {r.min():7.0f},{r.max():7.0f}  y:{y.min():7.0f},{y.max():7.0f}  have: {nanmin(have):7.0f},{nanmax(have):7.0f}')
            #print(f'r.min() = {r.min()}, r.max() = {r.max()}')
            #print(f'y.min() = {y.min()}, y.max() = {y.max()}')
            #print(f'have.min() = {have.min()}, have.max() = {have.max()}')

        return r,have

    # Set up for item on these axes: scatter of 2d data
    plotitem = plotaxes.new_plotitem(plot_type='1d_from_2d_data')
    #plotitem.show = False
    plotitem.map_2d_to_1d = surface
    #plotitem.plot_var = 0
    #plotitem.plotstyle = '-'
    plotitem.plotstyle = 'bo-'
    plotitem.kwargs = {'markersize':3}
    #plotitem.color = 'b'
    plotitem.amr_data_show = [1,0]  # which levels to plot data

    # Set up for item on these axes: scatter of 2d data
    plotitem = plotaxes.new_plotitem(plot_type='1d_from_2d_data')
    plotitem.map_2d_to_1d = surface
    #plotitem.plot_var = 0
    plotitem.plotstyle = '-'
    plotitem.color = 'r'
    plotitem.amr_data_show = [0,1]  # which levels to plot data

    # Set up for item on these axes: scatter of 2d data
    plotitem = plotaxes.new_plotitem(plot_type='1d_from_2d_data')
    plotitem.map_2d_to_1d = surface
    #plotitem.plot_var = 0
    plotitem.plotstyle = '-'
    plotitem.color = 'g'
    plotitem.amr_data_show = [0,0,1]  # which levels to plot data


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
