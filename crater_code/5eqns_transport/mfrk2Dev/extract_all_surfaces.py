from clawpack.clawutil.util import fullpath_import

# fix include path if necessary:
mfclaw_tools = fullpath_import('mfclaw_tools.py')

outdir = '4kDepthCraters_260126/_output300m44finer'
fname_nc = 'surface_RC300.nc'

mfclaw_tools.extract_surface(outdir, fname_nc)



