
# Makefile for Clawpack code in this directory.
# This version only sets the local files and frequently changed
# options, and then includes the standard makefile pointed to by CLAWMAKE.
CLAWMAKE = $(CLAW)/clawutil/src/Makefile.common

# See the above file for details and a list of make options, or type
#   make .help
# at the unix prompt.


# Adjust these variables if desired:
# ----------------------------------

CLAW_PKG = amrclaw                  # Clawpack package to use
EXE = xamr                          # Executable to create
SETRUN_FILE = setrun.py             # File containing function to make data
OUTDIR = _outputc                           # Directory for output
SETPLOT_FILE = setplot.py           # File containing function to set plots
PLOTDIR = _plotsc                            # Directory for plots

OVERWRITE ?= True                   # False ==> make a copy of OUTDIR first

# Environment variable FC should be set to fortran compiler, e.g. gfortran

# Compiler flags can be specified here or set as an environment variable
FFLAGS ?= -O2 -fopenmp -fbounds-check -std=legacy -ffpe-trap='invalid,overflow,zero' 
LFLAGS = -g -fopenmp 

# ---------------------------------
# package sources for this program:
# ---------------------------------

AMRLIB = $(CLAW)/amrclaw/src/2d
include $(AMRLIB)/Makefile.amr_2d

# ---------------------------------------
# package sources specifically to exclude
# (i.e. if a custom replacement source 
#  under a different name is provided)
# ---------------------------------------

EXCLUDE_MODULES = \

EXCLUDE_SOURCES = \

# ----------------------------------------
# List of custom sources for this program:
# ----------------------------------------

MODULES = \

SOURCES = \
  qad.f \
  advanc.f \
  b4step2.f90 \
  qinit.f \
  cellave2.f \
  mapc2p.f \
  fss.f \
  rhoe_prho_ideal.f  \
  rhoe_prho_sg.f \
  prmtoc_local.f \
  zeroin_local.f \
  bc2amr.f90  \
  flag2refine2.f90  \
  grav_steady_sg.f  \
  dgrav_steady_sg.f  \
  ctomfd.f \
  zero1d_aux.f \
  estdt.f \
  setprob.f \
  setprob_local.f \
  setgrd.f \
  ginit.f \
  gfixup.f \
  step2.f90 \
  setaux.f90 \
  outmsh.f \
  outgrd.f \
  out1.f \
  out2.f \
  q2qlqr_interp_tvd2_qp.f \
  ctoqe.f \
  tfluct2.f \
  philim.f \
  bc2_lr.f \
 ./rpn2.f \
 ./src2.f \
 ./hllc.f \
 ./rpt2.f \
 ./state2_rotate3.f \
 ./roe.f \
  ./stepgrid.f \
  ./flux2.f \
 ./filpatch.f90  \

#  b4step2.f90 \
#  tick.f \
#  ./valout.f90 \
#  intfil.f90 \

#-------------------------------------------------------------------
# Include Makefile containing standard definitions and make options:
include $(CLAWMAKE)

