
# Makefile for Clawpack code in this directory.
# This version only sets the local files and frequently changed
# options, and then includes the standard makefile pointed to by CLAWMAKE.
CLAWMAKE = $(CLAW)/clawutil/src/Makefile.common

# See the above file for details and a list of make options, or type
#   make .help
# at the unix prompt.

MFCLAW = /Users/rjl/git/Forks/multifluid-dispersive-wave-collab/crater_code/5eqns_transport/mfrk2Dev

# Adjust these variables if desired:
# ----------------------------------

CLAW_PKG = amrclaw                  # Clawpack package to use
EXE = xamr                          # Executable to create
SETRUN_FILE = setrun.py             # File containing function to make data
OUTDIR = _output
SETPLOT_FILE = setplot.py           # File containing function to set plots
PLOTDIR = _plots                            # Directory for plots

OVERWRITE ?= True                   # False ==> make a copy of OUTDIR first

# Environment variable FC should be set to fortran compiler, e.g. gfortran

# Compiler flags can be specified here or set as an environment variable
FFLAGS = -O2 -fopenmp -fbounds-check -std=legacy -ffpe-trap='invalid,overflow,zero' 
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
  qinit.f \
  $(MFCLAW)/qad.f \
  $(MFCLAW)/advanc.f \
  $(MFCLAW)/b4step2.f90 \
  $(MFCLAW)/cellave2.f \
  $(MFCLAW)/mapc2p.f \
  $(MFCLAW)/fss.f \
  $(MFCLAW)/rhoe_prho_ideal.f  \
  $(MFCLAW)/rhoe_prho_sg.f \
  $(MFCLAW)/prmtoc_local.f \
  $(MFCLAW)/zeroin_local.f \
  $(MFCLAW)/bc2amr.f90  \
  $(MFCLAW)/flag2refine2.f90  \
  $(MFCLAW)/grav_steady_sg.f  \
  $(MFCLAW)/dgrav_steady_sg.f  \
  $(MFCLAW)/ctomfd.f \
  $(MFCLAW)/zero1d_aux.f \
  $(MFCLAW)/estdt.f \
  $(MFCLAW)/setprob.f \
  $(MFCLAW)/setprob_local.f \
  $(MFCLAW)/setgrd.f \
  $(MFCLAW)/ginit.f \
  $(MFCLAW)/gfixup.f \
  $(MFCLAW)/step2.f90 \
  $(MFCLAW)/setaux.f90 \
  $(MFCLAW)/outmsh.f \
  $(MFCLAW)/outgrd.f \
  $(MFCLAW)/out1.f \
  $(MFCLAW)/out2.f \
  $(MFCLAW)/q2qlqr_interp_tvd2_qp.f \
  $(MFCLAW)/ctoqe.f \
  $(MFCLAW)/tfluct2.f \
  $(MFCLAW)/philim.f \
  $(MFCLAW)/bc2_lr.f \
  $(MFCLAW)/rpn2.f \
  $(MFCLAW)/src2.f \
  $(MFCLAW)/hllc.f \
  $(MFCLAW)/rpt2.f \
  $(MFCLAW)/state2_rotate3.f \
  $(MFCLAW)/roe.f \
  $(MFCLAW)/stepgrid.f \
  $(MFCLAW)/flux2.f \
  $(MFCLAW)/filpatch.f90  \

#  b4step2.f90 \
#  tick.f \
#  ./valout.f90 \
#  intfil.f90 \

#-------------------------------------------------------------------
# Include Makefile containing standard definitions and make options:
include $(CLAWMAKE)

