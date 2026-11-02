## To compile, select the compile options below and enter 'make' into the terminal
## the command 'make clean' removes all compiled files and executable files for recompilation

# for gfortran Compiler
#======================
F90          = gfortran
F90LINKER    = gfortran
DEFS      =
#FFLAGS   	= -Og -pipe -g -fbacktrace -Wall -Wextra -pedantic -fcheck=all -Wconversion -fbounds-check
FFLAGS   = -O3 -pipe
INCLUDES  =
LFLAGS    = -llapack -lblas

# for ifort Compiler
#======================
#F90          = ifort
#F90LINKER    = ifort
#DEFS      =
##FFLAGS   =
#FFLAGS   = -O3
#INCLUDES  =
#LFLAGS    = $(FFLAGS)

# for nvfortran Compiler
#======================
# F90          = nvfortran
# F90LINKER    = nvfortran
# DEFS      =
# #FFLAGS   =
# FFLAGS   = -fast
# INCLUDES  =
# LFLAGS    = $(FFLAGS)

OBJECTS = \
	  WENO4_mod.o \
	  BB_integrate_mod.o \
	  call_twostr.o \
	  RDI1MACH.o \
	  ERRPACK.o \
	  LINPAK.o \
	  LAPACK.o \
	  BDREF.o \
	  DISORT.o \
	  k_Rosseland_mod.o \
	  IC_mod.o \
	  ce_interp_mod.o \
	  ce_Burrows_mod.o \
	  ck_opacity_mod.o \
	  Manabe_dry_adjust_mod.o \
	  Pierrehumbert_dry_adjust_mod.o \
	  MLT_mod.o \
	  sw_direct_mod.o \
	  sw_adding_mod.o \
	  sw_SDA_mod.o \
	  sw_Toon89_mod.o \
	  sw_disort_ts_mod.o \
	  sw_disort_mod.o \
	  sw_SH2_mod.o \
	  sw_SH4_mod.o \
	  sw_AD_mod.o \
	  sw_Liou88_mod.o \
	  lw_AA_E_mod.o \
	  lw_AA_L_mod.o \
	  lw_sc_linear_mod.o \
	  lw_sc_scat_ric_mod.o \
	  lw_VIM_E_mod.o \
	  lw_VIM_L_mod.o \
	  lw_Toon89_mod.o \
	  lw_SH2_mod.o \
	  lw_SH4_mod.o \
	  lw_disort_ts_mod.o \
	  lw_disort_mod.o \
	  lw_AD_mod.o \
	  lw_DDA_mod.o \
	  lw_Fu97_mod.o \
	  lw_TAS_mod.o \
	  lw_DFE_mod.o \
	  FMS_RC.o

# executable statement
EXECS  = ../FMS_RC

# A bare .SUFFIXES: clears make's built-in suffix list before we set our own.
# Without it our list is *appended* to the built-ins, and the built-in Modula-2
# rule (.mod.o) shadows .f90.o for every target whose stem matches a gfortran
# generated <module>.mod file. Make then believes those objects have no
# prerequisites, silently skips recompiling them, and links stale objects.
.SUFFIXES:
.SUFFIXES: .o .f90 .F .f

default: FMS_RC

FMS_RC:  $(OBJECTS)
	 $(F90LINKER) $(FFLAGS) $(OBJECTS) -o $(EXECS) $(LFLAGS)

clean:
	rm -f *.o *.mod *__genmod.f90 $(EXECS)

.f90.o:
	$(F90) $(FFLAGS) $(DEFS) -c $<

.F.o:
	$(F90) $(FFLAGS) $(DEFS) -c $<

.f.o:
	$(F90) $(FFLAGS) $(DEFS) -c $<	

# DISORT 4.0.99 is distributed in single precision. In a gas-giant column the
# deep layers carry fluxes of order 1e5 W m-2 while the internal flux is of
# order 1 W m-2, so single precision loses the signal we are trying to
# resolve and the integration eventually goes unstable. Promote the whole
# DISORT toolchain to double precision instead of editing 600 kB of F77.
# (call_twostr.F was hand-converted to REAL*8 for the same reason.)
# gfortran flags; for nvfortran use -r8, for ifort use -real-size 64.
DISORT_FLAGS = $(FFLAGS) -fdefault-real-8 -fdefault-double-8

DISORT.o: DISORT.f
	$(F90) $(DISORT_FLAGS) $(DEFS) -c $<
BDREF.o: BDREF.f
	$(F90) $(DISORT_FLAGS) $(DEFS) -c $<
ERRPACK.o: ERRPACK.f
	$(F90) $(DISORT_FLAGS) $(DEFS) -c $<
LINPAK.o: LINPAK.f
	$(F90) $(DISORT_FLAGS) $(DEFS) -c $<
LAPACK.o: LAPACK.f
	$(F90) $(DISORT_FLAGS) $(DEFS) -c $<
RDI1MACH.o: RDI1MACH.f
	$(F90) $(DISORT_FLAGS) $(DEFS) -c $<
