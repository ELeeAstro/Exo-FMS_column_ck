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
	  call_twostr.o \
	  k_Rosseland_mod.o \
	  IC_mod.o \
    ce_interp_mod.o \
		ce_Burrows_mod.o \
	  ck_opacity_mod.o \
	  dry_conv_adj_mod.o \
		MLT_mod.o \
		sw_direct_mod.o \
		sw_adding_mod.o \
		sw_SDA_mod.o \
		sw_Toon_mod.o \
		sw_SH2_mod.o \
		sw_SH4_mod.o \
		sw_disort_ts_mod.o \
		lw_AA_E_mod.o \
		lw_AA_L_mod.o \
		lw_sc_linear_mod.o \
		lw_VIM_mod.o \
		lw_Toon_mod.o \
		lw_disort_ts_mod.o \
    FMS_RC.o

#		DISORT.o \
# 	sw_disort_mod.o \
#		lw_disort_mod.o \


# executable statement
EXECS  = ../FMS_RC

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
