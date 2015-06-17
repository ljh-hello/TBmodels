#LOAD SCIFOR LIBRARY:
#FC=gfortran
#EXE=hc_graphene
#EXE=hc_haldane
#EXE=stripe_graphene
#EXE=bhz_2x2_fcc
#EXE=bhz_afm_fcc
#EXE=bhz_fcc
#EXE=tki_fcc_3d
#EXE=tki_fcc_2d
EXE=bhz_fcc_stripe

DIREXE=$(HOME)/.bin

.SUFFIXES: .f90

BRANCH=$(shell git rev-parse --abbrev-ref HEAD)
VER = 'character(len=41),parameter :: revision = "$(REV)"' > revision.inc

#MKLARGS=-lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -lm
#FFLAG +=-fpp -D_$(FPP)

INCS= -I/opt/scifor/gnu/include -I/opt/dmft_tools/gnu/include
ARGS= -ldmftt -lscifor -lfftpack -llapack -lblas -lminpack

all:
	@echo "======================== $@ ============================="
	$(FC) $(FFLAG) $(EXE).f90 -o $(DIREXE)/$(EXE) $(INCS) $(ARGS)

clean: 
	@echo "Cleaning:"
	@rm -fv $(DIREXE)/$(EXE)
	@rm -f *.mod *.o *~
#########################################################################
