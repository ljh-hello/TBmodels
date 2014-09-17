#LOAD SCIFOR LIBRARY:
#======================================================
#include sfmake.inc
#======================================================
#EXE=hc_graphene
#EXE=hc_haldane
#EXE=stripe_graphene
#EXE=bhz_2x2_fcc
#EXE=bhz_afm_fcc
EXE=bhz_fcc
DIREXE=$(HOME)/.bin

.SUFFIXES: .f90

BRANCH=$(shell git rev-parse --abbrev-ref HEAD)
VER = 'character(len=41),parameter :: revision = "$(REV)"' > revision.inc

MKLARGS=-lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -lm

#FFLAG +=-fpp -D_$(FPP)

ARGS=-lscifor $(MKLARGS)

all:
	@echo "======================== $@ ============================="
	$(FC) $(FFLAG) $(EXE).f90 -o $(DIREXE)/$(EXE) $(ARGS)

clean: 
	@echo "Cleaning:"
	@rm -fv $(DIREXE)/$(EXE)
	@rm -f *.mod *.o *~
#########################################################################
