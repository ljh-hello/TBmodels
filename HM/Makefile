.SUFFIXES: .f90

BRANCH=$(shell git rev-parse --abbrev-ref HEAD)
VER = 'character(len=41),parameter :: revision = "$(REV)"' > revision.inc

MKLARGS=-lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -lm

ARGS=-lscifor $(MKLARGS)

EXE=hm_2bands_hyb_fcc3d
DIREXE=$(HOME)/.bin


all:
	@echo "======================== $@ ============================="
	${FC} ${OPT} -static-intel ${EXE}.f90 -o $(DIREXE)/${EXE} $(ARGS)

clean: 
	@echo "Cleaning:"
	@rm -f *.mod *.o *~
#########################################################################
