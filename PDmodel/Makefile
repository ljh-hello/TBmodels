.SUFFIXES: .f90

BRANCH=$(shell git rev-parse --abbrev-ref HEAD)
VER = 'character(len=41),parameter :: revision = "$(REV)"' > revision.inc

MKLARGS=-lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -lm

#FFLAG +=-fpp -D_$(FPP)

ARGS=-lscifor $(MKLARGS)

DIREXE=$(HOME)/.bin

all: pd-model_fcc pd-model_bethe pd-model-hyb_fcc

pd-model_fcc:
	@echo "======================== $@ ============================="
	${FC} ${OPT} -static-intel $@.f90 ${SFLIBS} -o $(DIREXE)/$@ $(ARGS)

pd-model-hyb_fcc:
	@echo "======================== $@ ============================="
	${FC} ${OPT} -static-intel $@.f90 ${SFLIBS} -o $(DIREXE)/$@ $(ARGS)

pd-model_bethe:
	@echo "======================== $@ ============================="
	${FC} ${OPT} -static-intel $@.f90 ${SFLIBS} -o $(DIREXE)/$@ $(ARGS)

clean: 
	@echo "Cleaning:"
	@rm -f *.mod *.o *~
#########################################################################
