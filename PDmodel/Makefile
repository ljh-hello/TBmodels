#LOAD SCIFOR LIBRARY:
#======================================================
include sfmake.inc
#======================================================

DIREXE=$(HOME)/.bin

all: pd-model_fcc pd-model_bethe

pd-model_fcc:
	@echo "======================== $@ ============================="
	${FC} ${OPT} -static $@.f90 ${SFLIBS} -o $(DIREXE)/$@

pd-model_bethe:
	@echo "======================== $@ ============================="
	${FC} ${OPT} -static $@.f90 ${SFLIBS} -o $(DIREXE)/$@

clean: 
	@echo "Cleaning:"
	@rm -fv $(DIREXE)/tdd-pam
	@rm -f *.mod *.o *~
#########################################################################
