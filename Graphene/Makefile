#LOAD SCIFOR LIBRARY:
#======================================================
include sfmake.inc
#======================================================
#EXE=hc_graphene
#EXE=hc_haldane
EXE=stripe_graphene
DIREXE=$(HOME)/.bin

all:
	@echo "======================== $@ ============================="
	${FC} ${STD} $(EXE).f90 ${SFLIBS} -o $(DIREXE)/$(EXE)

clean: 
	@echo "Cleaning:"
	@rm -fv $(DIREXE)/$(EXE)
	@rm -f *.mod *.o *~
#########################################################################
