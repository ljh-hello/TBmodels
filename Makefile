include make.inc

all:
	@echo ""
	@echo "compiling $(EXE).f90 .... "
	@echo ""
	$(FC) $(FFLAG) $(EXE).f90 -o $(DIREXE)/$(EXE) $(INCARGS) $(ARGS)
	@echo "Done"

debug:
	@echo ""
	@echo "compiling $(EXE).f90 .... "
	@echo ""
	$(FC) $(DFLAG) $(EXE).f90 -o $(DIREXE)/$(EXE) $(INCARGS) $(ARGS)
	@echo "Done"
clean: 
	@echo "Cleaning:"
	@rm -f *.mod *.o *~
#########################################################################
