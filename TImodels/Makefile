include make.inc

all:
	@echo ""
	@echo "compiling $(EXE).f90 .... "
	@echo ""
	$(FC) $(FFLAG) $(EXE).f90 -o $(DIREXE)/$(EXE) $(INCS) $(ARGS)

clean: 
	@echo "Cleaning:"
	@rm -fv $(DIREXE)/$(EXE)
	@rm -f *.mod *.o *~
#########################################################################
