##$ COMPILER: suppprted compilers are ifort, gnu >v4.7
##$ use mpif90 for parallel compiler
FC=ifort

##$ PLATFORM: supported platform are intel, gnu
##$ if using a parallel compiler (mpif90) check
##$ platform first using the cmd: $mpif90 -show
PLAT=intel


##$ LOCATION OF THE scifor and dmft_tools DIRECTORIES
##$ is placed in different locations edit INCARGS here below
LIBDIR=/opt


##$ CHOOSE THE MODEL BY SELECTING THE PROGRAM DRIVER
#EXE=hc_graphene
#EXE=hc_haldane
#EXE=stripe_graphene


#EXE=bhz_2d
EXE=bhz_2d_edge
#EXE=bhz_3d
#EXE=bhz_afm_2d
#EXE=bhz_afm_3d

#EXE=tki_3d
#EXE=tki_2d



##$ SET THE LOCATION OF YOU PROGRAM DRIVER (default is ./drivers)
DIR =drivers


##$ SET THE LOCATION WHERE TO PLACE THE EXECUTABLE (default is $HOME/.bin)
DIREXE=$(HOME)/.bin


ifeq ($(PLAT),intel)
FFLAG +=
endif

ifeq ($(PLAT),gnu)
INCARGS=-I$(LIBDIR)/scifor/gnu/include -I$(LIBDIR)/dmft_tools/gnu/include
FFLAG +=-ffree-line-length-none $(INCARGS)
endif




##$ CHOOSE LINKING OPTIONS:
##$ 
##$ If you intend to use mkl:
##$ 
#MKLARGS=-lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -lm
#ARGS=-ldmftt -lscifor $(MKLARGS) 
##$ 
##$ ELSE:
##$ 
ARGS=  -ldmftt -lscifor_all




##$ REVISION SOFTWARE VARIABLES
##$ 
REV=$(shell git rev-parse HEAD)
BRANCH=_$(shell git rev-parse --abbrev-ref HEAD)
REV=$(shell git rev-parse HEAD)
VER = 'character(len=41),parameter :: revision = "$(REV)"' > revision.inc

ifeq ($(BRANCH),_master)
BRANCH=
endif


##$ Extends the implicit support of the Makefile to .f90 files
.SUFFIXES: .f90
