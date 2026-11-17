## Test programs for the electrolyte equilibrium solver (see README.md).
##
## Build the model first with the ASP makefile (ASP repository, build/ASPMakefile.linux.gfortran),
## then, from this directory:
##
##   make ASP_BUILD=<ASP build directory>
##
## ASP_BUILD is the directory that holds the model's ObjectFiles/ and ModuleFiles/.  This makes
##   asp_feedump  the ASP model with the solver-state dump attached (ASP_FEE_DUMP=<file> to dump)
##   FEEReplay    replays a dump with the solver the deck selects (writes FEEReplay.txt)
## The model build itself does not use anything in this directory.

ASP_BUILD ?= ../../build
OBJDIR     = $(ASP_BUILD)/ObjectFiles/
MODDIR     = $(ASP_BUILD)/ModuleFiles/

FC  = gfortran -m64 -DGFORTRAN -fallow-argument-mismatch
CC  = gcc
F90FLAGS = -O3 -m64 -frounding-math -fsignaling-nans -gdwarf-3 -fdefault-real-8 -fdefault-double-8 \
           -g -fbacktrace -ffree-form -ffree-line-length-none -I$(MODDIR) -J.

## The model objects, as listed in the ASP makefile, without the main program
MODELOBJ = $(addprefix $(OBJDIR), CDFNor.o indexx.o PPND16.o dgaus8.o DMiLay.o MIEV0.o MieSupport.o \
           ModelParameters.o InfrastructuralCode.o GridPointFields.o Time.o Chemistry.o \
           Gibbs-DuhemEvaluationFunction.o LagrangianAerosols.o Condensation.o soa_unidriv.o \
           soa_unidriva.o unifac.o OutputRoutines.o Coagulation.o ReadInputDeck.o StepASP.o \
           EulerianBoxDriver.o opkda2.o opkda1.o opkdmain.o)

all: asp_feedump FEEReplay

FEEHarness.o: FEEHarness.f90 $(OBJDIR)LagrangianAerosols.o
	$(FC) -c $(F90FLAGS) $<

FEEReplay.o: FEEReplay.f90 FEEHarness.o
	$(FC) -c $(F90FLAGS) $<

FEEDumpAttach.o: FEEDumpAttach.c
	$(CC) -c -O2 $<

asp_feedump: FEEHarness.o FEEDumpAttach.o $(MODELOBJ) $(OBJDIR)DevelopmentDriver.o
	$(FC) -o $@ $(OBJDIR)DevelopmentDriver.o FEEHarness.o FEEDumpAttach.o $(MODELOBJ)

FEEReplay: FEEReplay.o FEEHarness.o $(MODELOBJ)
	$(FC) -o $@ FEEReplay.o FEEHarness.o $(MODELOBJ)

clean:
	rm -f *.o *.mod asp_feedump FEEReplay

.PHONY: all clean
