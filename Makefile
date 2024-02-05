#!/bin/bash

CXX = g++

CXXFLAGS += -O3 -fPIC -std=c++11

# LHAPDF
#LHAPDFINCS = $(shell lhapdf-config --cppflags)
#LHAPDFLIBS = $(shell lhapdf-config --ldflags)

# APFEL++
APFELPPINCS = $(shell apfelxx-config --cppflags)
APFELPPLIBS = $(shell apfelxx-config --ldflags)

# HOPPET
HOPPETINCS = $(shell hoppet-config --cxxflags)
HOPPETLIBS = $(shell hoppet-config --libs)

# Now set up the compiler and link flags and libs
CXXFLAGS += $(APFELPPINCS) $(LHAPDFINCS) $(HOPPETINCS)
LDFLAGS  += $(APFELPPINCS) $(LHAPDFINCS) $(HOPPETINCS)

CLIBS += $(APFELPPLIBS) $(LHAPDFLIBS) $(HOPPETLIBS) -L/usr/local/Cellar/gcc/13.2.0/lib/gcc/current/ -lgfortran

install : all
all : StructureFunctionsJoint TabulateStructureFunctions ScaleVariations

StructureFunctionsJoint: StructureFunctionsJoint.o
	$(CXX) $(LDFLAGS) -o $@ $< $(CLIBS)

TabulateStructureFunctions: TabulateStructureFunctions.o
	$(CXX) $(LDFLAGS) -o $@ $< $(CLIBS)

ScaleVariations: ScaleVariations.o
	$(CXX) $(LDFLAGS) -o $@ $< $(CLIBS)

.SUFFIXES : .cc .o .f .c .f90

.cxx.o:
	$(CXX) $(CXXFLAGS) -c $<

.f.o:
	$(FC) $(LDFLAGS) -c $<

.f90.o:
	$(FC) $(LDFLAGS) -c $<

clean:
	rm -rf *.lo *.o *.la StructureFunctionsJoint TabulateStructureFunctions ScaleVariations *~

