ROOTCFLAGS    = $(shell root-config --cflags)
#ROOTCFLAGS    = -pthread -m32
ROOTLIBS      = $(shell root-config --libs)
ROOTGLIBS     = #$(shell root-config --glibs)
#ROOTINCS      = $(shell root-config --incdir)

BAT_ROOT = $(HOME)/datanas/athena/15.6.9/external/BAT/0.3.2/i686-slc4-gcc34
BAT_LIBS = -L$(BAT_ROOT)/lib -lBAT
BAT_INCS = -I$(BAT_ROOT)/include

CXX           = g++ -m32
CXXFLAGS      = -I$(ROOTINCS) -O2 -Wall -fPIC
LD            = g++
LDFLAGS       = #-g
SOFLAGS       = -shared

CXXFLAGS     += $(ROOTCFLAGS)
LIBS          = $(BAT_LIBS) $(ROOTLIBS) -lMinuit -lReflex -lCintex
GLIBS         = $(ROOTGLIBS)

#PYTHONINCS    = -I$(HOME)/local/include/python2.5
#PYTHONLIBS    = -L$(HOME)/local/lib/python2.5

SOURCES       = CrossSectionFinder.cpp runComBAT.cxx
OBJECTS       = CrossSectionFinder.o 

all: lib project


dictionary.cxx: selection.xml CrossSectionFinder.h
	@echo Generating dictionary 
	genreflex CrossSectionFinder.h --deep $(ROOTINCS) $(BAT_INCS) -s selection.xml -o $@


CrossSectionFinder.o: CrossSectionFinder.cpp CrossSectionFinder.h dictionary.cxx
	$(CXX) $(CXXFLAGS) $(BAT_INCS) -c $< -o $@

lib: dictionary.cxx CrossSectionFinder.o
	$(CXX) $(CXXFLAGS) -shared -Wl,-soname,libCrossSectionFinder CrossSectionFinder.h $(LIBS) $(BAT_LIBS) \
	CrossSectionFinder.o dictionary.cxx -o libCrossSectionFinder.so

clean:
	rm -f *~ *.o *.so *.o~ *.gch core dictionary.cxx

project: $(OBJECTS) 
	$(CXX) $(LDFLAGS) $(LIBS) $(OBJECTS) $(BAT_LIBS) runComBAT.cxx  -o runComBAT

