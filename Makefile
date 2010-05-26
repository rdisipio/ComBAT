ROOTCFLAGS    = $(shell root-config --cflags)
#ROOTCFLAGS    = -pthread -m32
ROOTLIBS      = $(shell root-config --libs)
ROOTGLIBS     = #$(shell root-config --glibs)
ROOTINCS      = $(shell root-config --incdir)
#ROOTINCS      = $(ROOTSYS)/include
#ROOTLIBS      = -L$(ROOTSYS)/lib -lCore -lCint -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lTree -lRint -lPostscript -lMatrix
-lPhysics -l$

CXX           = g++
CXXFLAGS      = -I$(ROOTINCS) -O2 -Wall -fPIC
LD            = g++
LDFLAGS       = -g
SOFLAGS       = -shared

CXXFLAGS     += $(ROOTCFLAGS)
LIBS          = $(ROOTLIBS)
GLIBS         = $(ROOTGLIBS)

PYTHONINCS    = -I$(HOME)/local/include/python2.5
PYTHONLIBS    = -L$(HOME)/local/lib/python2.5

BAT_ROOT = ~/datanas/athena/15.6.9/external/BAT/0.3.2/i686-slc4-gcc34
BAT_LIBS = $(BAT_ROOT)/lib
BAT_INCS = $(BAT_ROOT)/include

