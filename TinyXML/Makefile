.cpp.o:
	$(CXX) -c $(CXXFLAGS) $(GDBFLAGS) $<

all: tinyxml.o tinystr.o tinyxmlparser.o tinyxmlerror.o

TinyXML: tinyxml.cpp tinystr.cpp tinyxmlparser.cpp tinyxmlerror.cpp
	$(CXX) -c $(CXXFLAGS) $(GDBFLAGS) $<

clean:
	rm -f *.o *.~
