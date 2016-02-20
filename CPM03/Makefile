
CXX      = g++
CXXFLAGS = -O3 -Wall -DDCOVER=5 -DTIMING

PROGS  = test-suffixsort

.SUFFIXES: .cpp .o

.cpp.o:
	$(CXX) $(CXXFLAGS) $(INCL) -c $<

all: $(PROGS)

OBJS1 =  test-suffixsort.o difference_cover.o
test-suffixsort: $(OBJS1)
	$(CXX) $(CXXFLAGS) $(INCL) -o $@ $(OBJS1)


#--------------------------------------------------------------------------
# cleaning up
#--------------------------------------------------------------------------
clean: depend
	rm -f *.o core


#--------------------------------------------------------------------------
# dependencies
#--------------------------------------------------------------------------

OBJS = $(OBJS1)
SOURCES = $(OBJS:.o=.cpp)

depend:
	$(CXX) -M $(SOURCES) > dependencies.mk

-include dependencies.mk
