package = moder
version = 0.9.0
tarname = $(package)
distdir = $(tarname)-$(version)
prefix=/usr/local
datarootdir=$(prefix)/share
docdir=$(datarootdir)/doc/$(package)

ifdef BOOSTROOT
BOOSTINC=-I $(BOOSTROOT)/include
BOOSTLIB=-L $(BOOSTROOT)/lib
PROGOPT=$(BOOSTROOT)/lib/libboost_program_options.a
else
BOOSTINC=
BOOSTLIB=
PROGOPT=-lboost_program_options -lboost_system -lboost_filesystem
endif

CXXFLAGS= -DPACKAGE_VERSION=\"$(version)\" -Wall -Wno-sign-compare -g $(BOOSTINC)

NOOPENMP?=0
ifeq ($(NOOPENMP),0)
	CXXFLAGS += -fopenmp #-Wno-unknown-pragmas
else
	CXXFLAGS += -Wno-unknown-pragmas
endif

DEBUG?=0
ifeq ($(DEBUG),1)
#    CFLAGS =-g3 -gdwarf2 -DDEBUG
	CXXFLAGS +=-O0
else
ifeq ($(DEBUG),2)
	CXXFLAGS+=-O0 -pg
else
	CXXFLAGS+=-O3
endif
endif

#OPTIMIZE=-O3

# -I /usr/include/x86_64-linux-gnu/c++/4.7/

PROGRAMS=moder

CXX=g++
#CXX=$(HOME)/usr/bin/g++

LDFLAGS=$(BOOSTLIB) $(PROGOPT) #-L$(HOME)/usr/lib #-B/usr/bin/ld.gold

all: CPM03 $(PROGRAMS)


.PHONY: CPM03
CPM03:
	@$(MAKE) -C CPM03 all

test: test/test_suffix_array_wrapper


install: all
	install -d $(prefix)/bin
	install -d $(docdir)
	install -m 0755 moder $(prefix)/bin
	install -m 0644 README $(docdir)

dist: $(distdir).tar.gz

$(distdir).tar.gz: FORCE $(distdir)
	tar chf - $(distdir) | gzip -9 -c > $(distdir).tar.gz
	rm -rf $(distdir)

FORCE:
	-rm $(distdir).tar.gz &> /dev/null
	-rm -rf $(distdir) &> /dev/null

distcheck: $(distdir).tar.gz
	gzip -cd $+ | tar xvf -
	$(MAKE) -C $(distdir) all check clean
	rm -rf $(distdir)
	@echo "*** Package $(distdir).tar.gz\
          ready for distribution."

check: all
	./moder --multinomial --prior addone --cob 0-0 data/TFAP2A-head-1000.seq GGGCA > /dev/null
	@echo "*** ALL TESTS PASSED ***"

$(distdir):
	mkdir -p $(distdir)
	mkdir -p $(distdir)/CPM03
	mkdir -p $(distdir)/data	
	cp README $(distdir)
	cp COPYING $(distdir)
	cp Makefile $(distdir)
	cp heatmap.R $(distdir)
	cp moder.cpp $(distdir)
	cp common.cpp $(distdir)
	cp probabilities.cpp $(distdir)
	cp parameters.cpp $(distdir)
	cp matrix_tools.cpp $(distdir)
	cp my_assert.cpp $(distdir)
	cp combinatorics.cpp $(distdir)
	cp multinomial_helper.cpp $(distdir)
	cp bndm.cpp $(distdir)
	cp orientation.cpp $(distdir)
	cp data.cpp $(distdir)
	cp iupac.cpp $(distdir)
	cp suffix_array_wrapper.cpp $(distdir)
	cp iupac.hpp $(distdir)
	cp matrix.hpp $(distdir)
	cp vectors.hpp $(distdir)
	cp timing.hpp $(distdir)
	cp matrix_tools.hpp $(distdir)
	cp common.hpp $(distdir)
	cp data.hpp $(distdir)
	cp cob.hpp $(distdir)
	cp orientation.hpp $(distdir)
	cp probabilities.hpp $(distdir)
	cp parameters.hpp $(distdir)
	cp combinatorics.hpp $(distdir)
	cp my_assert.hpp $(distdir)
	cp multinomial_helper.hpp $(distdir)
	cp suffix_array_wrapper.hpp $(distdir)
	cp bndm.hpp $(distdir)
	cp CPM03/checker.hpp $(distdir)/CPM03
	cp CPM03/COPYING $(distdir)/CPM03
	cp CPM03/difference_cover.cpp $(distdir)/CPM03
	cp CPM03/difference_cover.hpp $(distdir)/CPM03
	cp CPM03/doubling.hpp $(distdir)/CPM03
	cp CPM03/Makefile $(distdir)/CPM03
	cp CPM03/partition.hpp $(distdir)/CPM03
	cp CPM03/README $(distdir)/CPM03
	cp CPM03/stringsort.hpp $(distdir)/CPM03
	cp CPM03/suffixsort.hpp $(distdir)/CPM03
	cp CPM03/test-suffixsort.cpp $(distdir)/CPM03
	cp CPM03/timing.hpp $(distdir)/CPM03
	cp data/TFAP2A-head-1000.seq $(distdir)/data


# programs


MODER_OBJS=moder.o common.o  probabilities.o parameters.o matrix_tools.o my_assert.o combinatorics.o\
	multinomial_helper.o bndm.o orientation.o data.o iupac.o CPM03/difference_cover.o suffix_array_wrapper.o
moder: $(MODER_OBJS)
	$(CXX) $(CXXFLAGS) $+ -o $@ $(LDFLAGS)





# test programs

TEST_SUFFIX_ARRAY_WRAPPER_OBJS=CPM03/difference_cover.o suffix_array_wrapper.o iupac.o common.o test/test_suffix_array_wrapper.o
test/test_suffix_array_wrapper: $(TEST_SUFFIX_ARRAY_WRAPPER_OBJS)
	$(CXX) $(CXXFLAGS) $(TEST_SUFFIX_ARRAY_WRAPPER_OBJS) -o test/test_suffix_array_wrapper $(LDFLAGS)

.PHONY: dist FORCE distcheck


########################
#
# main compilation units
#
########################


# additional dependencies
my_assert.o: my_assert.hpp



#########################
#
# other compilation units
#
#########################

other_units=matrix_tools.o bndm.o common.o  alignment.o parameters.o probabilities.o combinatorics.o my_assert.o  orientation.o multinomial_helper.o data.o iupac.o suffix_array_wrapper.o

$(other_units): %.o: %.hpp   # all non-main units depend on their respective header files


.PHONY: clean
clean:
	rm -f *.o $(PROGRAMS)
	@$(MAKE) -C CPM03 clean

.PHONY: showprograms
showprograms:
	@echo $(PROGRAMS) | tr ' ' '\n'
