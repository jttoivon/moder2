package = moder
version = 0.9.3
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

CXXFLAGS= -std=gnu++11 -DPACKAGE_VERSION=\"$(version)\" -Wall -Wno-sign-compare -g $(BOOSTINC)

NOOPENMP?=0
ifeq ($(NOOPENMP),0)
	CXXFLAGS += -fopenmp #-Wno-unknown-pragmas
else
	CXXFLAGS += -Wno-unknown-pragmas
endif

DEBUG?=0
OBJDIR=.
PRGPREFIX=
ifeq ($(DEBUG),1)
#    CFLAGS =-g3 -gdwarf2 -DDEBUG
	CXXFLAGS +=-O0
	OBJDIR=debug
	PRGPREFIX=debug_
	RESULT:=$(shell mkdir -p $(OBJDIR))
else
ifeq ($(DEBUG),2)
	CXXFLAGS+=-O0 -pg
	OBJDIR=debug
	PRGPREFIX=debug_
	RESULT:=$(shell mkdir -p $(OBJDIR))
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
export CXXFLAGS
export LDFLAGS

all: CPM03 $(addprefix $(OBJDIR)/$(PRGPREFIX), $(PROGRAMS))


.PHONY: CPM03
CPM03:
	@$(MAKE) -C CPM03 all

# .PHONY: test
# test:
# 	$(MAKE) -C test all
# test/test_suffix_array_wrapper




install: all
	install -d $(prefix)/bin
	install -d $(docdir)
	install -m 0755 moder $(prefix)/bin
	install -m 0644 README.md $(docdir)

dist: $(distdir).tar.gz

$(distdir).tar.gz: FORCE $(distdir)
	tar chf - $(distdir) | gzip -9 -c > $(distdir).tar.gz
	rm -rf $(distdir)

FORCE:
	-rm $(distdir).tar.gz &> /dev/null
	-rm -rf $(distdir) &> /dev/null

distcheck: $(distdir).tar.gz
	gzip -cd $+ | tar xvf -
	$(MAKE) -C $(distdir) all
	$(MAKE) -C $(distdir) check
	$(MAKE) -C $(distdir) clean
	rm -rf $(distdir)
	@echo "*** Package $(distdir).tar.gz\
          ready for distribution."

check: all
	./moder --prior addone --cob 0-0 data/TFAP2A-head-1000.seq GGGCA > /dev/null
	@echo "*** ALL TESTS PASSED ***"

$(distdir):
#	rm -rf $(distdir)
	mkdir -p $(distdir)
	mkdir -p $(distdir)/CPM03
	mkdir -p $(distdir)/data	
	cp README.md $(distdir)
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
	cp huddinge.cpp $(distdir)
	cp kmer_tools.cpp $(distdir)
	cp iupac.cpp $(distdir)
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
	cp type.hpp $(distdir)
	cp iupac.hpp $(distdir)
	cp huddinge.hpp $(distdir)
	cp kmer_tools.hpp $(distdir)
	cp unordered_map.hpp $(distdir)
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
	cp CPM03/dependencies.mk $(distdir)/CPM03
	cp data/TFAP2A-head-1000.seq $(distdir)/data


# programs


MODER_OBJS=all_pairs_huddinge.o common.o  probabilities.o parameters.o matrix_tools.o my_assert.o combinatorics.o\
	multinomial_helper.o bndm.o orientation.o data.o iupac.o suffix_array_wrapper.o kmer_tools.o huddinge.o
$(PRGPREFIX)moder: $(addprefix $(OBJDIR)/, $(MODER_OBJS)) CPM03/difference_cover.o
	$(CXX) $(CXXFLAGS) $(addprefix $(OBJDIR)/, $(MODER_OBJS)) CPM03/difference_cover.o -o $@ $(LDFLAGS)


ALL_PAIRS_HUDDINGE_OBJS=all_pairs_huddinge.o common.o  probabilities.o parameters.o matrix_tools.o my_assert.o combinatorics.o\
	multinomial_helper.o bndm.o orientation.o data.o iupac.o suffix_array_wrapper.o kmer_tools.o huddinge.o




# test programs

#TEST_SUFFIX_ARRAY_WRAPPER_OBJS=CPM03/difference_cover.o suffix_array_wrapper.o iupac.o common.o test_suffix_array_wrapper.o
#test/test_suffix_array_wrapper: $(TEST_SUFFIX_ARRAY_WRAPPER_OBJS)
#	$(CXX) $(CXXFLAGS) $(TEST_SUFFIX_ARRAY_WRAPPER_OBJS) -o test/test_suffix_array_wrapper $(LDFLAGS)

.PHONY: dist FORCE distcheck


####################################
#
# dependencies of compilation units
#
####################################

# Automatically create dependency files (*.d)
$(OBJDIR)/%.d: %.cpp
	@set -e; rm -f $@; \
         $(CXX) -MM $(CPPFLAGS) $< > $@.$$$$; \
         sed 's,\($*\)\.o[ :]*,'$(OBJDIR)'/\1.o $@ : ,g' < $@.$$$$ > $@; \
         rm -f $@.$$$$

-include $(addprefix $(OBJDIR)/, $(MODER_OBJS:.o=.d))


$(OBJDIR)/%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

#TESTCXXFLAGS=$(CXXFLAGS) -I. -O0
#test/%.o: test/%.cpp
#	$(CXX) $(TESTCXXFLAGS) -c $< -o $@





.PHONY: clean
clean:
	rm -f $(OBJDIR)/*.o $(addprefix $(PRGPREFIX), $(PROGRAMS))
	@$(MAKE) -C CPM03 clean

.PHONY: showprograms
showprograms:
	@echo $(PROGRAMS) | tr ' ' '\n'
