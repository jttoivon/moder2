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

PROGRAMS=multinomial simulate meme_em enrichment iterative_enrichment mdscanize scanner homo_dimer gapped_kmers kmers kmer_analysis matrix_distance submatrices search moder hetero_dimer p_value_main

CXX=g++
#CXX=$(HOME)/usr/bin/g++

LDFLAGS=$(BOOSTLIB) $(PROGOPT) #-L$(HOME)/usr/lib #-B/usr/bin/ld.gold

all: $(PROGRAMS)

test: test/test_suffix_array_wrapper

#SOURCES=main.cpp hello.cpp factorial.cpp
#OBJECTS=$(SOURCES:.cpp=.o)

install: moder
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
	$(MAKE) -C $(distdir) moder check clean
	rm -rf $(distdir)
	@echo "*** Package $(distdir).tar.gz\
          ready for distribution."

check: moder
	./moder --multinomial --prior addone --cob 0-0 data/TFAP2A-head-1000.seq GGGCA > /dev/null
	@echo "*** ALL TESTS PASSED ***"

$(distdir):
	mkdir -p $(distdir)
	mkdir -p $(distdir)/CPM03
	mkdir -p $(distdir)/data	
	cp README $(distdir)
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

MATRIX_DISTANCE_OBJS=matrix_distance.o common.o probabilities.o my_assert.o parameters.o \
	matrix_tools.o combinatorics.o
matrix_distance: $(MATRIX_DISTANCE_OBJS) matrix.hpp
	$(CXX) $(CXXFLAGS) $(MATRIX_DISTANCE_OBJS) -o $@ $(LDFLAGS)

HUDDINGE_DISTANCE_OBJS=huddinge.o huddinge_distance.o data.o common.o iupac.o orientation.o CPM03/difference_cover.o suffix_array_wrapper.o
huddinge_distance: $(HUDDINGE_DISTANCE_OBJS)
	$(CXX) $(CXXFLAGS) $+ -o $@ $(LDFLAGS)

MEME_OBJS=meme_em.o common.o  probabilities.o parameters.o matrix_tools.o my_assert.o combinatorics.o iupac.o #alignment.o  data.o
meme_em: matrix.hpp $(MEME_OBJS)
	$(CXX) $(CXXFLAGS) $(MEME_OBJS) -o $@ $(LDFLAGS)

MODER_OBJS=moder.o common.o  probabilities.o parameters.o matrix_tools.o my_assert.o combinatorics.o\
	multinomial_helper.o bndm.o orientation.o data.o iupac.o CPM03/difference_cover.o suffix_array_wrapper.o
moder: $(MODER_OBJS)
	$(CXX) $(CXXFLAGS) $+ -o $@ $(LDFLAGS) #suffix_array/SAu.a 

ENRICHMENT_OBJS=enrichment.o common.o  probabilities.o parameters.o matrix_tools.o \
	kmp.o data.o dependence.o my_assert.o lambda.o combinatorics.o
enrichment: $(ENRICHMENT_OBJS) matrix.hpp
	$(CXX) $(CXXFLAGS) $(ENRICHMENT_OBJS) -o $@ $(LDFLAGS)

ITERATIVE_ENRICHMENT_OBJS=iterative_enrichment.o common.o matrix_tools.o probabilities.o \
	parameters.o kmp.o data.o dependence.o my_assert.o combinatorics.o
iterative_enrichment: $(ITERATIVE_ENRICHMENT_OBJS) matrix.hpp
	$(CXX) $(CXXFLAGS) $(ITERATIVE_ENRICHMENT_OBJS) -o $@ $(LDFLAGS)

MULTINOMIAL_OBJS=matrix_tools.o common.o parameters.o kmp.o data.o\
	dependence.o bndm.o probabilities.o iupac.o\
	my_assert.o combinatorics.o lambda.o kmer_tools.o multinomial_helper.o packed_string.o CPM03/difference_cover.o \
	suffix_array_wrapper.o orientation.o
multinomial: dinucleotide.o multinomial.o $(MULTINOMIAL_OBJS) ahocorasick/libahocorasick.a
	$(CXX) $(CXXFLAGS) $+ -o $@ $(LDFLAGS) -lpthread ahocorasick/libahocorasick.a

libmultinomial.la: $(MULTINOMIAL_OBJS:.o=.lo) ahocorasick/libahocorasick.a
	libtool --mode=link $(CXX) $(CXXFLAGS) $(MULTINOMIAL_OBJS:.o=.lo)  -o libmultinomial.la $(LDFLAGS) -lpthread  -rpath /home/jttoivon/selex/src


SUBMATRICES_OBJS=matrix_tools.o common.o parameters.o data.o\
	probabilities.o submatrices.o my_assert.o combinatorics.o
submatrices: matrix.hpp my_assert.hpp $(SUBMATRICES_OBJS)
	$(CXX) $(CXXFLAGS) $(SUBMATRICES_OBJS)  -o $@ $(LDFLAGS) -lpthread

GAPPED_KMERS_OBJS=gapped_kmers.o common.o parameters.o data.o matrix_tools.o \
	probabilities.o my_assert.o lambda.o stats.o combinatorics.o kmer_tools.o orientation.o
gapped_kmers: matrix.hpp my_assert.hpp $(GAPPED_KMERS_OBJS)
	$(CXX) $(CXXFLAGS) $(GAPPED_KMERS_OBJS)  -o $@ $(LDFLAGS) -lpthread

KMERS_OBJS=kmers.o common.o parameters.o data.o matrix_tools.o orientation.o \
	probabilities.o my_assert.o dinucleotide.o ahocorasick/libahocorasick.a \
	kmer_tools.o combinatorics.o lambda.o iupac.o CPM03/difference_cover.o suffix_array_wrapper.o
kmers: matrix.hpp my_assert.hpp $(KMERS_OBJS)
	$(CXX) $(CXXFLAGS) $(KMERS_OBJS)  -o $@ $(LDFLAGS) -lpthread

SEARCH_OBJS=search.o common.o parameters.o data.o \
	my_assert.o bndm.o suffix_array_wrapper.o \
	kmer_tools.o lambda.o iupac.o CPM03/difference_cover.o orientation.o
search: my_assert.hpp $(SEARCH_OBJS)
	$(CXX) $(CXXFLAGS) $(SEARCH_OBJS)  -o $@ $(LDFLAGS)

KMER_ANALYSIS_OBJS=kmer_analysis.o common.o parameters.o data.o matrix_tools.o \
	probabilities.o my_assert.o ahocorasick/libahocorasick.a \
	kmer_tools.o combinatorics.o bndm.o lambda.o huddinge.o packed_string.o iupac.o orientation.o 
kmer_analysis: matrix.hpp my_assert.hpp $(KMER_ANALYSIS_OBJS)
	$(CXX) $(CXXFLAGS) $(KMER_ANALYSIS_OBJS)  -o $@ $(LDFLAGS) -lpthread

HOMO_DIMER_OBJS=homo_dimer.o matrix_tools.o common.o parameters.o kmp.o data.o alignment.o \
	bndm.o probabilities.o cooccurrence.o cooccurrence_utils.o permutation_test.o\
	my_assert.o combinatorics.o log.o cob.o esko.o lambda.o kmer_tools.o vectors.o orientation.o\
	dynamic_programming_common.o bernoulli.o sum_method.o
homo_dimer: matrix.hpp alignment.hpp my_assert.hpp fvector.hpp common.hpp $(HOMO_DIMER_OBJS)
	$(CXX) $(CXXFLAGS) $(HOMO_DIMER_OBJS)  -o $@ $(LDFLAGS) -lpthread

HETERO_DIMER_OBJS=hetero_dimer.o matrix_tools.o common.o parameters.o kmp.o data.o alignment.o \
	bndm.o probabilities.o cooccurrence.o cooccurrence_utils.o permutation_test.o\
	my_assert.o combinatorics.o log.o cob.o esko.o lambda.o kmer_tools.o vectors.o orientation.o\
	dynamic_programming_common.o bernoulli.o sum_method.o
hetero_dimer: matrix.hpp alignment.hpp my_assert.hpp fvector.hpp common.hpp $(HETERO_DIMER_OBJS)
	$(CXX) $(CXXFLAGS) $(HETERO_DIMER_OBJS)  -o $@ $(LDFLAGS) -lpthread

SCANNER_OBJS=scanner.o my_assert.o matrix_tools.o common.o probabilities.o parameters.o combinatorics.o \
	orientation.o nucleosome.o cob.o iupac.o
scanner: matrix.hpp my_assert.hpp $(SCANNER_OBJS)
	$(CXX) $(CXXFLAGS) $(SCANNER_OBJS)  -o $@ $(LDFLAGS)

P_VALUE_MAIN_OBJS=my_assert.o matrix_tools.o common.o probabilities.o parameters.o combinatorics.o \
	iupac.o data.o p_value.o orientation.o p_value_main.o
p_value_main: $(P_VALUE_MAIN_OBJS)
	$(CXX) $(CXXFLAGS) -std=c++11 $(P_VALUE_MAIN_OBJS)  -o $@ $(LDFLAGS)

MDSCANIZE_OBJS=mdscanize.o matrix_tools.o common.o parameters.o kmp.o data.o \
	dependence.o bndm.o probabilities.o my_assert.o combinatorics.o
mdscanize: matrix.hpp $(MDSCANIZE_OBJS)
	$(CXX) $(CXXFLAGS) $(MDSCANIZE_OBJS)  -o $@ $(LDFLAGS)

SIMULATE_OBJS=simulate.o matrix_tools.o common.o cooccurrence_utils.o iupac.o orientation.o
simulate: matrix.hpp $(SIMULATE_OBJS)
	$(CXX) $(CXXFLAGS) $(SIMULATE_OBJS) -o $@ $(LDFLAGS)

# test programs

TEST_SUFFIX_ARRAY_WRAPPER_OBJS=CPM03/difference_cover.o suffix_array_wrapper.o iupac.o common.o test/test_suffix_array_wrapper.o
test/test_suffix_array_wrapper: $(TEST_SUFFIX_ARRAY_WRAPPER_OBJS)
	$(CXX) $(CXXFLAGS) $(TEST_SUFFIX_ARRAY_WRAPPER_OBJS) -o test/test_suffix_array_wrapper $(LDFLAGS)

.PHONY: dist FORCE distcheck

.PHONY: ahocorasickdummy
ahocorasick/libahocorasick.a: ahocorasickdummy
	@$(MAKE) -C ahocorasick

.PHONY: suffixarraydummy
suffix_array/SAu.a: suffixarraydummy
	@$(MAKE) -C suffix_array
	@$(MAKE) -C CPM03

#include ahocorasick.make

########################
#
# main compilation units
#
########################

main_units=matrix_distance.o simulate.o multinomial.o gapped_kmers.o kmers.o kmer_analysis.o scanner.o meme_em.o enrichment.o iterative_enrichment.o homo_dimer.o submatrices.o search.o

# $(main_units): %.o: %.hpp   # all main units depend on their respective header files

# additional dependencies
matrix_distance.o: matrix.hpp
simulate.o: matrix.hpp common.hpp
my_assert.o: my_assert.hpp
kmer_analysis.o: kmer_tools.hpp unordered_map.hpp

multinomial.o: suffix_array_wrapper.hpp

p_value.o: p_value.cpp
	$(CXX) $(CXXFLAGS) -std=c++11 -c $<

#########################
#
# other compilation units
#
#########################

other_units=permutation_test.o dependence.o matrix_tools.o kmp.o bndm.o cob.o vectors.o stats.o cooccurrence.o cooccurrence_utils.o common.o esko.o sum_method.o alignment.o parameters.o bernoulli.o dynamic_programming_common.o probabilities.o combinatorics.o dinucleotide.o lambda.o my_assert.o kmer_tools.o nucleosome.o orientation.o multinomial_helper.o huddinge.o packed_string.o data.o iupac.o suffix_array_wrapper.o

$(other_units): %.o: %.hpp   # all non-main units depend on their respective header files


%.sika: %.alku
	echo keno $< $@

%.lo: %.cpp
	libtool --mode=compile $(CXX) $(CXXFLAGS) -O -c $<

CPM03/difference_cover.lo: CPM03/difference_cover.cpp
	libtool --mode=compile $(CXX) $(CXXFLAGS) -O -c $< -o CPM03/difference_cover.o

libraryinstall: libmultinomial.la
	libtool --mode=install install -c libmultinomial.la /home/jttoivon/selex/src/libmultinomial.la

.PHONY: clean libraryinstall

clean:
	rm -f *.o $(PROGRAMS)
	@$(MAKE) -C ahocorasick clean
	@$(MAKE) -C suffix_array clean

.PHONY: showprograms

showprograms:
	@echo $(PROGRAMS) | tr ' ' '\n'
