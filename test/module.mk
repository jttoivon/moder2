##########################################
#
# test programs
#
##########################################

TEST_PROGRAMS2=test_suffix_array_wrapper test_kmer_tools test_packed_string test_dinucleotide test_multinomial
TEST_PROGRAMS=$(addprefix test/, $(TEST_PROGRAMS2))


TESTCXXFLAGS=$(CXXFLAGS) -I.
TESTLDFLAGS=$(LDFLAGS) -lboost_unit_test_framework

.PHONY: test_individual
test_individual: $(TEST_PROGRAMS)
	for p in $(TEST_PROGRAMS) ; do echo Running $$p; ./$$p ; done

.PHONY: test
test: test/test_main
	test/test_main

.PHONY: test_clean
test_clean:
	rm -f test/*.o $(TEST_PROGRAMS)


test/%.o: test/%.cpp
	$(CXX) $(TESTCXXFLAGS) -c $< -o $@

test/test_suffix_array_wrapper.o: suffix_array_wrapper.hpp

##########################################
#
# single executable
#
##########################################

# TEST_ALL_OBJS=dinucleotide.o iupac.o common.o data.o parameters.o probabilities.o my_assert.o multinomial_helper.o suffix_array_wrapper.o matrix_tools.o orientation.o kmer_tools.o CPM03/difference_cover.o combinatorics.o seed_basic.o packed_string.o huddinge.o
# TEST_ALL_SOURCES=test/test_main.cpp test/test_dinucleotide.cpp test/test_suffix_array_wrapper.cpp test/test_multinomial.cpp test/test_kmer_tools.cpp test/test_packed_string.cpp
# test/test_main: $(TEST_ALL_OBJS) $(TEST_ALL_SOURCES:.cpp=.o)
# 	$(CXX) $(TESTCXXFLAGS) $+ -o test/test_main  $(TESTLDFLAGS)

TEST_ALL_OBJS=dinucleotide.o iupac.o common.o data.o parameters.o probabilities.o my_assert.o multinomial_helper.o suffix_array_wrapper.o matrix_tools.o orientation.o kmer_tools.o CPM03/difference_cover.o combinatorics.o huddinge.o bndm.o
TEST_ALL_SOURCES=test/test_main.cpp test/test_dinucleotide.cpp test/test_probabilities.cpp test/test_suffix_array_wrapper.cpp
test/test_main: $(TEST_ALL_OBJS) $(TEST_ALL_SOURCES:.cpp=.o)
	$(CXX) $(TESTCXXFLAGS) $+ -o test/test_main  $(TESTLDFLAGS)



##########################################
#
# multiple executables
#
##########################################


TEST_DINUCLEOTIDE_OBJS=dinucleotide.o iupac.o common.o data.o parameters.o probabilities.o my_assert.o multinomial_helper.o suffix_array_wrapper.o matrix_tools.o orientation.o kmer_tools.o CPM03/difference_cover.o combinatorics.o
test/test_dinucleotide: $(TEST_DINUCLEOTIDE_OBJS) test/test_dinucleotide.cpp
	$(CXX) -DSTAND_ALONE $(TESTCXXFLAGS) $+ -o $@ $(TESTLDFLAGS)


# TEST_SUFFIX_ARRAY_WRAPPER_OBJS=CPM03/difference_cover.o suffix_array_wrapper.o iupac.o common.o orientation.o
# test/test_suffix_array_wrapper: $(TEST_SUFFIX_ARRAY_WRAPPER_OBJS) test/test_suffix_array_wrapper.cpp
# 	$(CXX) -DSTAND_ALONE $(TESTCXXFLAGS) $+ -o $@ $(TESTLDFLAGS)

# TEST_KMER_TOOLS_OBJS=kmer_tools.o common.o orientation.o iupac.o data.o probabilities.o parameters.o my_assert.o \
#                      combinatorics.o
# test/test_kmer_tools: $(TEST_KMER_TOOLS_OBJS) test/test_kmer_tools.cpp
# 	$(CXX) -DSTAND_ALONE  $(TESTCXXFLAGS) $+ -o $@ $(TESTLDFLAGS)

# TEST_PACKED_STRING_OBJS=packed_string.o common.o orientation.o iupac.o data.o probabilities.o parameters.o my_assert.o combinatorics.o
# test/test_packed_string: $(TEST_PACKED_STRING_OBJS) test/test_packed_string.cpp
# 	$(CXX) -DSTAND_ALONE $(TESTCXXFLAGS) $+ -o $@ $(TESTLDFLAGS)



# TEST_MULTINOMIAL_OBJS=dinucleotide.o  iupac.o common.o data.o parameters.o probabilities.o my_assert.o multinomial_helper.o suffix_array_wrapper.o matrix_tools.o orientation.o kmer_tools.o CPM03/difference_cover.o combinatorics.o
# test/test_multinomial: $(TEST_MULTINOMIAL_OBJS) test/test_multinomial.cpp
# 	$(CXX) -DSTAND_ALONE  $(TESTCXXFLAGS) $+ -o $@ $(TESTLDFLAGS)
