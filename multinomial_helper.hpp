/*

    MODER is a program to learn DNA binding motifs from SELEX datasets.
    Copyright (C) 2016, 2017  Jarkko Toivonen,
    Department of Computer Science, University of Helsinki

    MODER is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    MODER is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along
    with this program; if not, write to the Free Software Foundation, Inc.,
    51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

*/
#include "matrix.hpp"
#include "suffix_array_wrapper.hpp"

#include <boost/tuple/tuple.hpp>
#include <string>
#include <vector>
#include <boost/unordered_map.hpp>


extern int palindromic_index_limit;
extern bool use_palindromic_correction;
extern int extra_flank;
extern int cluster_threshold;

enum counting_type {
  all_occurrences,                      // no restrictions
  sequence_contains_one,                // count only those sequences that contain exactly one occurrence
  neighbourhood_contains_one,           // count only those occurrences that have no neighbours (a neighbour is an intersecting site)
  choose_one_per_cluster     // count only one occurrence from a cluster
};

static std::map< counting_type, const char * > counting_type_to_string = {
   {all_occurrences, "all_occurrences"},
   {sequence_contains_one, "sequence_contains_one"},
   {neighbourhood_contains_one, "neighbourhood_contains_one"},
   {choose_one_per_cluster, "choose_one_per_cluster"}
};

extern counting_type data_counting;
extern counting_type background_counting;

// This actually maps a string to vector of (position,nucleotide) pairs.
typedef boost::unordered_map<std::string, std::vector<boost::tuple<int, int> > > string_to_tuple_type;

boost::tuple<std::string,int>
most_common_pattern_multimer(const std::vector<std::string>& sequences, int k, std::string seed,
			     bool contains_N, int hamming_radius=0);

int
conflict_free_palindromic_index(int hamming_radius);


boost::tuple<std::string,int>
most_common_pattern_monomer(const std::vector<std::string>& sequences, int k, std::string seed = "",
			    int hamming_radius=0);


// Do not reject sequences with multiple occurrences of query strings
// compute the counts for the multinomial1 matrix
dmatrix
find_snips_multimer(const std::string& seed, const std::vector<std::string>& sequences, int hamming_distance);


boost::tuple<dmatrix,int,int>
find_snips_multimer_helper(const std::string& seed, const std::vector<std::string>& sequences);


string_to_tuple_type
get_n_neighbourhood(const std::string&seed, int n);

std::vector<std::pair<std::string, std::vector<boost::tuple<int, int> > > >
get_n_neighbourhood_in_vector(const std::string&seed, int n);

boost::tuple<dmatrix,int>
find_multinomial_n_suffix_array(const std::string& seed, const std::vector<std::string>& sequences, const suffix_array& sa, int n, bool use_multimer);

boost::tuple<dmatrix,int>
find_multinomial_n_background(const std::string& seed, const std::vector<std::string>& sequences, const std::vector<double>& bg,
			      int n, bool use_multimer);

dmatrix
align_all(const std::vector<std::string>& sequences);
