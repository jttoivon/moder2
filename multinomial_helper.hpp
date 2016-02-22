/*

    MODER is a program to learn DNA binding motifs from SELEX datasets.
    Copyright (C) 2016  Jarkko Toivonen

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

typedef boost::unordered_map<std::string, std::vector<boost::tuple<int, int> > > string_to_tuple_type;


// Do not reject sequences with multiple occurrences of query strings
// compute the counts for the multinomial1 matrix
dmatrix
find_snips_multimer(const std::string& seed, const std::vector<std::string>& sequences, int hamming_distance);


boost::tuple<dmatrix,int,int>
find_snips_multimer_helper(const std::string& seed, const std::vector<std::string>& sequences);


string_to_tuple_type
get_n_neighbourhood(const std::string&seed, int n);


boost::tuple<dmatrix,int>
find_multinomial_n_suffix_array(const std::string& seed, const std::vector<std::string>& sequences, const suffix_array& sa, int n, bool use_multimer);

dmatrix
align_all(const std::vector<std::string>& sequences);
