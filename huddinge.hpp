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
#include <map>
#include <set>
#include <string>
#include <climits>
#include <cstdio>

#include <boost/multi_array.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/foreach.hpp>
#include <boost/algorithm/string.hpp>

typedef boost::tuple<int, int> err_t;
typedef std::map<std::string, err_t > str_dist_container;


typedef std::vector<std::string> path_type;
typedef boost::tuple<std::string, double> mypair;
typedef std::vector<mypair> counted_path_type;

std::vector<int>
huddinge_alignment(const std::string& x, const std::string& y);



std::vector<counted_path_type>
find_paths(const std::string& x, const std::string& y);

std::vector<counted_path_type>
check_for_longer_paths(std::string x, std::string y, int dmax, int wmin, int L, int reference_length);


double
get_scaled_count(const std::string& x, int reference_length);

int
huddinge_distance(const std::string& x, const std::string& y);

bool
one_contiguous_gap(std::string s);

int
min3(int a, int b, int c);

int
defined_bases(const std::string& x);



bool
cerror(char a, char b);

class huddinge_neighbourhood
{
public:

  huddinge_neighbourhood(const std::string& x_, int h_, int L_, int min_kmer_len_, int max_kmer_len_) :
    x(boost::to_upper_copy(x_)), h(h_), L(L_), U(boost::extents[L+1][x.length()+1]),
    min_kmer_len(min_kmer_len_), max_kmer_len(max_kmer_len_)
  {
    //    printf("Array U has shape %zu x %zu\n", 
    //	   U.shape()[0], U.shape()[1]);
    U[0][0][""] = boost::make_tuple(0,0);
    //    printf("Size of U[%i][%i] is %zu\n", 0, 0, U[0][0].size());
  }

  std::vector<boost::tuple<std::string, int> >
  compute(bool allow_only_one_gap=false);

private:

  str_dist_container
  insert(int i, int j);

  str_dist_container
  del(int i, int j);

  str_dist_container
  match(int i, int j);

  err_t
  my_get(const str_dist_container& h, const std::string& key);



  std::string x;
  int h;
  int L;
  boost::multi_array<str_dist_container, 2> U;
  int min_kmer_len;
  int max_kmer_len;
};


