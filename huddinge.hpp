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

/*
class huddinge
{
public:

  std::vector<int>
  local_search(const std::string& pattern_, const std::string& text_)
  {
    pattern = pattern_;
    text = text_;
    plen = pattern.length();
    tlen = text.length();
    assert(plen <= tlen);

    err.resize(boost::extents[plen+1][tlen+1]);
    err[0][0] = err_t(0,0);

    for (int j=1; j < tlen+1; ++j)
      err[0][j] = insert(0, j);

    for (int i=1; i < plen+1; ++i)
      err[i][0] = del(i, 0);

    for (int i=1; i < plen+1; ++i) {
      for (int j=1; j < tlen+1; ++j) {
	err_t a = insert(i,j);
	err_t b = del(i,j);
	err_t c = match(i,j);

	int xerr = min3(a.get<0>(), b.get<0>(), c.get<0>());
	int yerr = min3(a.get<1>(), b.get<1>(), c.get<1>());
	assert(xerr != INT_MAX and yerr != INT_MAX);
	err[i][j] = err_t(xerr, yerr);
      }
    }

    std::vector<int> result(tlen+1);
    for (int j=0; j < tlen+1; ++j)
      result[j] = std::max(err[plen][j].get<0>(), err[plen][j].get<1>());
    return result;
  }

private:


  err_t
  insert(int i, int j)
  {
    if (i==0)
      return err_t(0,0);
    if (i==plen) {
      int errx, erry;
      boost::tie(errx, erry) = err[i][j-1];
      return err_t(errx + (text[j-1] != 'n'), erry);
    }
    return err_t(INT_MAX,INT_MAX);
  }

  err_t
  del(int i, int j) {
    if (j==0 or j==tlen) {
      int errx, erry;
      boost::tie(errx, erry) = err[i-1][j];
      return err_t(errx, erry + (pattern[i-1] != 'n'));
    }
    return err_t(INT_MAX,INT_MAX);
  }

  err_t
  match(int i, int j)
  {
    if (i==0 or j==0)
      return err_t(INT_MAX,INT_MAX);
    int errx, erry;
    boost::tie(errx, erry) = err[i-1][j-1];
    return err_t(errx + cerror(text[j-1], pattern[i-1]), erry + cerror(pattern[i-1], text[j-1]));
  }


  std::string pattern;
  std::string text;
  int plen;
  int tlen;
  boost::multi_array<err_t, 2> err;

} ;

*/
