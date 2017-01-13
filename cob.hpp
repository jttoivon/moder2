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
#ifndef COB_HPP
#define COB_HPP

#include "matrix.hpp"

#include <boost/multi_array.hpp>
//typedef boost::multi_array<double, 2> array_type;

// This is for iterating the boost::multi_array
#define MY_FOREACH(x, container)							\
  for (int x=container.index_bases()[0] ; x < container.index_bases()[0] + (int)container.shape()[0]; ++x)



typedef boost::multi_array<long int, 2> int_cob_type;
typedef boost::multi_array<double, 2> double_cob_type;

typedef int_cob_type int_array_type;
typedef double_cob_type double_array_type;

template<typename NumberType>
boost::multi_array<NumberType, 2>
make_cob(int L, int k)
{
  typedef typename boost::multi_array<NumberType, 2> cob_type;

  assert(L > 0 && k > 0);
  assert(L >= k);

  typedef typename cob_type::extent_range range;
  typename cob_type::extent_gen extents;
 
  int d_min = -k + 1;
  int d_max = L - 2*k;

  return cob_type(extents[3][range(d_min,d_max+1)]);
}

template<typename NumberType>
boost::multi_array<NumberType, 2>
make_hetero_cob(int L, int k1, int k2)
{
  typedef typename boost::multi_array<NumberType, 2> cob_type;

  assert(L > 0 && k1 > 0 && k2 > 0);
  assert(L >= k1);
  assert(L >= k2);

  typedef typename cob_type::extent_range range;
  typename cob_type::extent_gen extents;
 
  int d_min = -std::min(k1, k2) + 1;
  int d_max = L - k1 - k2;

  return cob_type(extents[4][range(d_min,d_max+1)]);
}

class coordinate_transform;


dmatrix
my_submatrix(const dmatrix& m, int pos, int width);



class coordinate_transform
{
public:
  coordinate_transform(int L_, int pos_, int k_short_, int k_);

  dmatrix
  operator()(const dmatrix& cob, double initial_value=0.0) const; 


private:
  int L;
  int pos;
  int k_short;
  int k;
};



double_array_type
matrix_to_array(const dmatrix& m, int min_pos);

int_array_type
matrix_to_array(const matrix<long int>& m, int min_pos);


dmatrix
array_to_matrix(const double_array_type& a);

matrix<long int>
array_to_matrix(const int_array_type& a);


#endif // COB_HPP
