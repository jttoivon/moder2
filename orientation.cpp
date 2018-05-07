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
#include "orientation.hpp"
#include "matrix.hpp"
#include "common.hpp"

#include <boost/tuple/tuple.hpp>
#include <cassert>

boost::tuple<dmatrix,dmatrix>
get_matrices_according_to_hetero_orientation(int o, const dmatrix& m1, const dmatrix& m2, bool use_rna)
{
  dmatrix x, y;
  if (use_rna) {
    assert(o >= 0 and o <=1);
    switch (o) {
    case HT: x = m1; y = m2; break;                      // HT
    case RNA_TH: x = m2; y = m1; break;  // TH
    }
  }
  else {
    assert(o >= 0 and o <=3);

    switch (o) {
    case HT: x = m1; y = m2; break;                      // HT
    case HH: x = m1; y = reverse_complement(m2); break;  // HH
    case TT: x = reverse_complement(m1); y = m2; break;  // TT
    case TH: x = reverse_complement(m1); y = reverse_complement(m2); break;  // TH
    }
  }
  
  return boost::make_tuple(x, y);
}

boost::tuple<std::string,std::string>
get_seeds_according_to_hetero_orientation(int o, const std::string& s1, const std::string& s2, bool use_rna)
{
  std::string x, y;

  if (use_rna) {
    assert(o >= 0 and o <=1);
    switch (o) {
    case HT: x = s1; y = s2; break;                      // HT
    case RNA_TH: x = s2; y = s1; break;  // TH
    }
  }
  else {
    assert(o >= 0 and o <=3);
    
    switch (o) {
    case HT: x = s1; y = s2; break;                      // HT
    case HH: x = s1; y = reverse_complement(s2); break;  // HH
    case TT: x = reverse_complement(s1); y = s2; break;  // TT
    case TH: x = reverse_complement(s1); y = reverse_complement(s2); break;  // TH
    }
  }
  return boost::make_tuple(x, y);
}

int
orientation(int o1, int o2)
{
  assert(o1 == 0 || o1 == 1);
  assert(o2 == 0 || o2 == 1);

  if (o1 == o2)
    return HT;
  if (o1 == 0 && o2 == 1)
    return HH;
  if (o1 == 1 && o2 == 0)
    return TT;

  assert(false);
}

int
orientation2(int o1, int o2)
{
  assert(o1 == 0 || o1 == 1);
  assert(o2 == 0 || o2 == 1);

  if (o1 == 0 && o2 == 0)
    return HT;
  if (o1 == 0 && o2 == 1)
    return HH;
  if (o1 == 1 && o2 == 0)
    return TT;
  if (o1 == 1 && o2 == 1)
    return TH;

  assert(false);
}

const char* dna_orients[4] = {"HT", "HH", "TT", "TH"};
const char* rna_orients[2] = {"HT", "TH"};
const char** orients = dna_orients;

const int hetero_orientation_class::a[4][2] = { {128+2,16+4}, {128+1,32+4}, {64+2,16+8}, {64+1,32+8}};

const int homo_orientation_class::a[4] = { 8+2, 8+1, 4+2, 4+1};

hetero_orientation_class hetero_orientation;
homo_orientation_class get_homo_orientation;


