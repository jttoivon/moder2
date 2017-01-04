/*

    MODER is a program to learn DNA binding motifs from SELEX datasets.
    Copyright (C) 2016, 2017  Jarkko Toivonen

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
#ifndef ORIENTATION_HPP
#define ORIENTATION_HPP

#include "matrix.hpp"
#include "boost/tuple/tuple.hpp"

#include <vector>
#include <set>
#include <map>
#include <cassert>

enum {HT=0, HH=1, TT=2, TH=3};

extern const char* orients[];


class string_to_orientation_type
{
public:

  string_to_orientation_type() 
  {
    for (int i=0; i < 4; ++i) 
      m[orients[i]] = i;
  }

  int
  operator()(const std::string& s)
  {
    return m[s];
  }

private:
  std::map<std::string, int> m;
};

static string_to_orientation_type string_to_orientation;   // WHY THIS ISN'T IN THE CPP FILE??????????????????

int
orientation(int o1, int o2);

int
orientation2(int o1, int o2);  // four different orientations

boost::tuple<dmatrix,dmatrix>
get_matrices_according_to_hetero_orientation(int o, const dmatrix& m1, const dmatrix& m2);

boost::tuple<std::string,std::string>
get_seeds_according_to_hetero_orientation(int o, const std::string& s1, const std::string& s2);



// Bits in bitvector "abcdefgh"
//                    76543210
// correspond to following:
// a = TF1 in normal direction in position 1
// b = TF1 in reverse direction in position 1
// c = TF2 in normal direction in position 1
// d = TF2 in reverse direction in position 1
// e = TF1 in normal direction in position 2
// f = TF1 in reverse direction in position 2
// g = TF2 in normal direction in position 2
// h = TF2 in reverse direction in position 2
//
//  First   Second
//  Pos     Pos
//|       |       | 
// a b c d e f g h
//|   |   |   |   |
// TF1 TF2 TF1 TF2
// + - + - + - + -
class hetero_orientation_class
{
public:

  hetero_orientation_class() : v(256, -1)
  {
    /*
    v[128+2] = HT;  // strand
    v[16+4]  = HT;  // reverse strand

    v[128+1] = HH;  // strand
    v[32+4]  = HH;  // reverse strand

    v[64+2]  = TT;  // strand
    v[16+8]  = TT;  // reverse strand

    v[64+1]  = TH;  // strand
    v[32+8]  = TH;  // reverse strand
    */

    for (int o=0; o < 4; ++o) {
      v[a[o][0]] = o;
      v[a[o][1]] = o;
    }
  }


  std::vector<int>
  list(int bits) const       // returns a subset of {HT,HH,TT,TH}
  {
    /*
    std::set<int> result;
    if (t(128+2, i) or t(16+4, i))
      result.insert(HT);
    if (t(128+1, i) or t(32+4, i))
      result.insert(HH);
    if (t(64+2, i) or t(16+8, i))
      result.insert(TT);
    if (t(64+1, i) or t(32+8, i))
      result.insert(TH);

    return std::vector<int>(result.begin(), result.end());
    */
    std::vector<int> result;
    for (int o=0; o < 4; ++o)
      if (this->contains(bits, o))
	result.push_back(o);

    return result;
  }


  std::vector<std::pair<int,int> >
  list2(int bits) const       // returns a subset of {HT,HH,TT,TH}
  {
    std::vector<std::pair<int,int> > result;
    for (int o=0; o < 4; ++o)
      if (this->contains(bits, o, +1))
	result.push_back(std::make_pair(o, +1));
      else if (this->contains(bits, o, -1))
	result.push_back(std::make_pair(o, -1));

    return result;
  }

  bool
  contains(int bits, int orientation, int dir=0) const // dir == -1, 1, or 0   that is, directed or directionless
  {                                                 // Direction is +1 if TF1 is in the first position
    assert(0 <= orientation and orientation < 4);
    assert(bits >= 0);
    assert(bits < 256);

    if (dir == 1)
      return t(a[orientation][0], bits);
    else if (dir == -1)
      return t(a[orientation][1], bits);
    else if (dir == 0)
      return t(a[orientation][0], bits) or t(a[orientation][1], bits);
    else return false;
  }

  /*
  int
  operator()(int i) const
  {
    assert(0 <= i and i < 256);
    return v[i];
  }
  */
private:

  static const int a[4][2];   // first index is orientation, second is strand (0) or reverse strand (1)

  bool
  t(int mask, int bits) const
  {
    return (mask & bits) == mask;
  }

  std::vector<int> v;
};


class homo_orientation_class
{
public:

  homo_orientation_class() : v(256, -1)
  {
    v[8+2]  = HT;
    v[8+1]  = HH;
    v[4+2]  = TT;
    v[4+1]  = HT;
  }


  std::vector<int>
  list(int i) const
  {
    std::vector<int> result;
    if (t(8+2, i) or t(4+1, i))
      result.push_back(HT);
    if (t(8+1, i))
      result.push_back(HH);
    if (t(4+2, i))
      result.push_back(TT);

    return result;
  }

  bool
  contains(int i, int orientation) const
  {
    assert(0 <= orientation and orientation < 3);
    assert(i >= 0);
    assert(i < 16);

    if (orientation == HT)
      return t(a[HT], i) or t(a[TH], i);
    else 
      return t(a[orientation], i);
  }

  /*
  int
  operator()(int i) const
  {
    assert(0 <= i and i < 16);
    return v[i];
  }
  */
private:
  static const int a[4];

  bool
  t(int mask, int code) const
  {
    return (mask & code) == mask;
  }

  std::vector<int> v;
};

extern hetero_orientation_class hetero_orientation;
extern homo_orientation_class get_homo_orientation;

#endif // ORIENTATION_HPP
