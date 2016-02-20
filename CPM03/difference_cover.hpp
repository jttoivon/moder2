/*
 * Copyright 2003 Juha K"arkk"ainen
 *
 * This file is part of DC-sort
 *
 * DC-sort is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version. 
 *
 * DC-sort is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *
 */

//======================================================================
// This file defines two classes:
//
// class difference_cover implements basic difference covers
//     modulo any power of two
// class dc_sampler implements the difference cover based sampling
//     used in the suffix array construction
//
// Part of the implementation of class difference_cover is in
// difference_cover.cpp
//======================================================================

#ifndef DIFFERENCE_COVER_HPP
#define DIFFERENCE_COVER_HPP

#include <vector>
#include <utility>
#include <stdexcept>

//----------------------------------------------------------------------
// class difference_cover 
//
// A difference cover D modulo v is a set of integers in the range [0,v) 
// such that for all i in [0,v), there exists j,k in D such that 
// i = k-j (mod v).
//
// This class represents a difference cover modulo any power of two.
// The sizes of the difference covers for small v are:
//    v   1  2  4  8  16  32  64  128  256  512  1024  ...
//   |D|  1  2  3  4   5   7   9   13   20   28    40  ...
// For v<=128, the size of the difference cover is minimal.
// For any v, the size is less than sqrt(1.5*v)+6 (sqrt(v) is a lower bound).
//
// Notes:
// - All arithmetic is done with unsigned int, because it guarantees
//   correct modular arithmetic including conversion to/from other
//   integral types.  This is a main reason why supporting values of v
//   other than powers of two would be difficult.
//----------------------------------------------------------------------
class difference_cover {

public:

  // create a difference cover modulo 2^lv
  explicit difference_cover(unsigned lv); // defined in difference_cover.cpp

  // basic dimensions
  unsigned modulus() const { return 1<<logmod; }
  unsigned size() const { return coversize; }

  // division by modulus
  unsigned mod(unsigned i) const { return i&mask; }
  template <typename IntType>
  IntType div(IntType i) const { return i>>logmod; }

  // mapping [0,modulus) <-> [0,size)
  unsigned rank(unsigned pos) const { return coverrank[pos&mask]; }
  unsigned pos(unsigned rank) const { return cover[rank]; }

  // compute k such that i+k and j+k are both in the cover
  unsigned offset(unsigned i, unsigned j) const {
    return (diff2pos[(j-i)&mask]-i)&mask;
  }

  typedef std::vector<unsigned>::const_iterator iterator;
  iterator begin() const { return cover.begin(); }
  iterator end() const { return cover.end(); }

private:

  unsigned logmod;              // log v
  unsigned mask;                // v-1 = 11...1 in binary
  unsigned coversize;           // |D|

  std::vector<unsigned> cover;          // D
  std::vector<unsigned> coverrank;      // rank of a value in D
  std::vector<unsigned> diff2pos;       // lookup table for finding a
                                        // pair with a given difference
};



//----------------------------------------------------------------------
// class dc_sampler
//
// Represents the periodic sample positions chosen according 
// to a difference cover.
//
// The main duties are:
// 1. Fill the sparse suffix array with the sample positions.
// 2. Provide a mapping pack() from sample positions to packed positions.
//    A packed position pack(i) of a sample position i is used as 
//    the corresponding position in the inverse sparse suffix array.
//    pack(i) must satisfy:
//    a) pack(i) is unique (different for each sample position i) and 
//       in the range [0,samplesize)
//    b) The difference pack(i+v)-pack(i), where v is the length of
//       the period (modulus of the difference cover), is the same for
//       all i. This difference is called the packed period.
// 3. Find the representative pair for two arbitrary positions i and j.
//    The representative pair is a pair of positions in the inverse
//    sparse suffix array that can be used for comparing the suffixes 
//    i and j if they have a common prefix of length period-1. 
//    The representative pair is computed by finding k in [0,period)
//    such that i+k and j+k are both sample positions and returning 
//    the packed positions corresponding to i+k and j+k.
//----------------------------------------------------------------------    
template <typename OffsetType>
class dc_sampler {
public:

  dc_sampler(OffsetType s, OffsetType v) : dcover(v), fullsize_(s), 
     samplesize_(dcover.div(s)*dcover.size()+dcover.rank(s))
  {}

  OffsetType fullsize() const { return fullsize_; }
  OffsetType samplesize() const { return samplesize_; }
  OffsetType period() const { return dcover.modulus(); }
  OffsetType packed_period() const { return dcover.size(); }
	
  OffsetType pack(OffsetType i) const {
      return dcover.div(i) * dcover.size() + dcover.rank(i);
  }

  // representative pair
  std::pair<OffsetType,OffsetType> 
  pair(OffsetType i, OffsetType j) const {
    OffsetType k = dcover.offset(i, j);
    return std::make_pair(pack(i+k), pack(j+k));
  }

  // fills a range with the sample positions
  template <typename OutputIterator>
  OutputIterator fill(OutputIterator it) const {
    difference_cover::iterator j = dcover.begin();
    OffsetType i = 0;
    OffsetType pos = i + static_cast<OffsetType>(*j);
    while (pos < fullsize()) {
      *it++ = pos;
      if (++j == dcover.end()) {
	i += dcover.modulus();
	j = dcover.begin();
      }
      pos = i + static_cast<OffsetType>(*j);
    }
    return it;
  }

private:

  const difference_cover dcover;
  const OffsetType fullsize_;   // size of the range where samples are chosen
  const OffsetType samplesize_;
};

#endif // DIFFERENCE_COVER_HPP
