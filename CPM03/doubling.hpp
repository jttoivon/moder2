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
// Doubling algorithm for suffix array construction.
// Follows closely the algorithm described in 
//
//   N.J. Larsson and K. Sadakane. Faster Suffix Sorting.
//   Technical Report LU-CS-TR:99-214, Dept. of Computer Science,
//   Lund University, Sweden, 1999.
//======================================================================

#include "partition.hpp"

#include <iterator>
#include <algorithm>
#include <utility>

//----------------------------------------------------------------------
// functor isa_access
//
// Given a suffix array iterator computes the corresponding
// inverse suffix array entry to be used in comparisons.
//----------------------------------------------------------------------
template <typename ISAIterator, typename SAIterator>
class isa_access
  : public std::unary_function<SAIterator,
	   typename std::iterator_traits<ISAIterator>::value_type>
{
private:
  typedef typename std::iterator_traits<ISAIterator>::difference_type
          difference_type;
  const ISAIterator isa;
public:
  isa_access(ISAIterator s, difference_type lcp) : isa(s+lcp) {}
  typename std::iterator_traits<ISAIterator>::value_type
  operator() (SAIterator i) const {
    return isa[*i]; 
  }
};

//----------------------------------------------------------------------
// function update_group
//
// Update a group after doubling.
// WARNING: cannot handle an empty group
//----------------------------------------------------------------------
template <typename SAIterator, typename ISAIterator>
void update_group(SAIterator beg, SAIterator end,
		  SAIterator sa, ISAIterator isa)
{
  typedef typename std::iterator_traits<SAIterator>::value_type
    pos_type;
  pos_type groupnumber = std::distance(sa, end) - 1;
  if (end-beg > 1) {
    for (SAIterator i = beg; i != end; ++i) {
      isa[*i] = groupnumber;
    }
  } else {
    isa[*beg] = groupnumber;
    *beg = -1;
  }
}


//----------------------------------------------------------------------
// function selection_partition
//
// Performs doubling partitioning for small groups.
//----------------------------------------------------------------------
template <typename SAIterator, typename ISAIterator>
void selection_partition(SAIterator sa, ISAIterator isa,
			 SAIterator beg, SAIterator end, 
         typename std::iterator_traits<SAIterator>::difference_type lcp)
{
  typedef typename std::iterator_traits<SAIterator>::value_type
    pos_type;
  isa_access<ISAIterator, SAIterator> key(isa, lcp);
  while (end-beg > 1) {
    SAIterator min_end = beg+1;
    pos_type min_key = key(beg);
    for (SAIterator i = beg+1; i != end; ++i) {
      if (key(i) < min_key) {
	std::swap(*beg, *i);
	min_end = beg+1;
	min_key = key(beg);
      } else if (!(min_key < key(i))) {
	std::swap(*min_end, *i);
	++min_end;
      }
    }
    update_group(beg, min_end, sa, isa);
    beg = min_end;
  }
  if (beg != end) update_group(beg, end, sa, isa);
}


//----------------------------------------------------------------------
// function doubling_partition
//
// Given a group of suffixes with a common prefix of length lcp
// partition it to subgroups according to the next lcp characters.
// The next lcp characters are represented by a single value
// in the inverse suffix array.
//----------------------------------------------------------------------
template <typename SAIterator, typename ISAIterator>
void doubling_partition(SAIterator sa, ISAIterator isa,
			SAIterator beg, SAIterator end, 
          typename std::iterator_traits<SAIterator>::difference_type lcp)
{
  typedef typename std::iterator_traits<SAIterator>::value_type
    pos_type;

  while (end-beg > 10) {

    // ternary partition
    isa_access<ISAIterator, SAIterator> get_isa_entry(isa, lcp);
    SAIterator pivot = random_pivot(beg, end, std::less<pos_type>(), 
				    get_isa_entry);
    std::pair<SAIterator, SAIterator> midrange
      = ternary_partition(beg, end, pivot, std::less<pos_type>(), 
			  get_isa_entry);
    
    // recurse on <-part
    if (beg != midrange.first) {
      doubling_partition(sa, isa, beg, midrange.first, lcp);
    }
    
    // update =-part
    update_group(midrange.first, midrange.second, sa, isa);
    
    // loop on >-part
    beg = midrange.second;
  }

  selection_partition(sa, isa, beg, end, lcp);
}



//======================================================================
// function doubling_sort
//
// The main doubling algorithm
//======================================================================
template <typename SAIterator, typename ISAIterator>
void
doubling_sort(SAIterator sa, SAIterator sa_end,
	      ISAIterator isa, ISAIterator isa_end,
         typename std::iterator_traits<ISAIterator>::difference_type lcp)
{
  typedef typename std::iterator_traits<SAIterator>::difference_type 
    difference_type;

  difference_type size = std::distance(sa, sa_end);
  while (*sa > -size) {

    SAIterator beg = sa;
    difference_type nrunlen = 0;
    while (beg != sa_end) {
      difference_type val = *beg;
      if (val < 0) {
	// jump over already sorted part
	beg -= val;
	nrunlen += val;
      } else {
	// combine preceding sorted parts
	if (nrunlen < 0) {
	  *(beg+nrunlen) = nrunlen;
	  nrunlen = 0;
	}
	// partition unsorted group
	SAIterator end = sa + isa[val] + 1;
	doubling_partition(sa, isa, beg, end, lcp);
	beg = end;
      }
    }
    if (nrunlen < 0) {
      *(beg+nrunlen) = nrunlen;
    }

    lcp *= 2;
  }
}
