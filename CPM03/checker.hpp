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
// Implements functions that check the correctness
// of a suffix array.
//
// INPUT: 
//  1. StringIterator str, StringIterator str_end
//     The range [str, str_end) is the string.
//     Type requirements:
//      * StringIterator is a model of Random Access Iterator
//      * If comparison operator is not provided, StringIterator's
//        value type is a model of Strict Weakly Comparable
//     Preconditions:
//      * [str, str_end) is a valid range
//     Notes:
//      * The algorithms make minimum assumptions on the type of
//        characters. In particular, characters are not copied,
//        printed or compared for (in)equality.
//  2. SArrayIterator sa, SArrayIterator sa_end
//     The range [sa, sa_end) is the suffix array.
//     Type requirements:
//      * SArrayIterator is a model of Random Access Iterator
//      * SArrayIterator's value type is an integral type that
//        can represent all values in the range [0,str_end-str].
//     Preconditions:
//      * [sa, sa_end) is a valid range
//  3. (optional) CharOrder less
//     Order comparison operator for characters
//     Type requirements:
//      * CharOrder is a model of Strict Weak Ordering
//      * StringIterator's value type is convertiple to
//        CharOrder's argument type
//
// OUTPUT: Throws an exception unless the input satisfies the
//         the following conditions:
//  1. [sa, sa_end) contains a permutation of the values [0,str_end-str)
//     Note that this implies that sa_end-sa == str_end-str.
//  2. For all i,j in [0,str_end-str), i occurs before j in [sa, sa_end)
//     if and only if the suffix [str+i,str_end) is lexicographically
//     smaller than [str+j,str_end).
//
//======================================================================

#ifndef SUFFIXARRAY_CHECKER_HPP
#define SUFFIXARRAY_CHECKER_HPP

#include <stdexcept>
#include <sstream>
#include <algorithm>
#include <functional>

//======================================================================
// function range_check_suffix_array
//
// The first stage of suffix array checking common to all versions.
// Makes the following checks:
//  1. the size of the suffix array is the same as the length 
//     of the string
//  2. all the entries are in the range [0,size-1)
//======================================================================
template <typename StringIterator, typename SArrayIterator>
void
range_check_suffix_array(StringIterator str, StringIterator str_end,
			 SArrayIterator sa, SArrayIterator sa_end) 
{
  typedef typename std::iterator_traits<SArrayIterator>::value_type 
    pos_type;
  
  // test 1
  if (sa_end-sa != str_end-str) {
    std::ostringstream os;
    os << "SUFFIX ARRAY CHECKER: "
       << "suffix array size " << sa_end-sa
       << " does not match the string length " << str_end-str;
    throw std::logic_error(os.str());
  }
  if (sa == sa_end) return;

  // test 2
  pos_type size = static_cast<pos_type>(str_end - str);
  for (SArrayIterator i = sa; i != sa_end; ++i) {
    if (*i < 0 || *i >= size) {
      std::ostringstream os;
      os << "SUFFIX ARRAY CHECKER: "
	 << "suffix array entry sa[" << i-sa << "]=" << *i
	 << " out of the range [0," << size << ")";
      throw std::out_of_range(os.str());
    }
  }
}


//======================================================================
// function check_suffix_array_trivial
//
// The simple but (potentially very) slow suffix array checker.
// METHOD: check correct order by string comparisons
// COMPLEXITY 
// time: total length of the distinguishing prefixes
//       O(n^2) in the worst case
// space: O(1)
//======================================================================
template <typename StringIterator, typename SArrayIterator,
	  typename CharOrder>
void
check_suffix_array_trivial(StringIterator str, StringIterator str_end,
			   SArrayIterator sa, SArrayIterator sa_end,
			   CharOrder less)
{
  range_check_suffix_array(str, str_end, sa, sa_end);
  if (sa == sa_end) return;

  for (SArrayIterator i = sa+1; i != sa_end; ++i) {
    if (!std::lexicographical_compare(str+*(i-1), str_end, str+*i, str_end,
				      less)) {
      std::ostringstream os;
      os << "SUFFIX ARRAY CHECKER (TRIVIAL): "
	 << "suffixes in wrong order\n"
	 << "sa[" << i-1-sa << "]=" << *(i-1)
	 << " and sa[" << i-sa << "]=" << *i;
      throw std::logic_error(os.str());
    }
  }
}


//======================================================================
// function check_suffix_array_using_inverse
//
// Fast but not very space efficient suffix array checker.
// METHOD: does the following tests
//   1. range_check_suffix_array
//   2. correct order by first character
//   3. for every pair of consecutive suffixes sa[i-1] and sa[i]
//      starting with the same first character:
//      check that the following suffixes sa[i-1]+1 and sa[i]+1
//      occur in the same relative order in the suffix array.
// Only the correct suffix array can pass these tests.
// The following suffixes are found using the inverse suffix
// array constructed by the algorithm.
// COMPLEXITY 
// time: linear
// space: n integers for the inverse suffix array
//======================================================================
template <typename StringIterator, typename SArrayIterator,
	  typename CharOrder>
void
check_suffix_array_using_inverse(
		 StringIterator str, StringIterator str_end,
		 SArrayIterator sa, SArrayIterator sa_end,
		 CharOrder less )
{
  range_check_suffix_array(str, str_end, sa, sa_end);
  if (sa == sa_end) return;

  typedef typename std::iterator_traits<SArrayIterator>::value_type 
    pos_type;
  pos_type size = static_cast<pos_type>(str_end - str);

  // construct inverse suffix array
  std::vector<pos_type> isa(size);
  SArrayIterator i;
  for (i = sa; i != sa_end; ++i) { 
    isa[*i] = static_cast<pos_type>(i-sa); 
  }

  SArrayIterator prev;
  for (prev = sa, i = sa+1; i != sa_end; prev = i++) {

    // test 2
    if (less(str[*i], str[*prev])) {
      std::ostringstream os;
      os << "SUFFIX ARRAY CHECKER (USING INVERSE): "
	 << "suffixes in wrong order\n"
	 << "first character of suffix sa[" << prev-sa << "]=" << *prev
	 << " is larger than the first character of suffix sa["
	 << i-sa << "]=" << *i;
      throw std::logic_error(os.str());
    } else if ( !less(str[*prev], str[*i]) ) {

      // test 3
      if (*i == size-1) {
	std::ostringstream os;
	os << "SUFFIX ARRAY CHECKER (USING INVERSE): "
	   << "suffixes in wrong order\n"
	   << "the suffix of length one is not the first among"
	   << " suffixes starting with the same character";
	throw std::logic_error(os.str());
      }
      if ( *prev != size-1 && isa[*i+1] < isa[*prev+1] ) {
	std::ostringstream os;
	os << "SUFFIX ARRAY CHECKER (USING INVERSE): "
	   << "suffixes in wrong order\n"
	   << "sa[" << prev-sa << "]=" << *prev
	   << " and sa[" << i-sa << "]=" << *i
	   << " have the same first character but\n"
	   << "sa[" << isa[*prev+1] << "]=" << sa[isa[*prev+1]]
	   << " and sa[" << isa[*i+1] << "]=" << sa[isa[*i+1]]
	   << " are in different relative order";
	throw std::logic_error(os.str());
      }
    }	
  }
}


//------------------------------------------------------------------------
// functor group_order
//
// Comparison of group iterators for binary searching in the function
// check_suffix_array_compact (just below).
//
// This might not work for all implementations of STL.  
// I could not find any fully generic way to create
// an SArrayIterator pointing to a given value, which I needed as 
// the query for binary searching using STL lower_bound.  This comparison 
// functor supports comparing an iterator to a value and allows using
// a value as a query.
//------------------------------------------------------------------------
template <typename StringIterator, typename SArrayIterator, 
	  typename CharOrder>
class group_order {
private:
  typedef typename std::iterator_traits<SArrayIterator>::value_type 
      offset_type;
  const StringIterator str;
  CharOrder less;
public:
  group_order(StringIterator s, CharOrder l) 
    : str(s), less(l) {}
  
  bool operator() (SArrayIterator a, SArrayIterator b) const
  { return less(str[*a], str[*b]); }
  bool operator() (SArrayIterator a, offset_type ob) const
  { return less(str[*a], str[ob]); }
  bool operator() (offset_type oa, SArrayIterator b) const
  { return less(str[oa], str[*b]); }
  bool operator() (offset_type oa, offset_type ob) const
  { return less(str[oa], str[ob]); }
};


//======================================================================
// function check_suffix_array_compact
//
// A suffix array checker using little extra space.
// METHOD: Does the following tests:
//   1. range_check_suffix_array
//   2. correct order by first character
//   3. for all characters c:
//        Let i1,i2,... be the suffixes that start with c
//        in this order (they are consecutive in the suffix array).
//        Check that suffixes i1+1,i2+1,... occur in the same
//        relative order in the suffix array.
// Only the correct suffix array can pass these tests.
//
// The suffixes i1+1,i2+1,... in the last test are found by scanning 
// the suffix array looking for suffixes preceded by the character c.
// All the characters can be processed in one pass by maintaining
// an array of iterators, one for each character, identifying the current
// position in the sequence i1,i2,... .
// There are two further space saving optimizations:
// * The array does not include iterators for characters that occur 
//   only once.
// * The array of iterators may be limited in size. If necessary,
//   the characters are divided into batches of this size, and
//   each group is processed in a separate pass. The size limit
//   is controlled by the macro CHECKER_COMPACT_LIMIT:
//   - if CHECKER_COMPACT_LIMIT > 0, the value gives the limit
//   - if CHECKER_COMPACT_LIMIT < 0, the limit is 
//     max(1, n/-CHECKER_COMPACT_LIMIT)
//   - otherwise there is no limit
//
// COMPLEXITY: 
// n = length of string
// k = number of distinct characters that occur more than once
// s = limit for the extra space
// if k <= s, time: O(n log k),  extra space: k iterators + O(1)
// if k > s,  time O((kn/s) log s), extra space: s iterators + O(1)
//======================================================================
template <typename StringIterator, typename SArrayIterator,
	  typename CharOrder>
void
check_suffix_array_compact(StringIterator str, StringIterator str_end,
			   SArrayIterator sa, SArrayIterator sa_end,
  			   CharOrder less)
{
  range_check_suffix_array(str, str_end, sa, sa_end);

  typedef typename std::iterator_traits<SArrayIterator>::value_type 
    pos_type;
  pos_type size = static_cast<pos_type>(str_end - str);
  if (size <= 1) return;
  
  // check that suffixes are correctly ordered by the first character
  // also count characters that occur more than once
  pos_type char_count = 0;
  bool newchar = true;
  SArrayIterator prev, i = sa; 
  for (prev = i++; i != sa_end; prev = i++) {
    if (less(str[*i], str[*prev])) {
      std::ostringstream os;
      os << "SUFFIX ARRAY CHECKER (COMPACT): "
	 << "suffixes in wrong order\n"
	 << "first character of suffix sa[" << prev-sa << "]=" << *prev
	 << " is larger than the first character of suffix sa["
	 << i-sa << "]=" << *i;
      throw std::logic_error(os.str());
    } else if (less(str[*prev], str[*i])) {
      newchar = true;
    } else if (newchar) {
      ++char_count;
      newchar = false;
    }
  }

  // create the array of pointers
#if CHECKER_COMPACT_LIMIT > 0
  pos_type space_limit = CHECKER_COMPACT_LIMIT;
#elif CHECKER_COMPACT_LIMIT < 0
  pos_type space_limit = std::max(1, size/(-(CHECKER_COMPACT_LIMIT)));
#else
  pos_type space_limit = size;
#endif
  std::vector<SArrayIterator> groups(std::min(char_count, space_limit));
  typedef typename std::vector<SArrayIterator>::iterator group_iterator;


  // test 3
  i = sa;
  prev = i++;
  newchar = true;
  while (char_count > 0) {

    // compute the next batch of characters to check
    group_iterator group = groups.begin();
    for ( ; i != sa_end; prev = i++) {
      if (group == groups.end()) break;   // this batch is full
      if (less(str[*prev], str[*i])) {
	newchar = true;
      } else if (newchar) {
	--char_count;
	*group++ = prev;
	newchar = false;
      }
    }
    // remove garbage from the end if the batch is not full
    // (only relevant to the last batch)
    groups.erase(group, groups.end());
    
    // check that suffix of length 1 is first in its group
    group = std::lower_bound(groups.begin(), groups.end(), size-1,
                group_order<StringIterator, SArrayIterator,CharOrder>(
							    str, less));
    if (group != groups.end() && !less(str[size-1], str[**group]) ) {
      if (**group != size-1) {
	std::ostringstream os;
	os << "SUFFIX ARRAY CHECKER (COMPACT): "
	   << "suffixes in wrong order\n"
	   << "the suffix of length one is not the first among"
	   << " suffixes starting with the same character";
	throw std::logic_error(os.str());
      }
      ++*group;
    }
    
    // check the rest
    for (SArrayIterator j = sa; j != sa_end; ++j) {
      if (*j > 0) {
	group = std::lower_bound(groups.begin(), groups.end(), *j-1,
                group_order<StringIterator,SArrayIterator,CharOrder>(
							   str, less));
	if (group != groups.end() && !less(str[*j-1], str[**group]) ) {
	  if (**group != *j-1) {
	    std::ostringstream os;
	    os << "SUFFIX ARRAY CHECKER (COMPACT): "
	       << "suffix in wrong position\n"
	       << "either suffix sa[" << (*group)-sa << "]=" << **group
	       << " among suffixes starting with the same character\n"
	       << "    or suffix sa[" << j-sa << "]=" << *j
	       << " among suffixes preceded by the same character";
	    throw std::logic_error(os.str());
	  }
	  ++*group;
	  if (*group == sa_end || less(str[*j-1], str[**group])) { --*group; }
	}
      }
    }
  }
}	
    


//======================================================================
// function check_suffix_array
//
// The default suffix array checker
//======================================================================
template <typename StringIterator, typename SArrayIterator, 
	  typename CharOrder>
void
check_suffix_array(StringIterator str, StringIterator str_end,
		   SArrayIterator sa, SArrayIterator sa_end,
		   CharOrder less)
{

#if defined CHECKER_COMPACT
  check_suffix_array_compact(str, str_end, sa, sa_end, less);
#elif defined CHECKER_TRIVIAL
  check_suffix_array_trivial(str, str_end, sa, sa_end, less);
#elif defined CHECKER_USING_INVERSE
  check_suffix_array_using_inverse(str, str_end, sa, sa_end, less);
#elif defined CHECKER_DOUBLE
  check_suffix_array_compact(str, str_end, sa, sa_end, less);
  check_suffix_array_using_inverse(str, str_end, sa, sa_end, less);
#else
  check_suffix_array_compact(str, str_end, sa, sa_end, less);
#endif

}

// a version with default comparison operator
template <typename StringIterator, typename SArrayIterator>
void
check_suffix_array(StringIterator str, StringIterator str_end,
		   SArrayIterator sa, SArrayIterator sa_end)
{
  typedef typename std::iterator_traits<StringIterator>::value_type
    char_type;
  check_suffix_array(str, str_end, sa, sa_end, std::less<char_type>());
}

#endif // SUFFIXARRAY_CHECKER_HPP
