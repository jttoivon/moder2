/*
 * Copyright 2003 Stefan Burkhardt, Juha K"arkk"ainen
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
// Multikey quicksort for suffixes
//======================================================================

#include "partition.hpp"

#include <algorithm>
#include <iterator>
#include <functional>
#include <vector>


//----------------------------------------------------------------------
// function string_compare
//
// Compare two strings until a mismatch but no longer than limit 
// characters.  The sign of the result indicates the order of the 
// strings and the absolute value the distance of the mismatch from 
// the limit. Zero means that no mismatch was found within limit 
// characters.
//
// WARNING: A check for the end of a string is user's responsibility.
//----------------------------------------------------------------------
template <typename StringIterator, typename DifferenceType,
	  typename CharCompare>
inline DifferenceType
string_compare(StringIterator s1, StringIterator s2, 
	       DifferenceType limit, CharCompare less) {
  for ( ; limit > 0; --limit, ++s1, ++s2) {
    if (less(*s1,*s2)) return -limit;
    if (less(*s2,*s1)) return limit;
  }
  return 0;
}
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
template <typename StringIterator, typename DifferenceType>
inline DifferenceType
string_compare(StringIterator s1, StringIterator s2, DifferenceType limit) 
{
  typedef typename std::iterator_traits<StringIterator>::value_type
    char_type;
  return string_compare(s1, s2, limit, std::less<char_type>());
}


//----------------------------------------------------------------------
// const int suffix_insertion_sort_threshold_
//
// Groups of suffixes with size less or equal to this value are 
// sorted by insertion sort.
//----------------------------------------------------------------------
static const int suffix_insertion_sort_threshold_ = 15;

//----------------------------------------------------------------------
// function suffix_insertion_sort
// 
// Sort suffixes using insertion sort.
// A difference from a straightforward insertion sort is that
// the lengths of longest common prefixes of the already sorted 
// suffixes are stored and used in later insertions.
//----------------------------------------------------------------------
template <typename StringIterator, typename SArrayIterator,
	  typename GroupHandler, typename CharCompare>
void
suffix_insertion_sort(StringIterator str, StringIterator str_end,
		      SArrayIterator sa, SArrayIterator sa_end,
     typename std::iterator_traits<StringIterator>::difference_type limit,
		      GroupHandler handle_unsorted_group,
		      CharCompare less)
{
  typedef typename std::iterator_traits<StringIterator>::difference_type 
    str_diff_type;
  typedef typename std::iterator_traits<SArrayIterator>::difference_type 
    sa_diff_type;
  typedef typename std::iterator_traits<SArrayIterator>::value_type 
    pos_type;
  
  static std::vector<str_diff_type> lcp(suffix_insertion_sort_threshold_);
  lcp.resize(sa_end-sa);
  lcp.back() = -1;
  
  for (sa_diff_type i = sa_end-sa-2; i >= 0; --i) {
    pos_type i_pos = sa[i];
    StringIterator i_str = str + i_pos;
    StringIterator j_str = str + sa[i+1];
    str_diff_type i_limit = std::min(limit, str_end - i_str);
    str_diff_type ij_limit = std::min(i_limit, str_end - j_str);
    str_diff_type dist = string_compare(i_str, j_str, ij_limit, less);
    if (dist < 0 || (dist == 0 && i_limit == ij_limit)) {
      lcp[i] = ij_limit + dist;
    } else {
      sa[i] = sa[i+1];
      str_diff_type i_lcp = ij_limit - dist;

      sa_diff_type j = i+1;
      for ( ; lcp[j] >= i_lcp; ++j) {
	if (lcp[j] == i_lcp) {
	  j_str = str+sa[j+1];
	  ij_limit = std::min(i_limit, str_end-j_str);
	  dist = string_compare(i_str+i_lcp, j_str+i_lcp, 
				ij_limit-i_lcp, less);
	  if (dist < 0 || (dist == 0 && i_limit == ij_limit)) {
	    lcp[j] = ij_limit + dist;
	    break;
	  } else {
	    i_lcp = ij_limit - dist;
	  }
	}
	sa[j] = sa[j+1];
	lcp[j-1] = lcp[j];
      }

      sa[j] = i_pos;
      lcp[j-1] = i_lcp;
    }
  }

  // find and process unsorted groups
  typename std::vector<pos_type>::iterator beg, end, last = lcp.end()-1;
  for (beg = lcp.begin(); beg < last; ++beg) {
    if (*beg == limit) {
      end = beg;
      do { ++end; } while (*end == limit);
      handle_unsorted_group(sa+(beg-lcp.begin()), sa+(end-lcp.begin())+1);
      beg = end;
    }
  }
}
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
template <typename StringIterator, typename SArrayIterator,
	  typename GroupHandler>
void
suffix_insertion_sort(StringIterator str, StringIterator str_end,
		      SArrayIterator sa, SArrayIterator sa_end,
    typename std::iterator_traits<StringIterator>::difference_type limit,
		      GroupHandler handle_unsorted_group)
{
  typedef typename std::iterator_traits<StringIterator>::value_type
    char_type;
  suffix_insertion_sort(str, str_end, sa, sa_end, limit,
			handle_unsorted_group, std::less<char_type>());
}


//----------------------------------------------------------------------
// functor char_access
//
// Compute the character that represents a suffix in comparisons.
//----------------------------------------------------------------------
template <typename StringIterator, typename SArrayIterator>
class char_access
  : public std::unary_function<SArrayIterator,
	   typename std::iterator_traits<StringIterator>::value_type>
{
private:
  const StringIterator str;
public:
  char_access(StringIterator s) : str(s) {}
  const typename std::iterator_traits<StringIterator>::value_type &
  operator() (SArrayIterator i) const {
    return str[*i]; 
  }
};


//======================================================================
// function suffix_mkqsort
// function suffix_mkqsort_noempty
//
// Sort suffixes using multikey quicksort.
// These two functions call each other to implement multikey quicksort.
// Either can be used as the initial call. suffix_mkqsort_noempty
// places more restrictions on the input but can be slightly faster.
//======================================================================

//----------------------------------------------------------------------
// function suffix_mkqsort
//
// Sort suffixes using multikey quicksort.
// This is the more general form.
//----------------------------------------------------------------------
template <typename StringIterator, typename SArrayIterator,
	  typename GroupHandler, typename CharCompare>
void
suffix_mkqsort(StringIterator str, StringIterator str_end,
	       SArrayIterator sa, SArrayIterator sa_end,
    typename std::iterator_traits<StringIterator>::difference_type limit,
	       GroupHandler handle_unsorted_group,
	       CharCompare less)
{		  
  typedef typename std::iterator_traits<SArrayIterator>::value_type pos_type;

  while (sa_end-sa > suffix_insertion_sort_threshold_ && limit > 0) {

    // check for empty suffix(es)
    pos_type length = static_cast<pos_type>(str_end-str);
    sa = std::partition(sa, sa_end, 
			std::bind2nd(std::equal_to<pos_type>(), length));

    // this would check only for one empty suffix
    //SArrayIterator empty = std::find(sa, sa_end, length);
    //if (empty != sa_end) {
    //  std::swap(*sa, *empty);
    //  ++sa;
    //}

    // partition
    char_access<StringIterator, SArrayIterator> key(str);
    SArrayIterator pivot = random_pivot(sa, sa_end, less, key);
    std::pair<SArrayIterator,SArrayIterator> midrange;
    midrange = ternary_partition(sa, sa_end, pivot, less, key);

    // recurse on <-part
    if (midrange.first - sa > 1) {
      suffix_mkqsort_noempty(str, str_end, sa, midrange.first, 
			     limit, handle_unsorted_group, less);
    }

    // recurse on >-part
    if (sa_end - midrange.second > 1) {
      suffix_mkqsort_noempty(str, str_end, midrange.second, sa_end, 
			     limit, handle_unsorted_group, less);
    }

    // loop on =-part
    sa = midrange.first;
    sa_end = midrange.second;

    // move to the next position in the suffixes
    ++str;
    --limit;

  }

  if (sa_end-sa > 1) {
    if (limit == 0) {
      handle_unsorted_group(sa, sa_end);
    } else {
      suffix_insertion_sort(str, str_end, sa, sa_end, limit, 
			    handle_unsorted_group, less);
    }
  }
}
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
template <typename StringIterator, typename SArrayIterator,
	  typename GroupHandler>
void
suffix_mkqsort(StringIterator str, StringIterator str_end,
	       SArrayIterator sa, SArrayIterator sa_end,
     typename std::iterator_traits<StringIterator>::difference_type limit,
	       GroupHandler handle_unsorted_group)
{
  typedef typename std::iterator_traits<StringIterator>::value_type
    char_type;
  suffix_mkqsort(str, str_end, sa, sa_end, limit,
		 handle_unsorted_group, std::less<char_type>());
}


//----------------------------------------------------------------------
// function suffix_mkqsort_noempty
//
// Sort suffixes using multikey quicksort.
// Used as a subroutine of suffix_mkqsort but can also be called 
// directly.  Places more restrictions on input (see below) but
// can be slightly faster than suffix_mkqsort.
//
// The additional restrictions:
// PRECONDITION: limit>0 and there are no empty suffixes
//----------------------------------------------------------------------
template <typename StringIterator, typename SArrayIterator,
	  typename GroupHandler, typename CharCompare>
void
suffix_mkqsort_noempty(StringIterator str, StringIterator str_end,
		       SArrayIterator sa, SArrayIterator sa_end,
    typename std::iterator_traits<StringIterator>::difference_type limit,
		       GroupHandler handle_unsorted_group,
		       CharCompare less)
{		  
  //  typedef typename std::iterator_traits<SArrayIterator>::value_type pos_type;

  while (sa_end-sa > suffix_insertion_sort_threshold_) {

    // partition
    char_access<StringIterator, SArrayIterator> key(str);
    SArrayIterator pivot = random_pivot(sa, sa_end, less, key);
    std::pair<SArrayIterator,SArrayIterator> midrange;
    midrange = ternary_partition(sa, sa_end, pivot, less, key);

    // recurse on <-part
    if (midrange.first - sa > 1) {
      suffix_mkqsort_noempty(str, str_end, sa, midrange.first, 
			     limit, handle_unsorted_group, less);
    }

    // recurse on =-part
    if (midrange.second - midrange.first > 1) {
      suffix_mkqsort(str+1, str_end, midrange.first, midrange.second, 
		     limit-1, handle_unsorted_group, less);
    }

    // loop on >-part
    sa = midrange.second;
  }

  if (sa_end-sa > 1) {
    suffix_insertion_sort(str, str_end, sa, sa_end, limit, 
			  handle_unsorted_group, less);
  }
}
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
template <typename StringIterator, typename SArrayIterator,
	  typename GroupHandler>
void
suffix_mkqsort_noempty(StringIterator str, StringIterator str_end,
		       SArrayIterator sa, SArrayIterator sa_end,
      typename std::iterator_traits<StringIterator>::difference_type limit,
		       GroupHandler handle_unsorted_group)
{
  typedef typename std::iterator_traits<StringIterator>::value_type
    char_type;
  suffix_mkqsort_noempty(str, str_end, sa, sa_end, limit,
			 handle_unsorted_group, std::less<char_type>());
}
    
