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
// A collection of flexible, generic, STL-like algorithm components 
// for comparison-based distribution (i.e., quicksort-like) sorting.
//======================================================================

#ifndef PARTITION_HPP
#define PARTITION_HPP

#include <algorithm>
#include <utility>
#include <functional>
#include <iterator>
#include <cstdlib>   // for rand()

//======================================================================
// Most of the algorithms here have three versions. For example a sort
// function can have the following forms:
//   sort(beg, end);
//   sort(beg, end, less);
//   sort(beg, end, less, key);
// The first two forms are standard STL. The third form adds an
// argument key, which is a unary function.  Where the second
// form would compare two iterators i and j as less(*i, *j),
// the third form does less(key(i), key(j)) instead. This allows
// moving some of the computation related to the comparison from
// less() into key().  When the same iterator is used in many
// comparisons (for example the pivot), key() needs to be evaluated
// only once. 
// 
// For example, to sort objects of type T with a member function
// size() returning int into increasing order of size() do
// (TODO: test/expand this, use vector of vectors)
//
//   sort(beg, end, std::less<unsigned>(), std::mem_fun(&T::size));
//
// key() may return a pointer or reference to (a member of) the object.
// All the algorithms assume that the value returned by key() stays
// valid only as long as the object is not moved.
//
// TYPE REQUIREMENTS
// For example:
// template <typename Iterator, typename Compare, typename KeyAccess>
// void sort(Iterator, Iterator, Compare, KeyAccess);
//  * KeyAccess is a model of Adaptable Unary Function
//  * Iterator is convertible to KeyAccess::argument_type
//  * KeyAccess::result_type is convertible to the argument type
//    of Compare.
//======================================================================

//----------------------------------------------------------------------
// functor copy_dereference
// functor const_dereference
// 
// Default functors to be used as key().
// (TODO: which one should be used as the default?)
// (TODO: is const_dereference::result_type a problem?)
//----------------------------------------------------------------------
template <typename Iterator>
class const_dereference
  : public std::unary_function<Iterator, 
          const typename std::iterator_traits<Iterator>::value_type &>
{
public:
  const typename std::iterator_traits<Iterator>::value_type &
  operator() (Iterator i) const { return *i; }
};

template <typename Iterator>
class copy_dereference 
  : public std::unary_function<Iterator, 
		 typename std::iterator_traits<Iterator>::value_type>
{
public:
  typename std::iterator_traits<Iterator>::value_type
  operator() (Iterator i) const { return *i; }
};

template <typename Iterator>
class default_key : public copy_dereference<Iterator> {};



//======================================================================
//
// Pivot selection
//
// (TODO: non-randomized pivot selectors)
//======================================================================

//----------------------------------------------------------------------
// function iter_median3
//
// Median of three
//----------------------------------------------------------------------
template <typename Iterator, typename Compare, typename KeyAccess>
Iterator
iter_median3(Iterator a, Iterator b, Iterator c, 
	     Compare less, KeyAccess key) {
  typedef typename KeyAccess::result_type key_type;
  key_type ka = key(a), kb = key(b), kc = key(c);
  if (less(ka, kb))
    if (less(kb, kc)) return b;
    else if (less(ka, kc)) return c;
    else return a;
  else if (less(ka, kc)) return a;
  else if (less(kb, kc)) return c;
  else return b;
};
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

template <typename Iterator, typename Compare>
Iterator
iter_median3(Iterator a, Iterator b, Iterator c, Compare less) {
  return iter_median3(a, b, c, less, default_key<Iterator>());
}

template <typename Iterator>
Iterator
iter_median3(Iterator a, Iterator b, Iterator c) {
  typedef typename std::iterator_traits<Iterator>::value_type key_type;
  return iter_median3(a, b, c, std::less<key_type>(), 
		      default_key<Iterator>());
}


//----------------------------------------------------------------------
// function random_pivot
//
// Choose pivot of as a median of three random elements
// (or as a single random element for small ranges)
//
// (TODO: better/faster random number generation?)
//----------------------------------------------------------------------
template <typename RandomAccessIterator, typename Compare, 
	  typename KeyAccess>
RandomAccessIterator
random_pivot (RandomAccessIterator beg, RandomAccessIterator end,
	      Compare less, KeyAccess key)
{
  static double scale_to_one = 1.0L/(RAND_MAX+1.0);

  typedef typename std::iterator_traits<RandomAccessIterator>::difference_type 
    difference_type;
  difference_type size = std::distance(beg, end);
  double scale_to_size = size * scale_to_one;
  RandomAccessIterator pivot;
  pivot = beg + static_cast<difference_type>(scale_to_size*std::rand());
  if (size > 50) {
    RandomAccessIterator pivot2, pivot3;
    pivot2 = beg + static_cast<difference_type>(scale_to_size*std::rand());
    pivot3 = beg + static_cast<difference_type>(scale_to_size*std::rand());
    pivot = iter_median3(pivot, pivot2, pivot3, less, key);
  }
  return pivot;
}
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

template <typename RandomAccessIterator, typename Compare>
RandomAccessIterator
random_pivot(RandomAccessIterator beg, RandomAccessIterator end,
		    Compare less)
{
  return random_pivot(beg, end, less,
		      default_key<RandomAccessIterator>());
}

template <typename RandomAccessIterator>
RandomAccessIterator
random_pivot(RandomAccessIterator beg, RandomAccessIterator end)
{
  typedef typename std::iterator_traits<RandomAccessIterator>::value_type
    key_type;
  return random_pivot(beg, end, std::less<key_type>(),
		      default_key<RandomAccessIterator>());
}




//======================================================================
// Partitioning
//======================================================================

//----------------------------------------------------------------------
// function binary_partition
//
// Partitions the input range [beg,end) into two parts [beg,mid) and
// [mid,end) with all elements smaller than pivot in [beg, mid) and
// all elements larger than pivot in [mid,end). Equal elements may
// be in either subrange. Returns mid.
// PRECONDITION: pivot is in the range [beg,end)
//----------------------------------------------------------------------
template <typename RandomAccessIterator, typename Compare, 
	  typename KeyAccess>
RandomAccessIterator
binary_partition(RandomAccessIterator beg, RandomAccessIterator end,
		 RandomAccessIterator pivot, Compare less, KeyAccess key)
{
  if (beg == --end) return beg; 

  // move the pivot to a safe place and setup sentinels
  std::swap(*beg, *pivot);
  if (less(key(end), key(beg))) {
    std::swap(*beg, *end);
    pivot = end;
  } else {
    pivot = beg;
  }

  // now the pivot won't move anymore
  typename KeyAccess::result_type pivot_key = key(pivot);

  while (true) {
    do { ++beg; } while (less(key(beg), pivot_key));
    do { --end; } while (less(pivot_key, key(end)));
    if (!(beg < end)) return beg;
    std::swap(*beg, *end);
  }
}
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

template <typename RandomAccessIterator, typename Compare>
RandomAccessIterator
binary_partition(RandomAccessIterator beg, RandomAccessIterator end,
		 RandomAccessIterator pivot, Compare less)
{
  return binary_partition(beg, end, pivot, less, 
			  default_key<RandomAccessIterator>());
}

template <typename RandomAccessIterator>
RandomAccessIterator
binary_partition(RandomAccessIterator beg, RandomAccessIterator end,
		 RandomAccessIterator pivot)
{
  typedef typename std::iterator_traits<RandomAccessIterator>::value_type
    key_type;
  return binary_partition(beg, end, pivot, std::less<key_type>(), 
			  default_key<RandomAccessIterator>());
}


//----------------------------------------------------------------------
// function ternary_partition
//
// Partitions the input range [beg,end) into three parts according to
// comparisons to the pivot:  | less | equal | greater |
// Returns the middle (equal) range.
// PRECONDITION: pivot is in the range [beg,end)
//
// METHOD
// The algorithm maintains the following invariant:
//
// | equal  |  less  |  unknown  | greater | equal |
//  ^        ^        ^         ^         ^         ^
// beg  first_less   left     right  last_greater  end
//
//----------------------------------------------------------------------
template <typename RandomAccessIterator, typename Compare,
	  typename KeyAccess>
std::pair<RandomAccessIterator,RandomAccessIterator>
ternary_partition(RandomAccessIterator beg, RandomAccessIterator end,
		  RandomAccessIterator pivot, Compare less,
		  KeyAccess key)
{
  RandomAccessIterator left = beg, right = end-1;
  RandomAccessIterator first_less = left, last_greater = right;

  // first move pivot to a safe place and setup sentinels
  std::iter_swap(pivot, left);
  if (less(key(right), key(left))) {
    std::iter_swap(left, right);
    pivot = right;
    ++left; --right; --last_greater;
  } else {
    pivot = left;
    ++left; ++first_less;
  }

  // now the pivot won't move anymore
  typename KeyAccess::result_type pivot_key = key(pivot);

  while(true) {
    while(less(key(left), pivot_key)) ++left;
    while (!less(pivot_key, key(left)) && left < right) {
      std::iter_swap(left++, first_less++);
      while(less(key(left), pivot_key)) ++left;
    }
    while(less(pivot_key, key(right))) --right;
    while (!less(key(right), pivot_key) && left < right) {
      std::iter_swap(right--, last_greater--);
      while(less(pivot_key, key(right))) --right;
    }
    if (!(left < right)) break;
    std::iter_swap(left++, right--);
  }

  // now the situation is this:
  //
  // | equal  |     less     |    greater   | equal |
  //  ^        ^            ^ ^            ^         ^
  // beg  first_less    right left    last_greater  end
  //
  // or this:
  //
  // | equal  |     less    |=|   greater   | equal |
  //  ^        ^             ^             ^         ^
  // beg  first_less       left       last_greater  end
  //                       right
  //
  // next: swap equal parts to the middle

  typename std::iterator_traits<RandomAccessIterator>::difference_type size;
  size = std::min(first_less-beg, left-first_less);
  std::swap_ranges(beg, beg+size, left-size);
  size = std::min(end-1-last_greater, last_greater-right);
  std::swap_ranges(right+1, right+1+size, end-size);
  return std::make_pair(beg+(left-first_less), end-(last_greater-right));
}
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

template <typename RandomAccessIterator, typename Compare>
std::pair<RandomAccessIterator,RandomAccessIterator>
ternary_partition(RandomAccessIterator beg, RandomAccessIterator end,
		  RandomAccessIterator pivot, Compare less)
{
  return ternary_partition(beg, end, pivot, less,
			   default_key<RandomAccessIterator>());
}

template <typename RandomAccessIterator>
std::pair<RandomAccessIterator,RandomAccessIterator>
ternary_partition(RandomAccessIterator beg, RandomAccessIterator end,
		  RandomAccessIterator pivot)
{
  typedef typename std::iterator_traits<RandomAccessIterator>::value_type
    key_type;
  return ternary_partition(beg, end, pivot, std::less<key_type>(),
			   default_key<RandomAccessIterator>());
}


//======================================================================
// Sorting
//======================================================================

//----------------------------------------------------------------------
// function insertion_sort
//
//----------------------------------------------------------------------
template <typename RandomAccessIterator, typename Compare, 
	  typename KeyAccess>
void
insertion_sort(RandomAccessIterator beg, RandomAccessIterator end, 
	       Compare less, KeyAccess key)
{
  if (end-beg <= 1) return;
  for (RandomAccessIterator i = beg+1; i != end; ++i) {
    typename KeyAccess::result_type i_key = key(i); // i is not moved ...
    RandomAccessIterator j = i-1;
    if (!less(i_key, key(j))) continue;
    for ( ; j != beg; --j) {
      if (!less(i_key, key(j-1))) break;
      std::swap(*(j-1), *j);
    }
    std::swap(*i, *j);                              // ... until here
  }
}
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

template <typename RandomAccessIterator, typename Compare>
void
insertion_sort(RandomAccessIterator beg, RandomAccessIterator end, 
	       Compare less)
{
  insertion_sort(beg, end, less, default_key<RandomAccessIterator>());
}

template <typename RandomAccessIterator>
void
insertion_sort(RandomAccessIterator beg, RandomAccessIterator end)
{
  typedef typename std::iterator_traits<RandomAccessIterator>::value_type
    key_type;
  insertion_sort(beg, end, std::less<key_type>(), 
		 default_key<RandomAccessIterator>());
}


//----------------------------------------------------------------------
// function quicksort
//
//----------------------------------------------------------------------
template <typename RandomAccessIterator, typename Compare, 
	  typename KeyAccess>
void
quicksort(RandomAccessIterator beg, RandomAccessIterator end, 
	  Compare less, KeyAccess key)
{
  while (end-beg > 15) {
    RandomAccessIterator pivot = random_pivot(beg, end, less, key);
    RandomAccessIterator mid = binary_partition(beg, end, pivot, less, key);
    quicksort(beg, mid, less, key);
    beg = mid;
  }
  insertion_sort(beg, end, less, key);
}
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

template <typename RandomAccessIterator, typename Compare>
void
quicksort(RandomAccessIterator beg, RandomAccessIterator end, Compare less)
{
  quicksort(beg, end, less, default_key<RandomAccessIterator>());
}

template <typename RandomAccessIterator>
void
quicksort(RandomAccessIterator beg, RandomAccessIterator end)
{
  typedef typename std::iterator_traits<RandomAccessIterator>::value_type
    key_type;
  quicksort(beg, end, std::less<key_type>(), 
	    default_key<RandomAccessIterator>());
}


#endif // PARTITION_HPP
