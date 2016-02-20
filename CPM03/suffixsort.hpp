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
// This is the main file of the difference cover based suffix array
// construction algorithm.
//======================================================================

#ifndef SUFFIXSORT_HPP
#define SUFFIXSORT_HPP

#include "doubling.hpp"
#include "stringsort.hpp"
#include "checker.hpp"
#include "difference_cover.hpp"
#include "timing.hpp"

#include <algorithm>
#include <vector>
#include <string>
#include <iterator>
#include <functional>
#include <stdexcept>
#include <sstream>

//======================================================================
// The algorithm is controlled by three macros:
// DCOVER: Use difference cover modulo 2^DCOVER.
// DOUBLING_ONLY: Run doubling for the full suffix array after
//     stringsorting it to depth 2^DCOVER.  This is also done 
//     if DCOVER is 0 or 1, since then the sample would contain
//     all suffixes anyway.
// STRINGSORT_ONLY: Do only stringsorting to depth 2^DCOVER-1.
//     This is also done if the length of string is less than 2^DCOVER.
//======================================================================

#ifndef DCOVER
#define DCOVER 5
#endif

#if DCOVER < 0 || DCOVER > 30
#error DCOVER must be in the range [0,30]
#elif DCOVER <= 1
#define DOUBLING_ONLY
#endif


//----------------------------------------------------------------------
// functor mark_unsorted_group
//
// When string sorting leaves a group of suffixes unsorted, 
// it marks that group by calling this functor.
//
// The marking is done by one's complementing the first and last
// entries of the group. Complemented values are assumed to become
// negative. Complementing is used instead of negation to be able
// mark the zero.
//----------------------------------------------------------------------
class mark_unsorted_group {
public:
  template <typename SArrayIterator>
  void operator() (SArrayIterator beg, SArrayIterator end) const { 
    *beg = ~*beg; 
    --end;
    *end = ~*end;  // note: group of size one stays unmarked
  }
};

//----------------------------------------------------------------------
// functor do_nothing
//
// Replaces mark_unsorted_group when nothing needs to be done
// to unsorted groups (this happens when doing stringsort only)
//----------------------------------------------------------------------
class do_nothing {
public:
  template <typename Iterator>
  void operator() (Iterator, Iterator) {}
};


//----------------------------------------------------------------------
// functor less_using_issa
//
// Comparison functor for the final sorting stage that sorts
// remaining unsorted groups using the inverse sparse suffix array.
//----------------------------------------------------------------------
template <typename ISSAIterator>
class less_using_issa {
private:
  typedef typename std::iterator_traits<ISSAIterator>::value_type 
          offset_type;
  const ISSAIterator issa;
  const dc_sampler<offset_type>& sample;
public:
  less_using_issa(ISSAIterator beg, const dc_sampler<offset_type>& s) 
    : issa(beg), sample(s) {}
  template <typename StringPosType>
  bool operator() (StringPosType a, StringPosType b) const {
    std::pair<offset_type,offset_type> rep = sample.pair(a, b);
    return issa[rep.first] < issa[rep.second];
  }
};




//======================================================================
// function: construct_suffix_array
//
// The main suffix array construction algorithm.
//
// Constructs the suffix array of string [str,str_end) 
// into range [sa,sa_end).
//
// Construction consists of the following steps:
// 1. Choose a sample of suffixes (using a difference cover)
//    and sort them to limited depth using string sorting.
// 2. Complete the sorting of the sample using the doubling algorithm.
//    The *inverse* suffix array of the sample is the result retained.
// 3. Sort the set of all suffixes to limited depth using string sorting.
// 4. Complete the sorting using the inverse sparse suffix array.
//
// Special cases:
// DOUBLING_ONLY is defined (or DCOVER is 0 or 1):
//    Do only steps 1 and 2 with the "sample" of all suffixes.
// STRINGSORT_ONLY is defined (or string length < 2^DCOVER):
//    Do only step 3.
// 
// TYPE REQUIREMENTS
//  - both iterator types are models of Random Access Iterator
//  - StringIterator::value_type is a model of Strict Weakly Comparable
//  - SArrayIterator::value_type is an integral type capable of
//    representing all values in the range [-|string|-1,|string|]
// PRECONDITIONS
//  - both ranges are valid
//  - |string| = |suffix array|
// POSTCONDITION
//  see checker.hpp
//======================================================================

template <typename StringIterator, typename SArrayIterator>
void
construct_suffix_array(StringIterator str, StringIterator str_end,
		       SArrayIterator sa, SArrayIterator sa_end)
{
  typedef typename std::iterator_traits<SArrayIterator>::value_type 
    offset_type;
  offset_type size = str_end-str;
  if (sa_end-sa != size) {
    std::ostringstream os;
    os << "construct_suffix_array:"
       << " suffix array size " << sa_end-sa
       << " is different from string length " << str_end-str;
    throw std::logic_error(os.str());
  }
  offset_type period = (1<<DCOVER);

#ifdef DOUBLING_ONLY
  bool doubling_only = true;
#else
  bool doubling_only = false;
#endif

#ifdef STRINGSORT_ONLY
  bool stringsort_only = true;
#else
  bool stringsort_only = false;
#endif

  bool very_short_string = (size < period);

  TIME_START2(total);
  TIME_START2(intermediate);

  if (stringsort_only || very_short_string) {

    //---------------------------------------------------------
    // Step 3 only: stringsort
    //---------------------------------------------------------

    TIME_CHECK2(intermediate); // Step 1pre
    TIME_CHECK2(intermediate); // Step 1
    TIME_CHECK2(intermediate); // Step 2pre
    TIME_CHECK2(intermediate); // Step 2

    SArrayIterator i;
    offset_type j;
    for(i = sa, j = 0; i<sa_end; ++i, ++j) { *i = j; }
    
#ifdef DEBUG
    std::cerr << "Sorting only to depth " << period-1 << endl;
#endif

    suffix_mkqsort_noempty(str, str_end, sa, sa_end, 
			   period-1, do_nothing());

    TIME_CHECK2(intermediate); // Step 3
    TIME_CHECK2(intermediate); // Step 4

  } else {
    
    //---------------------------------------------------------
    // Step 1: stringsort sample
    //---------------------------------------------------------
    
    // initialize the sample
    offset_type samplerange = size+1;  // Because of +1, the empty suffix
                                    // may be in the sample. In that case,
                                    // it is needed as a sentinel.
    int dcover = DCOVER;
    
    if (doubling_only) {
      --samplerange;    // No empty suffix
      dcover = 0;       // Sample that contains all suffixes
    }
    
    dc_sampler<offset_type> sample(samplerange, dcover);
    offset_type samplesize = sample.samplesize();
    SArrayIterator ssa = sa;                    // use suffix array to
    SArrayIterator ssa_end = sa + samplesize;   // store the sample
    SArrayIterator sample_end = sample.fill(ssa);
    if (sample_end != ssa_end) {
      std::ostringstream os;
      os << "construct_suffix_array:"
	 << " sample.fill produced " << ssa_end-ssa
	 << " samples but samplesize is " << str_end-str;
      throw std::logic_error(os.str());
    }

    TIME_CHECK2(intermediate); // Step 1pre
    
#ifdef DEBUG
    std::cerr << "Starting sample sort" << endl;
#endif
    
    // sort the sample to limited depth
    suffix_mkqsort(str, str_end, ssa, ssa_end, period,
		      mark_unsorted_group());

    TIME_CHECK2(intermediate); // Step 1

    //---------------------------------------------------------
    // Step 2: doubling sort sample
    //---------------------------------------------------------
   
#ifdef DEBUG
    std::cerr << "Initializing for doubling" << endl;
#endif
    
    // initialize inverse sparse suffix array
    typename std::vector<offset_type> issa(samplesize);
    typedef typename std::vector<offset_type>::iterator issa_iterator;
    
    SArrayIterator first = ssa;
    while (first != ssa_end) {
      if (*first >= 0) {
	SArrayIterator last;
	issa[sample.pack(*first)] = first - ssa;
	for (last=first+1; last != ssa_end && *last >= 0; ++last) {
	  issa[sample.pack(*last)] = last - ssa;
	}
	*first = -(last-first);
	first = last;
      } else {
	SArrayIterator last;
	*first = ~*first;
	for (last=first+1; *last >= 0; ++last) ;
	*last = ~*last;
	for ( ; first <= last; ++first) {
	  *first = sample.pack(*first);
	  issa[*first] = last - ssa;
	}
      }
    }

    TIME_CHECK2(intermediate); // Step 2pre
    
#ifdef DEBUG
    std::cerr << "Starting doubling" << endl;
#endif
    
    // fully sort the sample by doubling algorithm
    doubling_sort(ssa, ssa_end, issa.begin(), issa.end(), 
		  sample.packed_period());

    TIME_CHECK2(intermediate); // Step 2
    

    if (!doubling_only) {
      
      //---------------------------------------------------------
      // Step 3: stringsort full suffix array
      //---------------------------------------------------------
      
      // initialize full suffix array
      // (sparse suffix array is not needed anymore, only its inverse)
      SArrayIterator i;
      offset_type j;
      for(i = sa, j = 0; i!=sa_end; ++i, ++j) { *i = j; }
      
#ifdef DEBUG
      std::cerr << "Starting full suffix array sort" << endl;
#endif
      
      suffix_mkqsort_noempty(str, str_end, sa, sa_end, period-1,
			     mark_unsorted_group());

      TIME_CHECK2(intermediate); // Step 3
    
      //------------------------------------------------------------
      // Step 4: final sort using inverse suffix array of the sample
      //------------------------------------------------------------
      
#ifdef DEBUG
      std::cerr << "Starting final sort" << endl;
#endif
      
      SArrayIterator first = sa;
      while (true) {
	while (first != sa_end && *first >= 0) ++first;
	if (first == sa_end) break;
	*first = ~*first;
	SArrayIterator last = first + 1;
	while (*last >= 0) ++last;
	*last = ~*last;
	++last;
	quicksort(first, last, 
		  less_using_issa<issa_iterator>(issa.begin(), sample));
	first = last;
      }
      
      TIME_CHECK2(intermediate); // Step 4

    } else { // doubling_only
      
      // compute the suffix array from the inverse suffix array
      issa_iterator i;
      for (i = issa.begin(); i != issa.end(); ++i) {
      	sa[*i] = i - issa.begin();
      }
      
      TIME_CHECK2(intermediate); // Step 3
      TIME_CHECK2(intermediate); // Step 4
    
    }

  } // else branch of (stringsort_only || very_short_string)

  TIME_CHECK2(total);

#ifndef NO_CHECKING
  if (!stringsort_only || very_short_string) {
    check_suffix_array(str, str_end, sa, sa_end);
  }
#endif

  TIME_CHECK2(total);
  
}

#endif // SUFFIXSORT_HPP
