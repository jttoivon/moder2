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
// This file implements the computation of difference covers
//======================================================================

#include <algorithm>
#include <numeric>
#include <stdexcept>
#include <iostream>

#include "difference_cover.hpp"


//----------------------------------------------------------------------
// Precomputed difference covers. 
//
// The covers from 0 to 6 are from
//   Luk, Wong: Two new quorum base algorithms for distributed mutual
//   exclusion. In Proc. 17th International Conference on Distributed 
//   Computing Systems, pp. 100-106. IEEE, 1997.
// Covers 7 and 8 were found using an exhaustive search algorithm.
// All except cover 8 are known to be optimal.
//----------------------------------------------------------------------
namespace {
  const unsigned max_precomputed_cover = 8;
  const unsigned coversizes[max_precomputed_cover+1] = 
    {1,2,3,4,5,7,9,13,20};
  const unsigned cover0[] = {0};
  const unsigned cover1[] = {0,1};
  const unsigned cover2[] = {0,1,2};
  const unsigned cover3[] = {0,1,2,4};
  const unsigned cover4[] = {0,1,2,5,8};
  const unsigned cover5[] = {0,1,2,3,7,11,19};   //{0,7,8,10,14,19,23};
  const unsigned cover6[] = {0,1,2,5,14,16,34,42,59};
  const unsigned cover7[] = {0,1,3,7,17,40,55,64,75,85,104,109,117};      
  const unsigned cover8[] = {0,1,3,7,12,20,30,44,65,80,89,96,114,122,
			     128,150,196,197,201,219}; 
  const unsigned * covers[] = { cover0, cover1, cover2, cover3, cover4, 
				cover5, cover6, cover7, cover8 };
}


//----------------------------------------------------------------------
// constructor of difference_cover
//
// Constructs a difference cover modulo a power of two.
// (see difference_cover.hpp)
//
// Small difference covers have been precomputed (see above).
// Larger difference covers are computed by the algorithm in
//   C.J. Colbourn and A.C.H. Ling. Quorums from Difference Covers.
//   Information Processing Letters 75(1-2):9-12, 2000
//----------------------------------------------------------------------
difference_cover::difference_cover(unsigned lv)
  : logmod(lv), mask((1<<lv)-1),
    coverrank(1<<lv), diff2pos(1<<lv,1<<lv)
{
  if (logmod <= max_precomputed_cover) {
    
    // use precomputed covers
    coversize = coversizes[logmod]; 
    cover.assign(covers[logmod], covers[logmod] + coversize);
    
  } else {
    
    // use the Colbourn-Ling algorithm
    int v = (1<<logmod);
    int r = 0;
    while (24*r*r+36*r+13 < v) ++r;
    coversize = 6*r+4;
    cover.resize(coversize);
    std::vector<unsigned>::iterator i = cover.begin();
    *i = 0; ++i;
    std::fill(i, i+r,         1); i+=r;
    *i = r+1; ++i;
    std::fill(i, i+r,     2*r+1); i+=r;
    std::fill(i, i+2*r+1, 4*r+3); i+=2*r+1;
    std::fill(i, i+r+1,   2*r+2); i+=r+1;
    std::fill(i, i+r,         1); i+=r;
    if (i != cover.end()) {
      throw std::logic_error("error in the Coulbourn-Ling algorithm");
    }
    std::partial_sum(cover.begin(), cover.end(), cover.begin());
    
  }
  
#ifdef DEBUG
  std::cerr << "difference cover modulo " << (1<<logmod);
  std::cerr << " of size " << coversize << ": " << std::endl;
  std::copy(cover.begin(), cover.end(), 
	    std::ostream_iterator<unsigned>(std::cerr, " ") );
  std::cerr << std::endl;
#endif
  
  // compute the lookup tables
  unsigned i, j, modulus=1<<logmod;
  for (i=0, j=0; i<modulus; ++i) {
    coverrank[i] = j;
    if (j<coversize && cover[j] <= i) { ++j; }
  }
  for (i=coversize-1; i<coversize; --i) {
    for (j=0; j<coversize; ++j) {
      diff2pos[(cover[j]-cover[i])&mask] = cover[i];
    }
  }
  
  // check that this is a difference cover
  for (i=0; i<modulus; ++i) {
    if (diff2pos[i] == modulus) {
      throw std::logic_error("bad difference cover");
    }
#if DEBUG > 1
    std::cerr << i << " = " << ((diff2pos[i]+i)&mask);
    std::cerr << " - " << diff2pos[i];
    std::cerr << " (mod " << modulus << ")";
    std::cerr << std::endl;
#endif
  }
}
