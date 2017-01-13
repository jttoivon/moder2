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
#include "bndm.hpp"
#include "common.hpp"
#include "iupac.hpp"

#include <string>
#include <cstdlib>
#include <cstring>
#include <cassert>



// return number of occurences
// Pattern characters are 'ACGT'
// Joker characters are uipac codes
// N denotes anything
int
BNDM_with_joker(const std::string& text, const std::string& pattern) 
{
  const int ASIZE=256;
  const int WORD_SIZE=64;

  int counter=0;  // number of occurences
  
  int pattern_length = pattern.length();
  int text_length = text.length();
  std::vector<long long int> B(ASIZE, 0);
  long long int s, d;
  int i, pos, last;
  assert(pattern_length <= WORD_SIZE);

  /* Pre processing */
  s=1;
  for (i=pattern_length-1; i>=0; i--){  // iterate through all positions
    for (const char* pos=iupac_class(pattern[i]); *pos != '\0'; ++pos)
      B[(unsigned char)*pos] |= s;
    s <<= 1;
  }

  /* Searching phase */
  pos=0;
  while (pos <= text_length - pattern_length){
    i=pattern_length-1; last=pattern_length;
    d = ~0;
    while (i>=0 && d!=0) {
      d &= B[(unsigned char)text[pos+i]];
      i--;
      if (d != 0){
	if (i >= 0)
	  last = i+1;
	else
	  ++counter;
      }
      d <<= 1;
    }
    pos += last;
  }

  return counter;
}

