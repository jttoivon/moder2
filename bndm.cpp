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

