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
#include "iupac.hpp"
#include "data.hpp"
#include <cstring>
#include <algorithm>

std::string iupac_chars = "ACGTWSMKRYBDHVN";

bool
is_iupac_string(const std::string& str)
{
  for (int i=0; i < str.length(); ++i) {
    if (not iupac_class.is_iupac_code(str[i]))
      return false;
  }
  return true;
}


// c is an iupac code
// char_class is iupac code
// is c \subset char_class
bool
iupac_match(char c, char char_class)
{
  //std::string nucs = "ACGT";
  //assert(std::count(iupac_chars.begin(), iupac_chars.end(), c));  // Make sure that c is a nucleotide: A, C, G, or T
  assert(iupac_class.is_iupac_code(c));
  assert(iupac_class.is_iupac_code(char_class));

  int cbits = iupac_class.char_to_bits(c);
  if (cbits == (cbits & iupac_class.char_to_bits(char_class)))
    return true;
  else
    return false;
  /*
  const char* p = iupac_class(char_class);
  for (; *p != '\0'; ++p)
    if (*p == c)                               // Is c contained in the char_class
      return true;
  return false;
  */
}

// str is nucleotide sequence
// pattern is iupac sequence
bool
iupac_string_match(const std::string& str, const std::string& pattern)
{
  assert(str.length() == pattern.length());
  for (int i=0; i < str.length(); ++i)
    if (not iupac_match(str[i], pattern[i]))
      return false;
  return true;
}

std::string
complement_set(char char_class)
{
  std::string result;
  std::string nucs("ACGT");
  for (int i=0; i < 4; ++i) {
    if (not iupac_match(nucs[i], char_class))
      result.push_back(nucs[i]);
  }
  return result;
}

/*
bool
iupac_string_match(const char* str, const char* pattern)
{
  return iupac_string_match(std::string(str), std::string(pattern));
}
*/

iupac_class_type::char_class_t iupac_class_type::char_classes[15] =
  {{'A', "A"},
   {'C', "C"},
   {'G', "G"},
   {'T', "T"},
   {'W', "AT"},
   {'S', "CG"},
   {'M', "AC"},
   {'K', "GT"},
   {'R', "AG"},
   {'Y', "CT"},
   {'B', "CGT"},
   {'D', "AGT"},
   {'H', "ACT"},
   {'V', "ACG"},
   {'N', "ACGT"}};

iupac_class_type::char_bits_t iupac_class_type::char_bits[15] =
  {{'A', 8      },
   {'C',   4    },
   {'G',     2  },
   {'T',       1},
   {'W', 8    +1},
   {'S',   4+2  },
   {'M', 8+4    },
   {'K',     2+1},
   {'R', 8  +2  },
   {'Y',   4  +1},
   {'B',   4+2+1},
   {'D', 8  +2+1},
   {'H', 8+4  +1},
   {'V', 8+4+2  },
   {'N', 8+4+2+1}};


char complement(char c) 
{
  switch(c) {
  case 'A': return 'T';
  case 'C': return 'G';
  case 'G': return 'C';
  case 'T': return 'A';

  case 'W': return 'W';
  case 'S': return 'S';
  case 'M': return 'K';
  case 'K': return 'M';
  case 'R': return 'Y';
  case 'Y': return 'R';

  case 'B': return 'V';
  case 'D': return 'H';
  case 'H': return 'D';
  case 'V': return 'B';

  case 'N': return 'N';
    
  default:  return c;
  }
}

// Return a probability distribution.
// Probability is distributed evenly between nucleotides in the iupac class 'c'
dvector
iupac_probability(char c)
{
  dvector result(4, 0.0);
  const char* str = iupac_class(c);
  assert(str != 0);
  int size = strlen(str);
  assert(0 < size and size <= 4);

  double a = 1.0 / size;
  for (int i=0; i < size; ++i)
    result[to_int(str[i])] = a;

  return result;
}
