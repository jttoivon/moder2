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
#ifndef IUPAC_HPP
#define IUPAC_HPP

#include <string>
#include <cassert>
#include <cstring>
#include <vector>

typedef std::vector<double> dvector;

extern std::string iupac_chars;

class iupac_class_type
{
public:

  iupac_class_type()
  {
    for (int i=0; i < 256; ++i)
      char_to_class[i] = 0;

    int size = sizeof(char_classes)/sizeof(char_class_t);
    for (int i=0; i < size; ++i) {
      char_to_class[(unsigned char)char_classes[i].c] = char_classes[i].str;

      char_to_bits_[(unsigned char)char_bits[i].c] = char_bits[i].bits;
      bits_to_char_[char_bits[i].bits] = char_bits[i].c;
    }
  }

  int
  char_to_bits(char c) {
    return char_to_bits_[(unsigned char)c];
  }
  
  char
  bits_to_char(int i) {
    return bits_to_char_[i];
  }

  bool
  is_iupac_code(char c) const
  {
    return char_to_class[(unsigned char)c] != 0;
  }

  // function from iupac code to the corresponding subset of {A,C,G,T}
  const char*
  operator()(char c)
  {
    assert(char_to_class[(unsigned char)c] != 0);
    
    return char_to_class[(unsigned char)c];
  }

private:
  const char* char_to_class[256];  // If null then character is not an iupac character

  int char_to_bits_[256];
  char bits_to_char_[16];
  
  typedef struct {char c; const char* str;} char_class_t;
  typedef struct {char c; int bits;}             char_bits_t;

  static char_class_t char_classes[16];
  static char_bits_t  char_bits[16];
};


static iupac_class_type iupac_class;

dvector
iupac_probability(char c);

bool
iupac_match(char c, char char_class);

std::string
complement_set(char char_class);


bool
is_iupac_string(const std::string& str);

bool
iupac_string_match(const std::string& str, const std::string& pattern);

char complement(char c);

char complement_rna(char c);

// The most significant bit in the bit vector corresponds to position 0 in strings.
// If length of strings is n, then
// result & 1 == 1 iff str[n-1] does not match pattern[n-1]. 
template <typename T>
T
iupac_mismatch_positions(const std::string& str, const std::string& pattern)
{
  assert(str.length() == pattern.length());
  assert(str.length() * 2 <= sizeof(T)*8);
  T result = 0;
  for (int i=0; i < str.length(); ++i) {
    result <<= 1;
    result |= (iupac_match(str[i], pattern[i]) ? static_cast<T>(0) : static_cast<T>(1));
  }
  
  return result;
}

class sequences_of_iupac_string
{
public:
  sequences_of_iupac_string(const std::string& s)
    : k(s.length()), iupac_classes(k), base(k, 0), pattern(k, '-'), v(k, 0)
  {

    number_of_combinations = 1;
    for (int i=0; i < k; ++i) {
      iupac_classes[i] = iupac_class(s[i]);
      int size = strlen(iupac_classes[i]);
      base[i] = size - 1;
      number_of_combinations *= size;
    }
    reset();
  }

  int number_of_strings() const
  {
    return number_of_combinations;
  }
  
  void
  reset()
  {
    for (int i=0; i < k; ++i) {
      v[i] = 0;
      pattern[i] = iupac_classes[i][0];          // Initialize pattern to the first sequence of iupac string
    }
    v[k-1]=-1;  // Initialize
    current_string = 0;
  }

  // The loop goes through nucleotide sequences that belong to
  //given iupac sequence in alphabetical order
  std::string
  next()
  {
    if (current_string == number_of_combinations)
      return "";
    
    int i;
    for (i=k-1; v[i] == base[i]; --i) {
      v[i]=0;
      pattern[i] = iupac_classes[i][v[i]];
    }
    v[i]++;
    pattern[i] = iupac_classes[i][v[i]];

    ++current_string;
    return pattern;
  }
  
private:
  int k;
  int number_of_combinations;
  int current_string;
  std::vector<const char*> iupac_classes;
  std::vector<int> base;
  std::string pattern;
  std::vector<int> v;  // helper array
};

#endif // IUPAC_HPP
