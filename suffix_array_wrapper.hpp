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
#ifndef SUFFIX_ARRAY_WRAPPER_HPP
#define SUFFIX_ARRAY_WRAPPER_HPP

#include "CPM03/suffixsort.hpp"
#ifndef TIMING
#define TIMING
#endif
#include "timing.hpp"

#include <string>
#include <cassert>
#include <cstring>
#include <cstdio>


class suffix_array
{
public:
  typedef std::vector<unsigned char>::difference_type offset_type;
  typedef unsigned char uchar;

  suffix_array(bool use_rna=false) : text(""), sa(0)
  {
    init_array(use_rna);
  }

  suffix_array(const std::string& str, bool use_rna=false) : text(str), sa(str.length())
  {
    TIME_START(t);
    construct_suffix_array(text.begin(), text.end(),
			   sa.begin(), sa.end());
    TIME_PRINT("Constructing suffix array took %.2f seconds\n", t);
    init_array(use_rna);

    /*    
    int n=text.length();
    for (int i=0; i < n-1; ++i)
      assert(text.substr(sa[i]) < text.substr(sa[i+1]));
    */
  }

  void
  reset(const std::string& str) 
  {
    TIME_START(t);
    text = str;
    sa.resize(str.length());

    construct_suffix_array(text.begin(), text.end(),
			   sa.begin(), sa.end());
    TIME_PRINT("Constructing suffix array took %.2f seconds\n", t);

    /*    
    int n=text.length();
    for (int i=0; i < n-1; ++i)
      assert(text.substr(sa[i]) < text.substr(sa[i+1]));
    */
  }


  int
  count(const std::string& pattern) const
  {
    std::vector<offset_type> dummy;
    return count_helper(pattern, dummy, false);
  }

  int 
  locate(const std::string& pattern, std::vector<offset_type>& positions)
  {
    return count_helper(pattern, positions, true);
  }


  int 
  count_iupac(const std::string& pattern) const
  {
    std::vector<offset_type> dummy;
    return count_iupac_helper(pattern, 0, 0, text.length(), dummy, false);
  }

  int 
  locate_iupac(const std::string& pattern, std::vector<offset_type>& positions) const
  {
    assert(positions.size() == 0);
    return count_iupac_helper(pattern, 0, 0, text.length(), positions, true);
  }


  ~suffix_array() {
    
  }



private:

  offset_type 
  GetPos(offset_type position) const
  {
    return sa[position];
  }

  int ComparePos(offset_type position, const std::string& pattern) const
  {
    offset_type n=text.length();
    offset_type m=pattern.length();
    offset_type j=0;
    
    while (position+j < n && j < m && text[position+j]==pattern[j]) {
      j++;
    }
    if (j >= m)
      return 0;
    else if (position+j >= n)
      return +1;
    else if (text[position+j]>pattern[j])
      return -1;
    else
      return +1;
  }




/*
 -1 text >  pattern
  0 text == pattern
 +1 text <  pattern
*/
  int 
  ComparePos2(offset_type position, uchar pchar, int j) const
  {
    offset_type n=text.length();

    if (position+j >= n) /* text piece is shorter than pattern */
      return +1;
    else if (text[position+j] > pchar)
      return -1;
    else if (text[position+j] == pchar)
      return 0;
    else
      return +1;
  }

  typedef struct
  {
    char c;
    const char* str;
  } char_classes_t;


  static const char_classes_t char_classes_dna[16];
  static const char_classes_t char_classes_rna[16];


  const char* char_classes[256];

  void 
  init_array(bool use_rna = false)
  {
    if (use_rna) {
      int size = sizeof(char_classes_rna)/sizeof(char_classes_t);
      for(int i=0; i < size; ++i)
	char_classes[(unsigned char)char_classes_rna[i].c] = char_classes_rna[i].str;
    } else {
      int size = sizeof(char_classes_dna)/sizeof(char_classes_t);
      for(int i=0; i < size; ++i)
	char_classes[(unsigned char)char_classes_dna[i].c] = char_classes_dna[i].str;
    }
  }

  int 
  count_helper(const std::string& pattern, std::vector<offset_type>& positions, bool get_positions) const
  {
    if (text.length() == 0)
      return 0;
    offset_type L=0;
    offset_type R=text.length();
    offset_type Lraja;

    /* Search for the left boundary */
    if (ComparePos(GetPos(0), pattern) <= 0)  /* text is greater than or equal to pattern */
      Lraja=0;
    else {
      while (L<R-1) {  // while not consequent
        offset_type i = (R+L)/2;
        if (ComparePos(GetPos(i), pattern) <= 0) /* text is greater than or equal to pattern */
	  R = i;
        else
	  L = i;
      }
      Lraja=R;
    }
    // Post-condition: Lraja is the first suffix that is greater than or equal to the pattern, or
    //                 if not such suffix exists, Lraja is text.length()


    /* Search for the right boundary */
    L=0;
    R=text.length();
    offset_type Rraja;
    if (ComparePos(GetPos(0), pattern) < 0) /* text is greater than pattern */
      Rraja=0;
    else {
      while (L < R-1) {
	offset_type i = (R+L)/2;
	if (ComparePos(GetPos(i), pattern) >= 0)  /* text is less than or equal to pattern */
	  L = i;
	else
	  R = i;
      }
      Rraja=R;
    }
    // Post-condition: Rraja is the first suffix that is greater than the pattern, or
    //                 if not such suffix exists, Rraja is text.length()

    //printf("L is %lu, R is %lu\n", Lraja, Rraja); 

    if (get_positions) {
      for (int i=Lraja; i < Rraja; ++i)
	positions.push_back(GetPos(i));
    }

    return Rraja-Lraja;

  }


  /* Input range: [L,R) and [a,b), half-open interval */
  /* Output range: [Lraja,Rraja],  closed interval    */
  int count_iupac_helper(const std::string& pattern, offset_type pos, offset_type a, offset_type b, 
			 std::vector<offset_type>& positions, bool get_positions) const
  {
    if (a == b)
      return 0;

    int old_count = positions.size();

    /* offset_type L=0; */
    /* offset_type R=_index->n; */
    int length = pattern.length();

    uchar pchar = pattern[pos];
    int result = 0;
    int number_of_chars = strlen(char_classes[pchar]);

    for (int r = 0; r < number_of_chars; ++r) {
      uchar current_char = char_classes[pchar][r];


      /* Search for the left boundary */

      offset_type L=a;
      offset_type R=b;

      offset_type Lraja;
      if (ComparePos2(GetPos(L), current_char, pos) <= 0)   /* text is greater than or equal to pattern */
	Lraja=L;
      else {
	while (L<R-1) { // while not adjacent
	  offset_type i = (R+L)/2;
	  if (ComparePos2(GetPos(i), current_char, pos) <= 0) /* text is greater than or equal to pattern */
	    R = i;
	  else
	    L = i;
	}
	Lraja=R;
      }


      /* Search for the right boundary */
      offset_type Rraja;
      L=a;
      R=b;
      if (ComparePos2(GetPos(L), current_char, pos) < 0) /* text is greater than pattern */
	Rraja=L;
      else {
	while (L < R-1) { // while not adjacent
	  offset_type i=(L+R)/2;
	  if (ComparePos2(GetPos(i), current_char, pos) >= 0)  /* text is smaller or equal to pattern */
	    L = i;
	  else
	    R = i;                                           /* text greater than pattern */
	}
	Rraja=R;
      }

      if (pos + 1 == length) {
	result += Rraja-Lraja;
	if (get_positions) {
	  for (int i=Lraja; i < Rraja; ++i)
	    positions.push_back(GetPos(i));
	}
      }
      else
	result += count_iupac_helper(pattern, pos+1, Lraja, Rraja, positions, get_positions);
    
      /*a = Rraja + 1;*/
    }


    if (get_positions)
      assert(result == positions.size() - old_count);

    return result;
  }


  suffix_array(const suffix_array& dummy) {}  // disable
  std::string text;
  std::vector<offset_type> sa;  // suffix array
};


#endif //SUFFIX_ARRAY_WRAPPER_HPP
