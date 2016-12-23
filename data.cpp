/*

    MODER is a program to learn DNA binding motifs from SELEX datasets.
    Copyright (C) 2016  Jarkko Toivonen

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
#include "data.hpp"
#include "common.hpp"
#include "probabilities.hpp"

//#include "common.hpp"

#include <cassert>
#include <iostream>
#include <algorithm>
#include <map>
#include <set>
#include <cstdio>

#include <boost/foreach.hpp>

#ifdef __int128
typedef __int128 myint128;
#else
typedef __int128_t myint128;           // for old gcc compilers
#endif


// interpret string of nucleid acids as a number in base 4
// The first nucleotide of the string is in the most significant end
big_int
dna_to_number(const std::string& s)  
{
  assert(s.length() <= (sizeof(big_int)*4));

  big_int sum=0;
  for (int i=0; i < s.length(); ++i)
    sum = sum*4 + to_int(s[i]);
  
  return sum;
}

// converts integer back to a nucled acid sequence
// Sequence has l nucleotides
// The first nucleotide of the string is in the most significant end
std::string
number_to_dna(big_int i, int l)  
{
  assert(l <= sizeof(big_int)*4);
  assert(i>=0);

  char nucs[] = "ACGT";

  std::string s(l, 'X');
  while (l-- > 0) {
    s[l] = nucs[i & 3];
    i = i / 4;
  }
  return s;
}

// interpret string of nucleid acids A,C,G,T,N as a number in base 5
// The first nucleotide of the string is in the most significant end
big_int
extended_dna_to_number(const std::string& s)  
{
  
  assert(s.size() < 28);

  big_int sum=0;
  for (int i=0;i<s.size(); ++i)
    sum = sum*5 + to_int(s[i]);
  
  if (sum < 0) {
    std::cout << s << std::endl;
    throw "sika";
  }
  return sum;
}

// converts integer back to a nucled acid sequence including N
// Sequence has l nucleotides
// The first nucleotide of the string is in the most significant end
std::string
number_to_extended_dna(big_int i, int l)  
{
  assert(i>=0);

  char nas[] = "ACGTN";

  std::string s;
  while (l-- > 0) {
    int m = i % 5;
    s.push_back(nas[m]);
    i = i / 5;
  }
  reverse(s.begin(), s.end());
  return s;
}
















// interpret string of nucleid acids as a number in base 4
// The first nucleotide of the string is in the most significant end
big_int
vector_to_number(const std::vector<int>& s)  
{
  
  assert(s.size() < 32);

  big_int sum=0;
  for (int i=0;i<s.size(); ++i)
    sum = sum*4 + s[i];

  
  if (sum < 0) {
    throw "sika";
  }
  return sum;
}

// converts integer back to a nucled acid sequence
// Sequence has l nucleotides
std::vector<int>
number_to_vector(big_int i, int l)  
{
  assert(i>=0);

  std::vector<int> s;
  while (l-- > 0) {
    int m = i % 4;
    s.push_back(m);
    i = i / 4;
  }
  reverse(s.begin(), s.end());
  return s;
}

// interpret string of nucleid acids A,C,G,T,N as a number in base 5
// The first nucleotide of the string is in the most significant end
big_int
extended_vector_to_number(const std::vector<int>& s)  
{
  
  assert(s.size() < 28);

  big_int sum=0;
  for (int i=0;i<s.size(); ++i)
    sum = sum*5 + s[i];
  
  if (sum < 0) {
    throw "sika";
  }
  return sum;
}

// converts integer back to a nucled acid sequence including N
// Sequence has l nucleotides
std::vector<int>
number_to_extended_vector(big_int i, int l)  
{
  assert(i>=0);

  std::vector<int> s;
  while (l-- > 0) {
    int m = i % 5;
    s.push_back(m);
    i = i / 5;
  }
  reverse(s.begin(), s.end());
  return s;
}

int
N_count(const std::string& s)
{
  return std::count(s.begin(), s.end(), 'N');
}

int
non_N_count(const std::string& s)
{
  return s.length() - N_count(s);
}

bool
has_atmost_one_gap(const std::string& s)
{
  if (s.size() == 0)
    return true;

  if (toupper(s[0]) == 'N' or toupper(s[s.size()-1]) == 'N')
    throw "String contains n in the border";

  bool in_gap = false;
  int number_of_switches = 0;
  for (int i=0; i < s.size(); ++i) {
    if ((toupper(s[i]) == 'N') != in_gap) {
      ++number_of_switches;
      in_gap = !in_gap;
    }
  }

  return number_of_switches <= 2;
}

std::vector<std::string>
remove_duplicate_reads(const std::vector<std::string>& sequences_orig)
{
  std::map<std::string, int> read_count;
  BOOST_FOREACH(std::string s, sequences_orig) {
    std::string key = std::min(s, reverse_complement(s));
    ++read_count[key];
  }

  std::set<std::string> duplicates;
  std::string s;
  int count;
  BOOST_FOREACH(boost::tie(s, count), read_count) {
    if (count > 1)
      duplicates.insert(s);
  }
  printf("Dataset contains %zu duplicates\n", duplicates.size());
  std::vector<std::string> sequences;
  BOOST_FOREACH(std::string s, sequences_orig) {
    // if (duplicates.count(s) or duplicates.count(reverse_complement(s)))
    //   continue;

    std::string sr = reverse_complement(s);
    BOOST_FOREACH(std::string d, duplicates) {
      if (hamming_distance(s, d) <= 1 or hamming_distance(sr, d) <= 1)
	goto end;
    }
    
    sequences.push_back(s);
  end:
    ;
  }

  return sequences;
}


myint128
get_even_mask()
{
  myint128 x=0;
  int number_of_bits = sizeof(myint128) * 8;
  for (int i=0; i < number_of_bits / 2; ++i) {
    x <<= 2;
    x += 1;
  }
  return x;
}

myint128 evenmask = get_even_mask();
myint128 oddmask = evenmask << 1;

// Returns
//  0 if hd is zero
//  1 if hd is one
//  2 if hd is 2 or greater
int
hamming_distance_with_bits(myint128 x, myint128 y)
{
  myint128 w = x ^ y;
  myint128 result = ((w & oddmask) >> 1) | (w & evenmask);
  unsigned long long* p = reinterpret_cast<unsigned long long*>(&result);  // Break down into two 64 bit words
  int c = __builtin_popcountll(p[0]) + __builtin_popcountll(p[1]); // number of one bits in w
  return c;
}

// don't use, this is done better above
int
hamming_distance_with_bits_old(myint128 x, myint128 y)
{
  myint128 w = x ^ y;
  unsigned long long* p = reinterpret_cast<unsigned long long*>(&w);  // Break down into two 64 bit words
  int c = __builtin_popcountll(p[0]) + __builtin_popcountll(p[1]); // number of one bits in w
  if (c <= 1)
    return c;
  else if (c > 2)
    return 2;
  else {
    if ((w&(w>>1) & evenmask) != 0)
      return 1;
    else 
      return 2;
  }
}

/*
template <typename T>
T
dna_to_number(const std::string& s)  
{
  
  //  assert(s.length() < std::numeric_limits<T>::digits/2);

  T sum=0;
  for (int i=0; i < s.length(); ++i)
    sum = sum*4 + to_int(s[i]);

  return sum;
}
*/



// Take only unique set of sequences
// This implementation is slower, but works for reads longer than 64 bp
std::vector<std::string>
remove_duplicate_reads_generic(const std::vector<std::string>& sequences_orig, int hamming_distance)
{
  std::vector<std::string> sequences;
  if (sequences_orig.size() == 0)
    return std::vector<std::string >();

  std::map<std::string, int> read_count;
  BOOST_FOREACH(std::string s, sequences_orig) {
    ++read_count[s];
  }

  if (hamming_distance == 0) {
    std::string str;
    int count;

    BOOST_FOREACH(boost::tie(str, count), read_count)
      sequences.push_back(str);
    return sequences;
  }

  std::multimap<int, std::string, std::greater<int> > reads_multimap;
  std::string str;
  int count;
  BOOST_FOREACH(boost::tie(str, count), read_count) {
    reads_multimap.insert(std::make_pair(count, str));
  }
  std::vector<std::pair<int, std::string> > reads_vector;
  BOOST_FOREACH(boost::tie(count, str), reads_multimap) {
    reads_vector.push_back(std::make_pair(count, str));
  }

  int size = reads_vector.size();
  std::vector<bool> deleted(size, false);
  for (int i=0; i < size; ++i) {
    if (deleted[i])
      continue;
    std::string current_str = reads_vector[i].second;
    for (int j=i+1; j < size; ++j) {
      if (deleted[j])
	continue;
      if (::hamming_distance(current_str, reads_vector[j].second) <= hamming_distance) {
	deleted[j] = true;
      }
    }
    sequences.push_back(current_str);
  }
  return sequences;
}



// Take only unique set of sequences
std::vector<std::string>
remove_duplicate_reads_faster(const std::vector<std::string>& sequences_orig, int hamming_distance)
{
  
  std::vector<std::string> sequences;
  if (sequences_orig.size() == 0)
    return std::vector<std::string >();
  const int L = sequences_orig[0].length();
  if (L > 64)
    return remove_duplicate_reads_generic(sequences_orig, hamming_distance);
  
  std::map<myint128, int> read_count;
  BOOST_FOREACH(std::string s, sequences_orig) {
    std::string key = s;
    ++read_count[dna_to_number<myint128>(key)];
  }

  if (hamming_distance == 0) {
    myint128 id;
    int count;

    BOOST_FOREACH(boost::tie(id, count), read_count)
      sequences.push_back(number_to_dna<myint128>(id, L));
    return sequences;
  }

  std::multimap<int, myint128, std::greater<int> > reads_multimap;
  myint128 id;
  int count;
  BOOST_FOREACH(boost::tie(id, count), read_count) {
    reads_multimap.insert(std::make_pair(count, id));
  }
  std::vector<std::pair<int, myint128> > reads_vector;
  BOOST_FOREACH(boost::tie(count, id), reads_multimap) {
    reads_vector.push_back(std::make_pair(count, id));
  }

  int size = reads_vector.size();
  std::vector<bool> deleted(size, false);
  for (int i=0; i < size; ++i) {
    if (deleted[i])
      continue;
    myint128 current_id = reads_vector[i].second;
    for (int j=i+1; j < size; ++j) {
      if (deleted[j])
	continue;
      if (hamming_distance_with_bits(current_id, reads_vector[j].second) <= hamming_distance) {
	deleted[j] = true;
      }
    }
    sequences.push_back(number_to_dna<myint128>(current_id, L));
  }
  /*
  std::vector<std::string> sequences;
  std::vector<myint128> duplicates;
  std::vector<myint128> unique;
  myint128 b;
  int count;
  BOOST_FOREACH(boost::tie(b, count), read_count) {
    if (count > 1)
      duplicates.push_back(b);
    else
      unique.push_back(b);
    sequences.push_back(number_to_dna<myint128>(b, L));
  } 
  printf("Dataset contains %zu distinct sequences with multiple occurrences\n", duplicates.size());
  printf("Dataset contains %zu unique sequences\n", read_count.size() - duplicates.size());
  std::vector<std::string> sequences2;
  BOOST_FOREACH(std::string s, sequences_orig) {

    myint128 si  = dna_to_number<myint128>(s);
    BOOST_FOREACH(myint128 d, duplicates) {
      if (hamming_distance_with_bits(si, d) <= 1)
	goto end;
    }
    
    sequences2.push_back(s);
  end:
    ;
  }
  printf("If Hamming duplicates were removed, %zu sequences would remain", sequences2.size());
  */
  /*
  BOOST_FOREACH(myint128 si, unique) {

    std::string s = number_to_dna(si, L);
    std::string sr = reverse_complement(s);
    myint128 sri = dna_to_number<myint128>(sr);
    BOOST_FOREACH(myint128 d, duplicates) {
      if (hamming_distance_with_bits(si, d) <= 1 or hamming_distance_with_bits(sri, d) <= 1)
	goto end;
    }
    
    sequences.push_back(s);
  end:
    ;
  }
  */

  return sequences;
}


// returns number of scoring matrix hits
int
number_of_matches(const std::vector<std::string>& sequences, 
		  const matrix<double>& m, 
		  double threshold)
{
  int lines = sequences.size();
  assert(lines > 0);
  int L = sequences[0].length();
  int k = m.get_columns();
  score_function_t score_function;

  //  if (use_logodds_score)
  if (true)
    score_function = compute_logodds_probability;
  else
    score_function = compute_normal_probability;

  dmatrix matrix = m;
  dmatrix rev_comp_matrix = reverse_complement(matrix);
  int matching_sequences = 0;
  int matching_sites = 0;
  for (int i=0; i<lines; ++i) {
    const std::string& line = sequences[i];
    int limit = L-k+1;
    bool line_hit = false;
    for (int j=0; j < limit; ++j) {  // go through all possible starting positions
      double p1, p2;

      const std::string& site = line.substr(j, k);
      bool is_palindrome = is_palindromic(site);
      p1 = score_function(site, matrix);
      p2 = score_function(site, rev_comp_matrix);

      if (p1 >= threshold) {  // compare to cutoff
	//matches[i].push_back(occurrence(j,1,i));
	++matching_sites;
	line_hit = true;
	if (is_palindrome)  // count palindromes only once
	  continue;
      }
      
      if (p2 >= threshold) {  // compare to cutoff
	//matches[i].push_back(occurrence(j,-1,i));
	++matching_sites;
	line_hit = true;
      }
    }
    if (line_hit)
      ++matching_sequences;
  }

  return matching_sites;
}
