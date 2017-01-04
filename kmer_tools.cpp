/*

    MODER is a program to learn DNA binding motifs from SELEX datasets.
    Copyright (C) 2016, 2017  Jarkko Toivonen

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
#include "kmer_tools.hpp"
#include "common.hpp"
#include "data.hpp"
#include "parameters.hpp"

#include <cmath>
#include <cassert>

bool use_middle_gap = false;

namespace detail {
  reverse_4_2bit_string_type my_reverse_4_2bit_string;
}

//namespace boost {

// Hash function to make int128 work with boost::hash.
std::size_t 
hash_value(const myuint128& input)
{
    boost::hash<unsigned long long> hasher;
    unsigned long long hashVal = 0;
    unsigned long long mask = (1ull<<32) - 1;
    
    for(int i = 0; i < 3; ++i)
    {
    	hashVal *= 37;
    	hashVal += (mask & (input>>(32*i)));
    }

    return hasher(hashVal);
}

//}




/*

uint32_t
reverse32_2bitstring(uint32_t u, int k)
{
  unsigned char* p=reinterpret_cast<unsigned char*>(&u);
  for (int i=0; i<4; ++i)
    p[i]=my_reverse_4_2bit_string(p[i]);
  std::swap(p[0], p[3]);
  std::swap(p[1], p[2]);

  return u>>(32-2*k);
}

uint64_t
reverse64_2bitstring(uint64_t u, int k)
{
  unsigned char* p=reinterpret_cast<unsigned char*>(&u);
  for (int i=0; i<8; ++i)
    p[i]=my_reverse_4_2bit_string(p[i]);
  std::swap(p[0], p[7]);
  std::swap(p[1], p[6]);
  std::swap(p[2], p[5]);
  std::swap(p[3], p[4]);

  return u>>(64-2*k);
}
*/

/* Not needed anymore

big_int
reverse_complement_2bitstring_big_int(big_int c, int k)
{
  const int number_of_bits = sizeof(big_int)*8;
  assert(k <= number_of_bits/2);
  assert(k > 0);

  big_int result = reverse_2bitstring(c, k);

  big_int one = 1;
  
  big_int mask = (one<<(2*k))-1; // Don't replace one by 1 !!!!!!!!!!!!!!!!
  result = result^mask;  // do the complementation
  //assert(number_to_dna(result, l) == reverse_complement(number_to_dna(c, l)));
  return result;
}

*/

// get counts for each k-mer in data
void
get_kmer_counts(const std::vector<std::string>& sequences, int k,
		count_container_t& count, bool use_two_strands, bool count_palindromes_twice)
{

  assert(sequences.size() > 0);
  assert(k > 0 and k <= sizeof(big_int)*4);

#ifndef USE_HASH
  code_t size = pow(4,k);
  count.assign(size, 0);
#endif

  big_int limit = pow(4, k);
  //printf("k=%i limit=%i\n", k, limit);
  for (int i = 0; i < sequences.size(); ++i) {
    const std::string& line = sequences[i];
    const int L = line.length();
    const int positions = L-k+1;
    std::string s = line.substr(0, k);
			    
    big_int temp=dna_to_number(s);
    bool palindrome = is_palindromic(s);
    assert(temp < limit);
    ++count[temp];
    if (use_two_strands and (not palindrome or count_palindromes_twice))
      ++count[reverse_complement_2bitstring(temp, k)];
    for (int j=1; j < positions; ++j) {
      temp = (temp<<2) + to_int(line[k+j-1]);
      temp &= (limit-1);
      assert(temp < limit);
      ++count[temp];
      palindrome = (temp == reverse_complement_2bitstring(temp, k));
      if (use_two_strands and (not palindrome or count_palindromes_twice))
	++count[reverse_complement_2bitstring(temp, k)];
    }
  }
  //printf("Sizeof(count)=%lu\n", count.size());
}


kmers::kmers(const std::vector<std::string>& sequences, int maxk_) :
  maxk(maxk_), counts(maxk_+1), total_counts(maxk_+1) 
{
  for (int k=1; k <= maxk; ++k) {
    get_kmer_counts(sequences, k, counts[k], use_two_strands);
    total_counts[k] = sum(counts[k]);
    printf("Total_counts[%i]=%i\n", k, total_counts[k]);
  }
}

unsigned int 
kmers::count(int k, big_int code) const 
{
  assert(k >= 1 && k <= maxk);
  return counts[k][code];
}

unsigned int
kmers::count(int k, const std::string& str) const
{
  return count(k, dna_to_number(str));
}

double 
kmers::probability(int k, big_int code) const
{
  assert(k >= 1 && k <= maxk);
  return counts[k][code] / (double) total_counts[k];
}

double
kmers::probability(int k, const std::string& str) const
{
  return probability(k, dna_to_number(str));
}

namespace {

big_int
first_part(big_int c, int k, int gap, int pos)
{
  const int len1=pos;
  const int len2=k-pos;
  const big_int mask1 = (1<<(2*len1))-1;
  big_int first = (c >> (2*(gap+len2)) & mask1);
  return first;
}

big_int
second_part(big_int c, int k, int gap, int pos)
{
  const int len2=k-pos;
  const big_int mask2 = (1<<(2*len2))-1;
  big_int second = c & mask2;
  return second;
}

}

big_int
remove_gap_from_bitstring(big_int c, int k, int gap, int pos)
{
  const big_int mask = (1<<(2*k))-1;
  if (gap==0)
    return c & mask;
  big_int first = first_part(c, k, gap, pos);
  big_int second = second_part(c, k, gap, pos);

  big_int limit = 1 << (2*k);
  const int len2=k-pos;
  big_int result = (first << (2*len2)) | second;
  //printf("orig: %lli first: %lli second: %lli result: %lli\n", c, first, second, result);
  assert(result < limit);
  return result;
}

namespace {

template<typename T>
T
first_part_template(T c, int k, int gap, int pos)
{
  const int len1=pos;
  const int len2=k-pos;
  T one = 1;
  const big_int mask1 = (one<<(2*len1))-1;
  T first = (c >> (2*(gap+len2)) & mask1);
  return first;
}

template<typename T>
T
second_part_template(T c, int k, int gap, int pos)
{
  const int len2=k-pos;
  T one = 1;
  const big_int mask2 = (one<<(2*len2))-1;
  T second = c & mask2;
  return second;
}

}

template<typename T>
T
remove_gap_from_bitstring_template(T c, int k, int gap, int pos)
{
  T one = 1;
  const big_int mask = (one<<(2*k))-1;
  if (gap==0)
    return c & mask;
  big_int first = first_part_template(c, k, gap, pos);
  big_int second = second_part_template(c, k, gap, pos);

  big_int limit = one << (2*k);
  const int len2=k-pos;
  big_int result = (first << (2*len2)) | second;
  assert(result < limit);
  return result;
}


// this is almost equivalent to helper3, except that 'pos' argument works here
void
get_gapped_kmer_counts(const std::vector<std::string>& sequences, int k, int max_gap,
		       gapped_count_container_t& count)
{
  assert(sequences.size() > 0);
  assert(k > 0);
  assert(max_gap >= 0);
  const int L = sequences[0].length();
  assert(k + max_gap <= L);
  count.resize(boost::extents[max_gap+1][k+1]);
  for (int gap_len=0; gap_len <= max_gap; ++gap_len) {
    const int l = k+gap_len;
    const int positions = L-l+1;

    for (int i = 0; i < sequences.size(); ++i) {
      const std::string& line = sequences[i];
      const std::string& line2 = reverse_complement(line);
      
      for (int gap_pos=1; gap_pos < (gap_len == 0 ? 2 : k); ++gap_pos) {
	big_int id1=dna_to_number(line.substr(0, l));
	big_int id2=dna_to_number(line2.substr(0, l));
	++count[gap_len][gap_pos][remove_gap_from_bitstring(id1, k, gap_len, gap_pos)];
	++count[gap_len][gap_pos][remove_gap_from_bitstring(id2, k, gap_len, gap_pos)];
	for (int j=1; j < positions; ++j) {
	  id1 = (id1<<2) + to_int(line[l+j-1]);
	  id2 = (id2<<2) + to_int(line2[l+j-1]);
	  ++count[gap_len][gap_pos][remove_gap_from_bitstring(id1, k, gap_len, gap_pos)];
	  ++count[gap_len][gap_pos][remove_gap_from_bitstring(id2, k, gap_len, gap_pos)];
	}
      }  // gap_pos
    }
  }
}

template <typename T>
inline
T
bit_substr(T code, int L, int pos, int len)
{
  //  assert(L >= pos + len);
  T one = 1;
  T mask = (one << (2*len)) - 1;
  return (code >> (2*(L - (pos + len)))) & mask;
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

void
get_gapped_kmer_counts_fast(const std::vector<std::string>& sequences, int k, int max_gap,
			    gapped_count_container_t& count)
{
  int lines = sequences.size();
  assert(lines > 0);
  assert(k > 0);
  assert(max_gap >= 0);
  const int L = sequences[0].length();
  assert(k + max_gap <= L);
  count.resize(boost::extents[max_gap+1][k]);
#ifndef USE_HASH
  for (int gap=0; gap <= max_gap; ++gap)    // Initialize containers if not using hashing
    for (int pos=1; pos < k; ++pos)
      count[gap][pos].assign(pow(4, k), 0);
#endif
  typedef myuint128 int_type;
  std::vector<int_type> bits(lines);
  std::vector<int_type> bits2(lines);

  #pragma omp parallel for schedule(static)
  for (int i=0; i<lines; ++i) {
    bits[i] = dna_to_number<int_type>(sequences[i]);
    bits2[i] = dna_to_number<int_type>(reverse_complement(sequences[i]));
  }
    
  #pragma omp parallel for schedule(static)
  for (int gap_len=0; gap_len <= max_gap; ++gap_len) {
    const int l = k+gap_len;
    const int positions = L-l+1;

    int first, last;
    if (gap_len == 0)
      first = last = 1;
    else {
      if (use_middle_gap) {
	if (k % 2 == 0)
	  first = last = k / 2;
	else {
	  first = k / 2;
	  last = first + 1;
	}
      } else {
	first = 1;
	last = k-1;
      }
    }


    for (int i = 0; i < lines; ++i) {

      for (int gap_pos=first; gap_pos <= last; ++gap_pos) {
	for (int j=0; j < positions; ++j) {
	  int_type id1 = bit_substr(bits[i], L, j, l);
	  int_type id2 = bit_substr(bits2[i], L, j, l);
	  ++count[gap_len][gap_pos][remove_gap_from_bitstring_template(id1, k, gap_len, gap_pos)];
	  ++count[gap_len][gap_pos][remove_gap_from_bitstring_template(id2, k, gap_len, gap_pos)];
	}
      }  // gap_pos
    }
  }
}
