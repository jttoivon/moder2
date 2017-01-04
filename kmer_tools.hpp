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
#ifndef KMER_TOOLS_HPP
#define KMER_TOOLS_HPP


#include "type.hpp"

#include <vector>
#include <string>

#include <boost/unordered_map.hpp>
#include <boost/multi_array.hpp>
#include "unordered_map.hpp"


extern bool use_middle_gap;




/*
namespace std {
  std::size_t 
  hash_value(const myuint128& input);
}
*/

// Hash function to make int128 work with boost::hash.
struct hash128
    : std::unary_function<myuint128, std::size_t>
{
  std::size_t 
  operator()(const myuint128& input) const
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
};

#define USE_HASH 1


//typedef unsigned long long int code_t;   // These code bit representations of DNA sequences
//typedef myuint128 code_t;
typedef uint64_t code_t;

typedef unsigned int count_t;

///////////////////////////
//
// Define count_container_t
//
#ifdef USE_HASH
//typedef boost::unordered_map<code_t, count_t> count_container_t;
typedef my_unordered_map<code_t, count_t, hash128> count_container_t;
#else
typedef std::vector<count_t> count_container_t;   // index type of std::vector is size_t
#endif


typedef boost::multi_array<count_container_t, 2> gapped_count_container_t; 


void
get_kmer_counts(const std::vector<std::string>& sequences, int k,
		count_container_t& count, bool use_two_strands = true, bool count_palindromes_twice = false);

void
get_gapped_kmer_counts(const std::vector<std::string>& sequences, int k, int gaps,
		       gapped_count_container_t& count);

void
get_gapped_kmer_counts_fast(const std::vector<std::string>& sequences, int k, int max_gap,
			    gapped_count_container_t& count);


namespace detail {

  class reverse_4_2bit_string_type  // this functor reverses the 4-mer stored in a byte
  {
  public:

    reverse_4_2bit_string_type() : v(256) {        // initialise the v vector
      unsigned char mask[4] = {0xc0, 0x30, 0x0c, 0x03};
      for (int i=0; i < 256; ++i) {
	unsigned char x0=(i & mask[0])>>6;
	unsigned char x1=(i & mask[1])>>2;
	unsigned char x2=(i & mask[2])<<2;
	unsigned char x3=(i & mask[3])<<6;
	v[i] = x0 | x1 | x2 | x3;
      }
    }

    unsigned char
    operator()(unsigned char i) { return v[i]; }

  private:
    std::vector<unsigned char> v;
  };

  extern reverse_4_2bit_string_type my_reverse_4_2bit_string;

} // end namespace detail



/*
uint32_t
reverse32_2bitstring(uint32_t u, int k);

uint64_t
reverse64_2bitstring(uint64_t u, int k);
*/

template <typename T>
T
reverse_2bitstring(T u, int k)
{
  static_assert(std::is_unsigned<T>::value,
		"reverse_2bitstring requires unsigned parameter type");
  unsigned char* p=reinterpret_cast<unsigned char*>(&u);
  int bytes = sizeof(T);
  int bits = bytes * 8;
  
  for (int i=0; i < bytes; ++i)
    p[i]=detail::my_reverse_4_2bit_string(p[i]);
  for (int i=0; i < bytes/2; ++i)
    std::swap(p[i], p[bytes-i-1]);
    
  return u>>(bits-2*k);
}

template <typename T>
T
reverse_complement_2bitstring(T u, int k)
{
  static_assert(std::is_unsigned<T>::value,
		"reverse_complement_2bitstring requires unsigned parameter type");
  const int number_of_bits = sizeof(T)*8;
  assert(k <= number_of_bits/2);
  assert(k > 0);

  T result = reverse_2bitstring(u, k);

  const T one = 1;
  
  T mask = (one<<(2*k))-1;  // Don't replace one by 1 !!!!!!!!!!!!!!!!
  result = result^mask;     // do the complementation
  //assert(number_to_dna(result, l) == reverse_complement(number_to_dna(c, l)));
  return result;
}

//big_int
//reverse_complement_2bitstring_big_int(big_int c, int l);



big_int
remove_gap_from_bitstring(big_int c, int k, int gap, int pos);


namespace {

template <typename T>
int
mypopcount(T u)
{
  static_assert(std::is_unsigned<T>::value, "mypopcount requires unsigned parameter type");
  unsigned char* p = reinterpret_cast<unsigned char*>(&u);  // Break down into bytes
  int result = 0;
  for (int i=0; i < sizeof(u); ++i)
    result += __builtin_popcount(p[i]);
  return result;
}

template <>
int
mypopcount<uint32_t>(uint32_t u)
{
  return __builtin_popcount(u);
}

template <>
int
mypopcount<uint64_t>(uint64_t u)
{
  return __builtin_popcountll(u);
}
  
template <>
int
mypopcount<myuint128>(myuint128 u)
{
  unsigned long long* p = reinterpret_cast<unsigned long long*>(&u);  // Break down into two 64 bit words
  int c = __builtin_popcountll(p[0]) + __builtin_popcountll(p[1]); // number of one bits in w
  return c;
}


} // end nameless namespace


namespace detail {

  template <typename T>
  T
  get_even_mask()
  {
    T x=0;
    int number_of_bits = sizeof(T) * 8;
    for (int i=0; i < number_of_bits / 2; ++i) {
      x <<= 2;
      x += 1;
    }
    return x;
  }

} // end namespace detail


template <typename T>
int
hamming_distance_with_bits(T x, T y)
{
  static_assert(std::is_unsigned<T>::value,
		"hamming_distance_with_bits requires unsigned parameter type");
  // Note! Here we assume that the leftmost bits are cleared and contain no thrash.
  // We could also get length of string as parameter and mask out the extra bits, just to be sure, but this would slow things down.
  static const T evenmask = detail::get_even_mask<T>();
  static const T oddmask = evenmask << 1;
  T w = x ^ y;
  T result = ((w & oddmask) >> 1) | (w & evenmask);
  return mypopcount<T>(result);

}


class kmers
{
public:

  kmers(const std::vector<std::string>& sequences, int maxk_);

  unsigned
  count(int k, big_int code) const ;


  unsigned
  count(int k, const std::string& str) const ;

  double 
  probability(int k, big_int code) const ;


  double
  probability(int k, const std::string& str) const;

  int get_maxk() const 
  { return maxk; }

private:
  int maxk;
  std::vector<count_container_t> counts;
  std::vector<unsigned int> total_counts;
};



#endif // KMER_TOOLS_HPP
