#ifndef DATA_HPP
#define DATA_HPP

#include "matrix.hpp"

#include <string>
#include <vector>
#include <cassert>

typedef unsigned long long int big_int;


class my_to_int_type
{
public:

  my_to_int_type() 
  {
    std::string nucleotides = "ACGTN";
    std::string nucleotides_lowercase = "acgtn";
    
    for (int i=0; i < nucleotides.length(); ++i)
      to_int_array[(unsigned char)nucleotides[i]] = i;
    for (int i=0; i < nucleotides_lowercase.length(); ++i)
      to_int_array[(unsigned char)nucleotides_lowercase[i]] = i;
  }

  int
  operator()(char c)
  {
    return to_int_array[(unsigned char)c];
  }

private:
  int to_int_array[256];
};

static my_to_int_type to_int;




// interpret string of nucleid acids as a number in base 4
// The first nucleotide of the string is in the most significant end
big_int
dna_to_number(const std::string& s);

template <typename T>
T
dna_to_number(const std::string& s)  
{
  assert(s.length() <= (sizeof(T)*4));

  T sum=0;
  for (int i=0; i < s.length(); ++i)
    sum = sum*4 + to_int(s[i]);
  
  return sum;
}

// The first nucleotide of the string is in the most significant end
std::string
number_to_dna(big_int i, int l);  // converts integer back to a nucled acid sequence

template <typename T>
std::string
number_to_dna(T i, int l)  
{
  assert(i>=0);

  char nucs[] = "ACGT";

  std::string s(l, 'X');
  while (l-- > 0) {
    s[l] = nucs[i & 3];
    i = i / 4;
  }
  return s;
}


big_int
extended_dna_to_number(const std::string& s);

// converts integer back to a nucled acid sequence including N
// Sequence has l nucleotides
std::string
number_to_extended_dna(big_int i, int l);










// interpret string of nucleid acids as a number in base 4
// The first nucleotide of the string is in the most significant end
big_int
vector_to_number(const std::vector<int>& s);


// converts integer back to a nucled acid sequence
// Sequence has l nucleotides
std::vector<int>
number_to_vector(big_int i, int l);

// interpret string of nucleid acids A,C,G,T,N as a number in base 5
// The first nucleotide of the string is in the most significant end
big_int
extended_vector_to_number(const std::vector<int>& s);

// converts integer back to a nucled acid sequence including N
// Sequence has l nucleotides
std::vector<int>
number_to_extended_vector(big_int i, int l);


int
N_count(const std::string& s);

int
non_N_count(const std::string& s);

bool
has_atmost_one_gap(const std::string& s);

std::vector<std::string>
remove_duplicate_reads_faster(const std::vector<std::string>& sequences_orig, int hamming_distance);

// returns number of scoring matrix hits
int
number_of_matches(const std::vector<std::string>& sequences, 
		  const dmatrix& m, 
		  double threshold);

#endif
