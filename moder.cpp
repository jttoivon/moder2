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

///////////////////////////////
//
// (c) Jarkko Toivonen
//
// email: jarkko.toivonen@cs.helsinki.fi
//
///////////////////////////////

#define TIMING

//extern "C" {
//#include "suffix_array/interface.h"
//}

#ifndef PACKAGE_VERSION
#define PACKAGE_VERSION "unknown"
#endif

#include "iupac.hpp"
#include "matrix.hpp"
#include "vectors.hpp"
#include "timing.hpp"
#include "matrix_tools.hpp"
#include "common.hpp"
#include "probabilities.hpp"
#include "parameters.hpp"
#include "orientation.hpp"
#include "kmer_tools.hpp"
#include "multinomial_helper.hpp"
#include "suffix_array_wrapper.hpp"

#include <sys/stat.h>
#include <cstdio>
#include <cmath>
#include <ctime>
#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif
#include <fenv.h>

#include <libgen.h>   // for dirname and basename

#ifdef _OPENMP
#include <omp.h>
#endif

#include <string>
#include <iostream>
#include <fstream>
#include <numeric>
#include <queue>
#include <cfloat>

#include <boost/algorithm/string.hpp>
#include <boost/foreach.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/program_options.hpp>
namespace po = boost::program_options;

// clang does not seem to be able to handle 128 bit type (long double)
// in openmp constructs
#if defined(__clang__)
typedef double FloatType;
#else
typedef long double FloatType;
#endif

typedef boost::tuple<int,int> cob_combination_t;

typedef std::vector<boost::tuple<int, int, int, int, std::string> > overlapping_dimer_cases_t;
typedef std::vector<boost::tuple<int, int, int, int> > spaced_dimer_cases_t;

typedef boost::multi_array<dmatrix, 2> cob_of_matrices;




bool use_palindromic_correction=false;
bool use_multimer=true;
bool use_meme_init=false;
bool get_full_flanks=false;  // make motifs for each model that have width 2*L-motif_width

int max_iter = 50;  // was 300
int minimum_distance_for_learning = 4;
int global_dmax = 10;
//int global_max_dist_for_deviation = -1;
//int global_max_dist_for_deviation = 4;
int global_max_dist_for_deviation = 1000;
//int min_flank = 3;  // On each side of overlapping area of a dimer at least min_flank number of positions
                    // must be reserved. That is the minimum non-overlapping part on each side.
double ic_threshold = 0.40;
double learning_fraction = 0.02;
//int hamming_threshold = 4;

int character_count=0;
int digram_count=0;

std::vector<double> background_frequencies(4);
std::vector<double> background_probabilities(4);
dmatrix background_frequency_matrix(4,4);   // for Markov background
dmatrix background_probability_matrix(4,4); // for Markov background


prior<double> pseudo_counts;

double cob_cutoff = 0.001;  // if an element in a cob table is smaller than this constant,
                           // the element is excluded. Note! This works for spaced dimers as well

bool adjust_seeds = true;
bool use_multinomial=true;
bool local_debug = true;
bool extra_debug = false;   // Even more printing
bool allow_extension = false;
bool use_dimers = true;
bool seeds_given = false;
bool no_unique = false;
bool use_output = false; // whether to write model parameters to files
bool maximize_overlapping_seeds=true;
bool require_directional_seed = false;
bool avoid_palindromes = false; // If tf1==tf2, the orientation is HH or TT and the PPM \tau_ht1,ht2,o,d, the
                                // probability of a sequence is the same in both directions.
                                // If we recognize this situation, we can try to escape from this palindromicity
                                // by using temporarily only a single strand.
int hamming_distance_overlapping_seeds_N=0;
int hamming_distance_overlapping_seeds_OR=2;
int hamming_radius = 1; // For learning the model from sequences in the hamming neighbourhood of the seed

std::string unbound; // Filename where to store the sequences that had greatest probability under background model
std::vector<std::string> names;
std::string outputdir = ".";

typedef std::string (*combine_seeds_func_ptr)(const std::string& seed1, const std::string& seed2, int d);

enum { COMBINE_SEEDS_OR=0, COMBINE_SEEDS_AND=1, COMBINE_SEEDS_N=2};
int default_combine_method = COMBINE_SEEDS_N;
combine_seeds_func_ptr combine_seeds_func;

struct combine_seeds_method_type {
  const char* name;
  combine_seeds_func_ptr function;
};

std::string
combine_seeds_OR(const std::string& seed1, const std::string& seed2, int d);

std::string
combine_seeds_AND(const std::string& seed1, const std::string& seed2, int d);

std::string
combine_seeds_masked(const std::string& seed1, const std::string& seed2, int d);

combine_seeds_method_type combine_seeds_methods[] =
  {
    { "UNION", combine_seeds_OR},
    { "INTERSECTION", combine_seeds_AND},
    { "N", combine_seeds_masked}
  };

template <typename T>
class BinOp
{
public:
  typedef const T& (*type)(const T&, const T&);
};


//template <typename T, typename BinaryOperation>
template <typename T>
class myaccumulate
{
public:
  typedef typename BinOp<T>::type BinaryOperation;

  typedef const T& (*func_ptr)(const T&, const T&);


  myaccumulate(T init, BinaryOperation op_) : result(init), op(op_) {}

  void 
  operator()(T v)
  {
    result = op(result, v);
  }

  T get() const { return result; }

  void
  reset(T t = T())
  { result=t; }

private:
  T result;
  BinaryOperation op;
};

std::string
combine_seeds_OR(const std::string& seed1, const std::string& seed2, int d)
{
  assert(d < 0);
  int k1 = seed1.length();
  int k2 = seed2.length();
  std::string result(k1 + k2 + d, '-');

  for (int i=0; i < k1 + d; ++i)
    result[i] = seed1[i];

  for (int i=k1+d; i < k1; ++i)
    result[i] = iupac_class.bits_to_char(iupac_class.char_to_bits(seed1[i]) | iupac_class.char_to_bits(seed2[i-k1-d]));
  for (int i=k1; i < k1+k2+d; ++i)
    result[i] = seed2[i-k1-d];

  return result;
}


std::string
combine_seeds_AND(const std::string& seed1, const std::string& seed2, int d)
{
  assert(d < 0);
  int k1 = seed1.length();
  int k2 = seed2.length();
  std::string result(k1 + k2 + d, '-');

  for (int i=0; i < k1 + d; ++i)
    result[i] = seed1[i];

  for (int i=k1+d; i < k1; ++i)
    result[i] = iupac_class.bits_to_char(iupac_class.char_to_bits(seed1[i]) & iupac_class.char_to_bits(seed2[i-k1-d]));
  for (int i=k1; i < k1+k2+d; ++i)
    result[i] = seed2[i-k1-d];

  return result;
}

std::string
combine_seeds_masked(const std::string& seed1, const std::string& seed2, int d)
{
  assert(d < 0);
  int k1 = seed1.length();
  int k2 = seed2.length();
  std::string result(k1 + k2 + d, 'N');
  int x = 0; // Make the gap wider by x from both sides
  for (int i=0; i < k1 + d - x; ++i)   // left flank
    result[i] = seed1[i];

  for (int i=k1 + x; i < k1+k2+d; ++i) // right flank
    result[i] = seed2[i-k1-d];

  return result;
}



// gives descending order on the first position of the triple
class comp_first_desc
{
public:

  bool operator()(boost::tuple<double, int, int> x, 
		  boost::tuple<double, int, int> y) const {
    return boost::get<0>(x) > boost::get<0>(y);
  }

};


class gapped_kmer_context
{
public:

  std::string
  make_string(const std::vector<std::string>& sequences)
  {
    std::string str1=join(sequences, '#');
    if (use_two_strands) {
      str1.append("#");
      str1 += join_rev(sequences, '#');
    }

    return str1;
  }

  gapped_kmer_context(const std::vector<std::string>& sequences) : sa(make_string(sequences))
  {
    timer = 0.0;
  }


  unsigned long int
  count(const std::string& pattern) const
  {
    return sa.count_iupac(pattern);
  }

  // Finds a string with maximum count from a set of strings.
  // The set of strings is specified by an iupac sequence 's'
  boost::tuple<std::string,int>
  max(const std::string& s) const
  {
    TIME_START(s);
    int k = s.length();
    assert(is_iupac_string(s));

    //int max_count = 0;
    std::string arg_max=s;

    std::vector<int> v(k, 0);  // helper array

    std::vector<int> base(k, 0);
    std::vector<const char*> iupac_classes(k);
    std::string pattern(k, '-');
    int number_of_combinations = 1;
    for (int i=0; i < k; ++i) {
      iupac_classes[i] = iupac_class(s[i]);
      int size = strlen(iupac_classes[i]);
      base[i] = size - 1;
      pattern[i] = iupac_classes[i][0];            // Initialize pattern to the first sequence of iupac string
      number_of_combinations *= size;
    }
    v[k-1]=-1;  // Initialize
    // The loop goes through nucleotide sequences that belong to
    //given iupac sequence in alphabetical order
    std::multimap<unsigned long int, std::string> kmers;
    for (int j=0; j < number_of_combinations; ++j) {
      
      int i;
      for (i=k-1; v[i] == base[i]; --i) {
	v[i]=0;
	pattern[i] = iupac_classes[i][v[i]];
      }
      v[i]++;
      pattern[i] = iupac_classes[i][v[i]];

    /* Writes in numocc the number of occurrences of the substring 
       pattern[0..length-1] found in the text indexed by index. */

      unsigned long int number_of_occurrences;
      number_of_occurrences = count(pattern);
      kmers.insert(std::make_pair(number_of_occurrences, pattern));
      
      // if (number_of_occurrences > max_count) {
      // 	max_count = number_of_occurrences;
      // 	arg_max = pattern;
      // }
    }
    unsigned long int count=0;
    BOOST_REVERSE_FOREACH(boost::tie(count, pattern), kmers) {
      if (require_directional_seed and palindromic_index(pattern) <= 1) //is_palindromic(pattern))
	;
      else {
	break;
      }
	
    }
    timer += TIME_GET(s);

    return boost::make_tuple(pattern, count);
    //    return boost::make_tuple(arg_max, max_count);
  }

  // Allow mismatch ('N') in 'errors' positions in the subsequence s_orig[begin,begin+len)
  boost::tuple<std::string,int>
  max(const std::string& s_orig, int errors, int begin=0, int len=-1) const
  {
    assert(errors == 0 or errors == 1 or errors == 2);
    if (errors == 0)
      return max(s_orig);
    std::string s = s_orig;
    int k = s.length();
    if (len == -1)
      len = k;
    int max_count=0;
    std::string arg_max = s;
    char old_i = 'N';
    char old_j = 'N';
    std::string result;
    int count;
    int end = begin+len;
    if (errors == 1) {
      for (int i=begin; i < end; ++i) {
	if (s[i] == 'N')
	  continue;
	std::swap(s[i], old_i);
	boost::tie(result, count) = max(s);
	if (count > max_count) {
	  max_count = count;
	  arg_max = result;
	}
	std::swap(s[i], old_i);
      }
      return boost::make_tuple(arg_max, max_count);
    }
    else if (errors == 2) {
      for (int i=begin; i < end-1; ++i) {
	if (s[i] == 'N')
	  continue;
	std::swap(s[i], old_i);
	for (int j=i+1; j < end; ++j) {
	  if (s[j] == 'N')
	    continue;
	  std::swap(s[j], old_j);
	  boost::tie(result, count) = max(s);
	  if (count > max_count) {
	    max_count = count;
	    arg_max = result;
	  }
	  std::swap(s[j], old_j);
	}  // end for j
	std::swap(s[i], old_i);
      } // end for i
      return boost::make_tuple(arg_max, max_count);
    }
    return boost::make_tuple(std::string(""), 0);
  }

  /*
  // Old version, do not use
  std::string
  max(const std::string& s, int gap_pos, int gap_len) const
  {
    TIME_START(s);
    int k = s.length();
    assert(gap_pos > 0);
    assert(gap_pos + gap_len < k);
    for(int i=gap_pos; i < gap_pos + gap_len; ++i) {
      assert(s[i] == 'N');
    }

// Writes in numocc the number of occurrences of the substring 
//       pattern[0..length-1] found in the text indexed by index.

    int max_count = 0;
    std::string arg_max;
    unsigned long int number_of_occurrences;
    char nucs[] = "ACGT";

    std::vector<int> v(gap_len, 0);  // helper array
    std::string pattern = s;
    v[gap_len-1]=-1;  // Initialize
    // The loop goes throught the gap filler in the following order:
    // AAAA....AAA
    // AAAA....AAC
    // ...
    // TTTT....TTG
    // TTTT....TTT
    for (int j=0; j < pow(4, gap_len); ++j) {

      int i;
      for (i=gap_len-1; v[i] == 3; --i) {
	v[i]=0;
	pattern[gap_pos+i] = nucs[v[i]];
      }
      v[i]++;
      pattern[gap_pos+i] = nucs[v[i]];

      number_of_occurrences = count(pattern);

      if (number_of_occurrences > max_count) {
	max_count = number_of_occurrences;
	arg_max = pattern;
      }
    }

    timer += TIME_GET(s);

    return arg_max;
  }
*/

  ~gapped_kmer_context() 
  { 
    printf("Gapped kmer context spend %.2f seconds in total finding maxima\n", timer);
  }

private:
  suffix_array sa;
  mutable double timer;
};


// return -1 if Hamming distance is 0,
//        -2 if Hamming distance is greater than one,
//        pos of the mismatch if distance equals one
int
iupac_hamming_dist_helper(const std::string& str, const std::string& pattern)
{
  assert(str.length() == pattern.length());
  int result = -1;
  for (int i=0; i < str.length(); ++i) {
    if (not iupac_match(str[i], pattern[i])) {
      if (result == -1)
	result = i;
      else
	return -2;
    }
  }
  
  return result;
}


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

template <typename T, typename F>
matrix<T>
log2(const matrix<F>& orig)
{
  matrix<T> result(orig.dim());
  int rows=orig.get_rows();
  int cols=orig.get_columns();
  for (int i=0; i < rows; ++i)
    for (int j=0; j < cols; ++j)
      result(i, j) = log2l(orig(i, j));

  return result;
}

template <typename T>
T
log_sum(std::priority_queue<T, std::vector<T>, std::greater<T> >& queue)
{
  do {
    T q = queue.top(); queue.pop();
    T p = queue.top(); queue.pop();
    assert (p >= q);
    queue.push(p + log2l(1+exp2l(q-p)));
  } while (queue.size() > 1);
  
  return queue.top();
}

template <typename T>
T
log_sum(std::vector<T>& v)
{
  T p = v.back(); v.pop_back();
  do {
    T q = v.back(); v.pop_back();
    if (q > q)
      std::swap(p, q);
    p = p + log2l(1+exp2l(q-p));
  } while (v.size() > 0);
  
  return p;
}

typedef boost::multi_array<FloatType, 4> array_4d_type; // for Z variable, indices are (i,k,dir,j)
typedef boost::multi_array<FloatType, 5> array_5d_type; // for Z variable, indices are (i,o,d,dir,j)


struct cob_params_t
{
  typedef boost::multi_array<double, 2>::extent_range range;

  cob_params_t(int tf1_, int tf2_, 
	       const boost::multi_array<double, 2>& dimer_lambdas_,
	       const boost::multi_array<std::string, 2>& dimer_seeds_,
	       const boost::multi_array<dmatrix, 2>& overlapping_dimer_PWM_,
	       const std::vector<std::string>& fixed_seeds_, const std::vector<int>& L_, int dmin_, int dmax_,
	       int max_dist_for_deviation_)
    :  tf1(tf1_), tf2(tf2_), dimer_lambdas(dimer_lambdas_), dimer_seeds(dimer_seeds_),
       overlapping_dimer_PWM(overlapping_dimer_PWM_), L(L_), dmin(dmin_), dmax(dmax_),
       max_dist_for_deviation(max_dist_for_deviation_)
  {
    k1 = fixed_seeds_[tf1].length();
    k2 = fixed_seeds_[tf2].length();
    //    assert(dmin == -std::min(k1, k2) + min_flank);
    
    number_of_orientations = tf1 == tf2 ? 3 : 4;
    expected_overlapping_dimer_PWMs.resize(boost::extents[number_of_orientations][range(dmin, max_dist_for_deviation+1)]);
    lines = L.size();
    dimer_w.resize(boost::extents[range(dmin, dmax+1)]);
    dimer_m.resize(boost::extents[lines][range(dmin, dmax+1)]);

    deviation.resize(boost::extents[number_of_orientations][range(dmin, max_dist_for_deviation+1)]);
    for (int o=0; o < number_of_orientations; ++o) {
      for (int d=dmin; d <= max_dist_for_deviation; ++d) {
	//	deviation[o][d] = dmatrix(4, 2 - d);
	deviation[o][d] = dmatrix(4, k1 + k2 + d);  // This contains redundant flanks. Just to ease
	                                            // operating with matrices expected and observed
	                                            // which have the same dimensions.
      }
    }
  }

  std::string
  name() const 
  {
    return to_string("%i-%i", tf1, tf2);
  }

  void
  set_lengths() {
    //printf("Overlapping dimer cases:\n");
    //typedef BinOp<int>::type func_ptr;
    myaccumulate<int> acc(0, static_cast<BinOp<int>::type>(std::max));
    for (int d=dmin; d < std::min(0, dmax+1); ++d) {
      dimer_w[d] = overlapping_dimer_PWM[0][d].get_columns();
      for (int i=0; i < lines; ++i) {
	dimer_m[i][d] = L[i] - dimer_w[d] + 1;
	acc(dimer_m[i][d]);
      }
    }  // end for d

    overlapping_dimer_m_max = acc.get();
    
    //printf("Spaced dimer cases:\n");
    //typedef const int& (*func_ptr)(const int&, const int&);
    acc.reset(0);
    for (int d=0; d <= dmax; ++d){
      dimer_w[d] = k1 + d + k2;
      for (int i=0; i < lines; ++i) {
	dimer_m[i][d] = L[i] - dimer_w[d] + 1;  // One past maximum start pos for dimer
	acc(dimer_m[i][d]);
      }
    }  // end for d
      
    spaced_dimer_m_max = acc.get();
  }

  void
  initialize_Z(int lines) {
    int directions = use_two_strands ? 2 : 1;
    overlapping_dimer_Z.resize(boost::extents[lines][number_of_orientations][range(dmin,max_dist_for_deviation+1)][directions][overlapping_dimer_m_max]);
    spaced_dimer_Z.resize(boost::extents[lines][number_of_orientations][range(max_dist_for_deviation+1,dmax+1)][directions][spaced_dimer_m_max]);
    //    overlapping_dimer_Z.resize(boost::extents[lines][number_of_orientations][range(dmin,0)][directions][overlapping_dimer_m_max]);
    //    spaced_dimer_Z.resize(boost::extents[lines][number_of_orientations][range(0,dmax+1)][directions][spaced_dimer_m_max]);
  }

  void
  update_oriented_matrices(const std::vector<dmatrix>& fixed_PWM, const std::vector<std::string>& fixed_seed) {
    oriented_dimer_matrices.resize(4);
    oriented_dimer_seeds.resize(4);
    for (int o=0; o < 4; ++o) {
      oriented_dimer_matrices[o] =
	get_matrices_according_to_hetero_orientation(o, fixed_PWM[tf1], fixed_PWM[tf2]);
      oriented_dimer_seeds[o] =
	get_seeds_according_to_hetero_orientation(o, fixed_seed[tf1], fixed_seed[tf2]);
    }
  }

  void
  compute_expected_matrices(const std::vector<dmatrix>& fixed_PWM) {
    // Compute expected matrices
    for (int o=0; o < number_of_orientations; ++o) {
      for (int d=dmin; d <= max_dist_for_deviation; ++d) {
	//	dmatrix a, b;
	//	boost::tie(a, b) = get_matrices_according_to_hetero_orientation(o, fixed_PWM[tf1],fixed_PWM[tf2]);
	//	boost::tie(a, b) = oriented_dimer_matrices[o];
	dmatrix expected = matrix_product(oriented_dimer_matrices[o].get<0>(), oriented_dimer_matrices[o].get<1>(), d);
	// write_matrix(stdout, expected, to_string("Expected (unnormalized) dimer case matrix %s %s %i:\n", 
	// 					 name().c_str(), orients[o], d).c_str(), "%.6f");
	normalize_matrix_columns(expected);
	//    write_matrix(stdout, expected, to_string("Expected (normalized) dimer case matrix %s %i:\n", orients[o], d).c_str(), "%.6f");
	expected_overlapping_dimer_PWMs[o][d] = expected;
      } // end for d
    }   // end for o
  }

  void
  compute_deviation_matrices()
  {
    const int& w1 = k1;
    for (int o=0; o < number_of_orientations; ++o) {
      for (int d=dmin; d <= max_dist_for_deviation; ++d) {
	int first, last;
	// The interval [first,last] is either the overlap area or the gap.
	if (d < 0) {
	  first = w1 + d;
	  last = w1 - 1;
	} else {
	  first = w1;
	  last = w1 + d - 1;
	}
	for (int row = 0; row < 4; ++row) {
	  for (int column = first - 1; column <= last + 1; ++column) // The overlapping part with flanks of 1bp on both sides
	    deviation[o][d](row,column) = 
	      overlapping_dimer_PWM[o][d](row,column) - expected_overlapping_dimer_PWMs[o][d](row,column);
	}
      }  // end for d
	 
    }  // end for o
 
  }

  int tf1;
  int tf2;
  boost::multi_array<double, 2>      dimer_lambdas;
  boost::multi_array<std::string, 2> dimer_seeds;
  boost::multi_array<dmatrix, 2>     overlapping_dimer_PWM;
  boost::multi_array<dmatrix, 2>     deviation;

  //  const std::vector<std::string>&    fixed_seeds;

  int k1;
  int k2;
  const std::vector<int>& L;
  int lines;
  int dmin;
  int dmax;
  int max_dist_for_deviation;
  int number_of_orientations;

  int overlapping_dimer_m_max;                     // maximum width of a pwm
  int spaced_dimer_m_max;                     // maximum width of a pwm

  boost::multi_array<dmatrix, 2> expected_overlapping_dimer_PWMs;
  boost::multi_array<double, 1> dimer_w;
  boost::multi_array<double, 2> dimer_m;


  array_5d_type overlapping_dimer_Z;
  array_5d_type spaced_dimer_Z;

  std::vector<boost::tuple<dmatrix,dmatrix> > oriented_dimer_matrices;
  std::vector<boost::tuple<std::string,std::string> > oriented_dimer_seeds;
};

/*

// This computes the expected log likelihood of the INCOMPLETE data. that is: E log P(X|total model)
double
incomplete_data_maximum_log_likelihood(const std::vector<dmatrix>& PWM, 
		       const std::vector<double>& q, 
		       const dmatrix& q2, 
		       const std::vector<double>& lambda, 
		       double lambda_bg, std::vector<cob_params_t>& my_cob_params,
		       const std::vector<std::string>& sequences,
		       const std::vector<std::string>& sequences_rev)
{

  double result=0.0;
  int p=PWM.size();
  int L=sequences[0].length();
  int number_of_cobs = my_cob_params.size();

  std::vector<int> w(p);
  std::vector<int> m(p);
  for (int i=0; i < p; ++i) {
    w[i] = PWM[i].get_columns();
    m[i] = L-w[i]+1;
  }

  typedef boost::multi_array<dmatrix, 2> cob_of_matrices_t;
  typedef cob_of_matrices_t::extent_range range;
  std::vector<cob_of_matrices_t> overlapping_models;
  for (int r=0; r < number_of_cobs; ++r) {
    int tf1 = my_cob_params[r].tf1;
    int tf2 = my_cob_params[r].tf2;
    int number_of_orientations = tf1 == tf2 ? 3 : 4;
    int dmin = my_cob_params[r].dmin;
    int max_dist_for_deviation = my_cob_params[r].max_dist_for_deviation;
    //    int dmax = my_cob_params[r].dmax;
    
    overlapping_models.push_back(cob_of_matrices_t(boost::extents[number_of_orientations][range(dmin, max_dist_for_deviation+1)]));
    for (int o=0; o < number_of_orientations; ++o) {
      for (int d=dmin; d <= max_dist_for_deviation; ++d) {
	overlapping_models[r][o][d] = my_cob_params[r].expected_overlapping_dimer_PWMs[o][d] + my_cob_params[r].deviation[o][d];
      }
    }
  }
#pragma omp parallel for shared(use_two_strands) reduction(+:result) schedule(static)
  for (int i=0; i < sequences.size(); ++i) {
    const std::string& line = sequences[i];
    const std::string& line_rev = sequences_rev[i];
    double temp = 0.0;
    for (int k=0; k < p; ++k) {
      for (int j=0; j < m[k]; ++j) {
	temp += compute_probability(line, line_rev, j, 1, PWM[k], q, q2) * lambda[k] / (double)m[k] / 2.0;
	temp += compute_probability(line, line_rev, j, -1, PWM[k], q, q2) * lambda[k] / (double)m[k] / 2.0;
      }
    }

    for (int r=0; r < number_of_cobs; ++r) {
      int tf1 = my_cob_params[r].tf1;
      int tf2 = my_cob_params[r].tf2;
      int number_of_orientations = tf1 == tf2 ? 3 : 4;
      int dmin = my_cob_params[r].dmin;
      int dmax = my_cob_params[r].dmax;
      int max_dist_for_deviation = my_cob_params[r].max_dist_for_deviation;
      for (int o=0; o < number_of_orientations; ++o) {
	for (int d=max_dist_for_deviation+1; d <= dmax; ++d) {       // spaced dimers
	  int m = my_cob_params[r].dimer_m[d];
	  for (int j=0; j < m; ++j) {
	    temp += compute_dimer_probability(line, line_rev,
					      j, +1,
					      my_cob_params[r].oriented_dimer_matrices[o].get<0>(), 
					      my_cob_params[r].oriented_dimer_matrices[o].get<1>(), 
					      d, q, q2) * my_cob_params[r].dimer_lambdas[o][d] / (double)m / 2.0;
	    temp += compute_dimer_probability(line, line_rev,
					      j, -1,
					      my_cob_params[r].oriented_dimer_matrices[o].get<0>(), 
					      my_cob_params[r].oriented_dimer_matrices[o].get<1>(), 
					      d, q, q2) * my_cob_params[r].dimer_lambdas[o][d] / (double)m / 2.0;
	  } // end for j
	} // end for d

	for (int d=dmin; d <= max_dist_for_deviation; ++d) {        // overlapping dimers
	  int m = my_cob_params[r].dimer_m[d];
	  const dmatrix& model = overlapping_models[r][o][d];
	  for (int j=0; j < m; ++j) {
	    temp += compute_probability(line, line_rev, j, +1, model, q, q2) * my_cob_params[r].dimer_lambdas[o][d] / (double)m / 2.0;
	    temp += compute_probability(line, line_rev, j, -1, model, q, q2) * my_cob_params[r].dimer_lambdas[o][d] / (double)m / 2.0;
	  }
	}

      }  // end for o
    } // end for r

    temp += compute_background_probability(line, q, q2) * lambda_bg;
    result += log2(temp);
  } // end for i

  return result;
} // end incomplete_data_log_likelihood

*/

void
print_math_error(int retval)
{
  printf("Error! ");
  if (retval & FE_INVALID)
    printf("FE_INVALID ");
  if (retval & FE_DIVBYZERO)
    printf("FE_DIVBYZERO ");
  if (retval &FE_OVERFLOW )
    printf("FE_OVERFLOW ");
  if (retval & FE_UNDERFLOW)
    printf("FE_UNDERFLOW ");

}


// This computes the expected log likelihood of the COMPLETE data. that is: E log P(X,Z|total model)
double
complete_data_log_likelihood(const std::vector<dmatrix>& PWM, 
			     const std::vector<double>& q, 
			     const dmatrix& q2, 
			     const std::vector<double>& q_rev, 
			     const dmatrix& q2_rev, 
			     const std::vector<double>& lambda, 
			     double lambda_bg, std::vector<cob_params_t>& my_cob_params,
			     const array_4d_type& fixed_Z,
			     const std::vector<std::string>& sequences,
			     const std::vector<std::string>& sequences_rev,
			     const std::vector<int>& fixed_w,
			     const boost::multi_array<int, 2>& fixed_m)
{

  double result=0.0;
  int p=PWM.size();
  //  int L=sequences[0].length();
  int number_of_cobs = my_cob_params.size();

  FloatType l2 = log2l(2);
  FloatType log_lambda_bg = log2l(lambda_bg);
  std::vector<FloatType> log_q(4);
  std::vector<FloatType> log_q_rev(4);
  std::vector<double> log_lambda(p);
  std::vector<matrix<FloatType> > log_PWM;
  for (int i=0; i < 4; ++i) {
    log_q[i] = log2l(q[i]);
    log_q_rev[i] = log2l(q_rev[i]);
  }
  for (int i=0; i < p; ++i) {
    log_lambda[i] = log2l(lambda[i]);
    log_PWM.push_back(log2<FloatType>(PWM[i]));
  }
  boost::multi_array<matrix<FloatType>, 3> log_oriented_dimer_matrices(boost::extents[number_of_cobs][4][2]);
  typedef boost::multi_array<matrix<FloatType>, 2> cob_of_matrices_t;
  typedef cob_of_matrices_t::extent_range range;
  std::vector<cob_of_matrices_t> overlapping_models;
  for (int r=0; r < number_of_cobs; ++r) {
    int tf1 = my_cob_params[r].tf1;
    int tf2 = my_cob_params[r].tf2;
    int number_of_orientations = tf1 == tf2 ? 3 : 4;
    int dmin = my_cob_params[r].dmin;
    int max_dist_for_deviation = my_cob_params[r].max_dist_for_deviation;
    //    int dmax = my_cob_params[r].dmax;
    
    overlapping_models.push_back(cob_of_matrices_t(boost::extents[number_of_orientations][range(dmin, max_dist_for_deviation+1)]));
    for (int o=0; o < number_of_orientations; ++o) {
      log_oriented_dimer_matrices[r][o][0] = log2<FloatType>(my_cob_params[r].oriented_dimer_matrices[o].get<0>());
      log_oriented_dimer_matrices[r][o][1] = log2<FloatType>(my_cob_params[r].oriented_dimer_matrices[o].get<1>());
      for (int d=dmin; d <= max_dist_for_deviation; ++d) {
	overlapping_models[r][o][d] = log2<FloatType>(my_cob_params[r].expected_overlapping_dimer_PWMs[o][d] + my_cob_params[r].deviation[o][d]);
      }
    }
  }
#pragma omp parallel for shared(use_two_strands) reduction(+:result) schedule(static)
  for (int i=0; i < sequences.size(); ++i) {
    const std::string& line = sequences[i];
    const std::string& line_rev = sequences_rev[i];
    double temp = 0.0;
    double sum_of_Zs=0.0;

    /////////////////////
    // Monomer models
    
    for (int k=0; k < p; ++k) {
      for (int j=0; j < fixed_m[i][k]; ++j) {
	feclearexcept(FE_ALL_EXCEPT);
	temp +=
	  (compute_log_probability<FloatType>(line, line_rev, j, 1, log_PWM[k], log_q, q2)
	   + log_lambda[k] - log(fixed_m[i][k]) - l2)
	  * fixed_Z[i][k][0][j];
	if (use_two_strands)
	  temp +=
	    (compute_log_probability<FloatType>(line, line_rev, j, -1, log_PWM[k], log_q_rev, q2_rev)
	     + log_lambda[k] - log(fixed_m[i][k]) - l2)
	    * fixed_Z[i][k][1][j];
	int retval = fetestexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW | FE_UNDERFLOW);
	if (retval) {
	  print_math_error(retval);
	  printf("Sequence %i, monomer model %i, position %i\n", i, k, j);
	}
	sum_of_Zs += fixed_Z[i][k][0][j];
	if (use_two_strands)
	  sum_of_Zs += fixed_Z[i][k][1][j];
      }
    }

    for (int r=0; r < number_of_cobs; ++r) {
      int tf1 = my_cob_params[r].tf1;
      int tf2 = my_cob_params[r].tf2;
      int number_of_orientations = tf1 == tf2 ? 3 : 4;
      int dmin = my_cob_params[r].dmin;
      int dmax = my_cob_params[r].dmax;
      int max_dist_for_deviation = my_cob_params[r].max_dist_for_deviation;
      for (int o=0; o < number_of_orientations; ++o) {
	for (int d=max_dist_for_deviation+1; d <= dmax; ++d) {       // spaced dimers
	  int m = my_cob_params[r].dimer_m[i][d];
	  if (m <= 0)
	    continue;
	  FloatType log_m = log2l(m);
	  for (int j=0; j < m; ++j) {
	    if (my_cob_params[r].dimer_lambdas[o][d] > 0.0) {
	      feclearexcept(FE_ALL_EXCEPT);
	      temp +=
		(compute_log_dimer_probability<FloatType>(line, line_rev,
							  j, +1,
							  log_oriented_dimer_matrices[r][o][0], 
							  log_oriented_dimer_matrices[r][o][0], 
							  d, log_q, q2)
		 + log2l(my_cob_params[r].dimer_lambdas[o][d]) - log_m - l2)
		* my_cob_params[r].spaced_dimer_Z[i][o][d][0][j];
	      if (use_two_strands)
		temp +=
		  (compute_log_dimer_probability<FloatType>(line, line_rev,
							    j, -1,
							    log_oriented_dimer_matrices[r][o][0], 
							    log_oriented_dimer_matrices[r][o][1], 
							    d, log_q_rev, q2_rev)
		   + log2l(my_cob_params[r].dimer_lambdas[o][d]) - log_m - l2)
		  * my_cob_params[r].spaced_dimer_Z[i][o][d][1][j];
	      int retval = fetestexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW | FE_UNDERFLOW);
	      if (retval) {
		print_math_error(retval);
		printf("Sequence %i, r=%i, spaced dimer model %s %i, position %i\n", i, r, orients[o], d, j);
	      }
	    }
	    sum_of_Zs += my_cob_params[r].spaced_dimer_Z[i][o][d][0][j];
	    if (use_two_strands)
	      sum_of_Zs += my_cob_params[r].spaced_dimer_Z[i][o][d][1][j];
	  } // end for j
	} // end for d

	for (int d=dmin; d <= max_dist_for_deviation; ++d) {        // overlapping dimers
	  int m = my_cob_params[r].dimer_m[i][d];
	  if (m <= 0)
	    continue;
	  FloatType log_m = log2l(m);
	  const matrix<FloatType>& model = overlapping_models[r][o][d];
	  for (int j=0; j < m; ++j) {
	    if (my_cob_params[r].dimer_lambdas[o][d] > 0.0) {
	      feclearexcept(FE_ALL_EXCEPT);
	      temp +=
		(compute_log_probability<FloatType>(line, line_rev, j, +1, model, log_q, q2)
		 + log2l(my_cob_params[r].dimer_lambdas[o][d]) - log_m - l2)
		* my_cob_params[r].overlapping_dimer_Z[i][o][d][0][j];
	      if (use_two_strands)
		temp +=
		  (compute_log_probability<FloatType>(line, line_rev, j, -1, model, log_q_rev, q2_rev)
		   + log2l(my_cob_params[r].dimer_lambdas[o][d]) - log_m - l2)
		  * my_cob_params[r].overlapping_dimer_Z[i][o][d][1][j];
	      int retval = fetestexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW | FE_UNDERFLOW);
	      if (retval) {
		print_math_error(retval);
		printf("Sequence %i, r=%i, overlapping dimer model %s %i, position %i\n", i, r, orients[o], d, j);
	      }
	    }
	    sum_of_Zs += my_cob_params[r].overlapping_dimer_Z[i][o][d][0][j];
	    if (use_two_strands)
	      sum_of_Zs += my_cob_params[r].overlapping_dimer_Z[i][o][d][1][j];
	  }
	}

      }  // end for o
    } // end for r

    feclearexcept(FE_ALL_EXCEPT);
    if (lambda_bg > 0.0)
      temp += (compute_log_background_probability<FloatType>(line, log_q, q2) + log_lambda_bg) * (1.0 - sum_of_Zs);
    int retval = fetestexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW | FE_UNDERFLOW);
    if (retval) {
      print_math_error(retval);
      printf("Sequence %i, background model\n", i);
    }
    result += temp;
  } // end for i

  return result;
}




////////////////////////////////////////////
//
// Helper functions to compute expectations
//
////////////////////////////////////////////

// This is for monomer cases

void
expectation_Z_dir_j(boost::multi_array<FloatType, 4>& Z, int i, int k,
		    const std::string& line,
		    const std::string& line_rev,
		    int m, FloatType log_lambda, const matrix<FloatType>& log_PWM, 
		    const std::vector<FloatType>& log_bg_model,
		    const std::vector<FloatType>& log_bg_model_rev,
		    const dmatrix& bg_model_markov, const dmatrix& bg_model_markov_rev)
{
  FloatType log_m = log2l((long double)m);
  FloatType log_lambda_per_pos = log_lambda - log_m;
  if (use_two_strands)
    log_lambda_per_pos -= log2l(2.0);
  for (int j=0; j < m; ++j) {
    FloatType p1 = compute_log_probability<FloatType>(line, line_rev, j, 1, log_PWM, log_bg_model, bg_model_markov);

    // f_j
    feclearexcept(FE_ALL_EXCEPT);
    Z[i][k][0][j]=p1 + log_lambda_per_pos;
    int retval = fetestexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW | FE_UNDERFLOW);
    if (retval) {
      print_math_error(retval);
      printf("Expectation, sequence %i, pos %i, monomer model %i\n", i, j, k);
    }
    if (use_two_strands) {
      FloatType p2 = compute_log_probability<FloatType>(line, line_rev, j, -1, log_PWM, log_bg_model_rev, bg_model_markov_rev);

      feclearexcept(FE_ALL_EXCEPT);
      Z[i][k][1][j]=p2 + log_lambda_per_pos;   // reverse complement
      int retval = fetestexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW | FE_UNDERFLOW);
      if (retval) {
	print_math_error(retval);
	printf("Expectation, sequence %i, pos %i, monomer model %i\n", i, j, k);
      }
    }
  }
}

//This is for overlapping dimer cases
void
expectation_Z_dir_j_overlapping(boost::multi_array<FloatType, 5>& Z, int i, int o, int d,
				const std::string& line,
				const std::string& line_rev,
				int m, FloatType log_lambda, const matrix<FloatType>& log_PWM, 
				const std::vector<FloatType>& log_bg_model,
				const std::vector<FloatType>& log_bg_model_rev,
				const dmatrix& bg_model_markov, const dmatrix& bg_model_markov_rev)
{
  
  if (not std::isfinite(log_lambda) or log_lambda == 0.0) { // the dimer does not fit on this line
    for (int j=0; j < m; ++j) {
      Z[i][o][d][0][j]=0.0;
      if (use_two_strands)
	Z[i][o][d][1][j]=0.0;
    }
    return;
  }
  
  FloatType log_m = log2l(m);
  FloatType log_lambda_per_pos = log_lambda - log_m;
  if (use_two_strands)
    log_lambda_per_pos -= log2l(2.0);
  for (int j=0; j < m; ++j) {
    FloatType p1 = compute_log_probability<FloatType>(line, line_rev, j, 1, log_PWM, log_bg_model, bg_model_markov);

    // f_j
    feclearexcept(FE_ALL_EXCEPT);
    Z[i][o][d][0][j]=p1 + log_lambda_per_pos;
    int retval = fetestexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW | FE_UNDERFLOW);
    if (retval) {
      print_math_error(retval);
      printf("Expectation, sequence %i, pos %i, overlapping model %s %i\n", i, j, orients[o], d);
    }

    if (use_two_strands) {
      FloatType p2 = compute_log_probability<FloatType>(line, line_rev, j, -1, log_PWM, log_bg_model_rev, bg_model_markov_rev);

      feclearexcept(FE_ALL_EXCEPT);
      Z[i][o][d][1][j]=p2 + log_lambda_per_pos;   // reverse complement
      int retval = fetestexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW | FE_UNDERFLOW);
      if (retval) {
        print_math_error(retval);
        printf("Expectation, sequence %i, pos %i, overlapping model %s %i\n", i, j, orients[o], d);
      }

    }
  }
}

//This is for spaced dimer cases
void
expectation_Z_dir_j_spaced(boost::multi_array<FloatType, 5>& Z, int i,
			   int o, int d,
			   const std::string& line,
			   const std::string& line_rev,
			   int m, FloatType log_lambda,
			   const matrix<FloatType>& log_PWM1,
			   const matrix<FloatType>& log_PWM2,
			   const std::vector<FloatType>& log_bg_model,
			   const std::vector<FloatType>& log_bg_model_rev,
			   const dmatrix& bg_model_markov, const dmatrix& bg_model_markov_rev)
{
  if (not std::isfinite(log_lambda) or log_lambda == 0.0) {
    for (int j=0; j < m; ++j) {
      Z[i][o][d][0][j]=0.0;
      if (use_two_strands)
	Z[i][o][d][1][j]=0.0;
    }
    return;
  }
  FloatType log_m = log2l(m);
  FloatType log_lambda_per_pos = log_lambda - log_m;
  if (use_two_strands)
    log_lambda_per_pos -= log2l(2.0);
  for (int j=0; j < m; ++j) {
    FloatType p1 = compute_log_dimer_probability<FloatType>(line, line_rev, j, 1, log_PWM1, log_PWM2, d, log_bg_model, bg_model_markov);


    // f_j
    feclearexcept(FE_ALL_EXCEPT);
    Z[i][o][d][0][j]=p1 + log_lambda_per_pos;
    int retval = fetestexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW | FE_UNDERFLOW);
    if (retval) {
      print_math_error(retval);
      printf("Expectation, sequence %i, pos %i, spaced model %s %i\n", i, j, orients[o], d);
    }
    if (use_two_strands) {
      FloatType p2 = compute_log_dimer_probability<FloatType>(line, line_rev, j, -1, log_PWM1, log_PWM2, d, log_bg_model_rev, bg_model_markov_rev);

      feclearexcept(FE_ALL_EXCEPT);
      Z[i][o][d][1][j]=p2 + log_lambda_per_pos;   // reverse complement
      int retval = fetestexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW | FE_UNDERFLOW);
      if (retval) {
        print_math_error(retval);
        printf("Expectation, sequence %i, pos %i, spaced model %s %i\n", i, j, orients[o], d);
      }
    }
  }
}


typedef std::priority_queue<FloatType, std::vector<FloatType>, std::greater<FloatType> > queue_type;
//typedef std::vector<FloatType> queue_type;

void
sum_Z_dir_j(const boost::multi_array<FloatType, 5>& Z, int i, int o, int d, int m, queue_type& queue)
{
  for (int j=0; j < m; ++j) {
    if (Z[i][o][d][0][j] != 0.0) {
      queue.push(Z[i][o][d][0][j]);
      //queue.push_back(Z[i][o][d][0][j]);
    }
  }
  if (use_two_strands) {
    for (int j=0; j < m; ++j) {
      if (Z[i][o][d][1][j] != 0.0) {
	queue.push(Z[i][o][d][1][j]);
	//queue.push_back(Z[i][o][d][1][j]);
      }
    }
  }
  return;
}

void
sum_Z_dir_j(const boost::multi_array<FloatType, 4>& Z, int i, int k, int m, queue_type& queue)
{
  for (int j=0; j < m; ++j) {
    if (Z[i][k][0][j] != 0.0) {
      queue.push(Z[i][k][0][j]);
      //queue.push_back(Z[i][k][0][j]);
    }
  }
  if (use_two_strands) {
    for (int j=0; j < m; ++j) {
      if (Z[i][k][1][j] != 0.0) {
	queue.push(Z[i][k][1][j]);
	//queue.push_back(Z[i][k][1][j]);
      }
    }
  }
  return;
}




void
normalize_Z_dir_j(boost::multi_array<FloatType, 4>& Z, int i, int k, int m, FloatType divisor)
{
  for (int j=0; j < m; ++j) {
    if (Z[i][k][0][j] != 0.0)
      Z[i][k][0][j] = exp2l(Z[i][k][0][j] - divisor);
    

    assert(Z[i][k][0][j] >= 0);
    assert(Z[i][k][0][j] <= 1);
  }

  if (use_two_strands) {
    for (int j=0; j < m; ++j) {
      if (Z[i][k][1][j] != 0.0)
	Z[i][k][1][j] = exp2l(Z[i][k][1][j] - divisor);
	  
      assert(Z[i][k][1][j] >= 0);
      assert(Z[i][k][1][j] <= 1);
    }
  }
}

void
normalize_Z_dir_j(boost::multi_array<FloatType, 5>& Z, int i, int o, int d, int m, FloatType divisor)
{
  for (int j=0; j < m; ++j) {
    if (Z[i][o][d][0][j] != 0.0)
      Z[i][o][d][0][j] = exp2l(Z[i][o][d][0][j] - divisor);

    assert(Z[i][o][d][0][j] >= 0);
    assert(Z[i][o][d][0][j] <= 1);
  }

  if (use_two_strands) {
    for (int j=0; j < m; ++j) {
      if (Z[i][o][d][1][j] != 0.0)
	Z[i][o][d][1][j] = exp2l(Z[i][o][d][1][j] - divisor);
	  
      assert(Z[i][o][d][1][j] >= 0);
      assert(Z[i][o][d][1][j] <= 1);
    }
  }
}


void
get_new_weights(int j1, int dir, double z, int w, const std::string& seed,
		const std::string& line_orig, const std::string& line_rev_orig,
		dmatrix& weights, std::vector<double>& pred_flank, std::vector<double>& succ_flank,
		bool force_multinomial, bool compute_flanks)
{
  assert(seed.length() == w);
  typedef myuint128 bitstring_t;
  const std::string& line = dir == 1 ? line_orig : line_rev_orig;
  int L = line.length();
  if (dir == -1)
    j1 = L - j1 - w;
  bitstring_t mismatches = iupac_mismatch_positions<bitstring_t>(line.substr(j1, w), seed);
  int hd = mypopcount(mismatches);
  bitstring_t positions = 0;
  if (not force_multinomial or hd < hamming_radius)
    positions = ~static_cast<bitstring_t>(0);  // update all
  else if (hd == hamming_radius)  // update only mismatch positions
    positions = mismatches;
  else
    positions = 0;   // update nothing
  //  printf("HD is %i %i %s %s\n", hd, mismatches, print_bitvector(mismatches).c_str(), print_bitvector(positions).c_str());
  bitstring_t mask = static_cast<bitstring_t>(1)<<(w-1);
  for (int pos=0; pos < w; ++pos, mask>>=1) {
    if (positions & mask)
      weights(to_int(line[j1+pos]), pos) += z; // update columns of pwm marked by bit vector positions
  }
  if (allow_extension and compute_flanks and positions == ~0) {     // If dir == -1, then the caller must swap pred_flank and succ_flank
    if (j1 != 0)
      pred_flank[to_int(line[j1-1])] += z;
    if (j1 + w < L)
      succ_flank[to_int(line[j1+w])] += z;
  }
}

void
get_new_spaced_dimer_weights(int j1, int dir, double z, int d,
			     int w1, int w2,
			     const std::string& seed1, const std::string& seed2,
			     const std::string& line_orig, const std::string& line_rev_orig,
			     dmatrix& weights1, dmatrix& weights2, 
			     bool force_multinomial)
{
  assert(seed1.length() == w1);
  assert(seed2.length() == w2);
  assert(z >= 0);
  assert(z < 1.0);
  assert(w1 == weights1.get_columns());
  assert(w2 == weights2.get_columns());
  const std::string& line = dir == 1 ? line_orig : line_rev_orig;
  int L = line.length();
  int dimer_len = w1 + d + w2;
  if (dir == -1)
    j1 = L - j1 - dimer_len;
  int j2 = j1 + d + w1;  // position of the second leg

  // first part
  code_t mismatches = iupac_mismatch_positions<code_t>(line.substr(j1, w1), seed1);
  int hd = mypopcount(mismatches);
  code_t positions = 0;
  if (not force_multinomial or hd < hamming_radius)
    positions = ~0;  // update all
  else if (hd == hamming_radius)  // update only mismatch positions
    positions = mismatches;
  else
    positions = 0;   // update nothing
  code_t mask = 1ull<<(w1-1);
  for (int pos=0; pos < w1; ++pos, mask>>=1) {
    if (positions & mask)
      weights1(to_int(line[j1+pos]), pos) += z; // update all columns
  }


  // second part
  code_t mismatches2 = iupac_mismatch_positions<code_t>(line.substr(j2, w2), seed2);
  int hd2 = mypopcount(mismatches2);
  code_t positions2 = 0;
  if (not force_multinomial or hd2 < hamming_radius)
    positions2 = ~0;  // update all
  else if (hd2 == hamming_radius)  // update only mismatch positions
    positions2 = mismatches2;
  else
    positions2 = 0;   // update nothing
  mask=1ull<<(w2-1);
  for (int pos=0; pos < w2; ++pos,mask>>=1) {
    if (positions2 & mask)
      weights2(to_int(line[j2+pos]), pos) += z; // update all columns
  }

  
}

void
get_new_gap_weights(int j1, int dir, double z, int d,
		    int w1, int w2,
		    const std::string& seed1, const std::string& seed2,
		    const std::string& line_orig, const std::string& line_rev_orig,
		    dmatrix& weights1,
		    bool force_multinomial)//, bool compute_flanks)
{
  assert(seed1.length() == w1);
  assert(seed2.length() == w2);
  assert(z >= 0);
  assert(z < 1.0);
  const std::string& line = dir == 1 ? line_orig : line_rev_orig;
  int L = line.length();
  int dimer_len = w1 + d + w2;
  if (dir == -1)
    j1 = L - j1 - dimer_len;

  // first part
  std::string seed = seed1 + std::string(d, 'N') + seed2;
  seed[w1-1] = 'N';  // flanks of the gap are also N
  seed[w1+d] = 'N';
  typedef myuint128 bitstring_t;
  assert(seed.length() == dimer_len);
  bitstring_t mismatches = iupac_mismatch_positions<bitstring_t>(line.substr(j1, dimer_len), seed);
  int hd = mypopcount(mismatches);
  bitstring_t positions = 0;
  if (not force_multinomial or hd <= hamming_radius)
    positions = ~static_cast<bitstring_t>(0);  // update all
  else
    positions = 0;   // update nothing
  int first = w1-1;  // last position of the first half-site
  int last = w1+d;   // first position of the second half-site
  bitstring_t mask = static_cast<bitstring_t>(1)<<(d+w2);
  for (int pos=first; pos <= last; ++pos, mask>>=1) {
    if (positions & mask)
      weights1(to_int(line[j1+pos]), pos) += z; // update all columns
  }

}




void
get_new_weights_with_flanks(int j1, int dir, double z, int w, int Lmax, const std::string& seed,
		const std::string& line_orig, const std::string& line_rev_orig,
		dmatrix& weights,
		bool force_multinomial, bool compute_flanks)
{
  assert(seed.length() == w);
  typedef myuint128 bitstring_t;
  const std::string& line = dir == 1 ? line_orig : line_rev_orig;
  int L = line.length();
  if (dir == -1)
    j1 = L - j1 - w;
  int motif_pos = Lmax - w;        // Motif pos inside matrix 'weights' which has width 2*Lmax-w
  int seq_pos = motif_pos - j1; // Position of sequence inside matrix 'weights'
  
  bitstring_t mismatches = iupac_mismatch_positions<bitstring_t>(line.substr(j1, w), seed);
  int hd = mypopcount(mismatches);
  //  assert(hamming_distance(line.substr(j1, w), seed) == hd);
  bitstring_t positions = 0;
  if (not force_multinomial or hd < hamming_radius)
    positions = ~static_cast<bitstring_t>(0);  // update all
  else if (hd == hamming_radius)  // update only mismatch positions
    positions = mismatches;
  else
    positions = 0;   // update nothing
  //  printf("HD is %i %i %s %s\n", hd, mismatches, print_bitvector(mismatches).c_str(), print_bitvector(positions).c_str());
  bitstring_t mask = static_cast<bitstring_t>(1)<<(w-1);
  for (int pos=0; pos < w; ++pos, mask>>=1) {
    if (positions & mask)
      weights(to_int(line[j1+pos]), motif_pos + pos) += z; // update columns of pwm marked by bit vector positions
  }

  if (compute_flanks and (not force_multinomial or hd < hamming_radius)) {
    // Left flank
    for (int i = 0; i < j1; ++i)
      weights(to_int(line[i]), seq_pos + i) += z;
    // Right flank
    for (int i = j1+w; i < L; ++i)
      weights(to_int(line[i]), seq_pos + i) += z;
  }
}



void
get_new_spaced_dimer_weights_with_flanks(int j1, int dir, double z, int d,
					 int w1, int w2, int Lmax,
		    const std::string& seed1, const std::string& seed2,
		    const std::string& line_orig, const std::string& line_rev_orig,
		    dmatrix& weights,
		    bool force_multinomial,
		    bool compute_flanks)
{
  assert(seed1.length() == w1);
  assert(seed2.length() == w2);
  assert(z >= 0);
  assert(z < 1.0);
  const std::string& line = dir == 1 ? line_orig : line_rev_orig;
  int L = line.length();
  int dimer_len = w1 + d + w2;
  if (dir == -1)
    j1 = L - j1 - dimer_len;
  int motif_pos = Lmax - dimer_len;        // Motif pos inside matrix 'weights' which has width 2*Lmax-w
  int seq_pos = motif_pos - j1; // Position of sequence inside matrix 'weights'

  // first part
  std::string seed = seed1 + std::string(d, 'N') + seed2;
  seed[w1-1] = 'N';  // flanks of the gap are also N
  seed[w1+d] = 'N';
  typedef myuint128 bitstring_t;
  assert(seed.length() == dimer_len);
  bitstring_t mismatches = iupac_mismatch_positions<bitstring_t>(line.substr(j1, dimer_len), seed);
  int hd = mypopcount(mismatches);
  bitstring_t positions = 0;
  //  bitstring_t gap_positions = ((static_cast<bitstring_t>(1) << d) - 1) << w2;
  if (not force_multinomial or hd < hamming_radius)
    positions = ~static_cast<bitstring_t>(0);  // update all, including the gap
  else if (hd == hamming_radius)  // update only mismatch positions
    positions = mismatches;
  else
    positions = 0;   // update nothing
  bitstring_t mask = static_cast<bitstring_t>(1)<<(dimer_len-1);
  for (int pos=0; pos < dimer_len; ++pos, mask>>=1) {
    if (positions & mask)
      weights(to_int(line[j1+pos]), motif_pos + pos) += z; // update all columns
  }
  if (compute_flanks and (not force_multinomial or hd < hamming_radius)) {
    // Left flank
    for (int i = 0; i < j1; ++i)
      weights(to_int(line[i]), seq_pos + i) += z;
    // Right flank
    for (int i = j1+dimer_len; i < L; ++i)
      weights(to_int(line[i]), seq_pos + i) += z;
  }
}


boost::tuple<std::vector<dmatrix>, std::vector<cob_of_matrices> >
get_models_with_flanks(const std::vector<std::string>& sequences,
		       const std::vector<std::string>& sequences_rev,
		       const std::vector<std::string>& fixed_seed,
		       const std::vector<int>& fixed_w,
		       int Lmax,
		       const boost::multi_array<int, 2>& fixed_m,
		       const std::vector<cob_params_t>& my_cob_params,
		       const array_4d_type& fixed_Z)
{
  typedef boost::multi_array<double, 2>::extent_range range;

  int number_of_cobs = my_cob_params.size();
  int fixed_p = fixed_seed.size();
  int lines = sequences.size();
  

  std::vector<dmatrix> flank_fixed_PWM;

  std::vector<cob_of_matrices> flank_dimer_PWM;
      


  // Initialize the models to zero

  // spaced models
  for (int r = 0; r < number_of_cobs; ++r) {
    int no = my_cob_params[r].number_of_orientations;
    //    int max_dist_for_deviation = my_cob_params[r].max_dist_for_deviation;
    int dmax = my_cob_params[r].dmax;
    flank_dimer_PWM.push_back(cob_of_matrices(boost::extents[no][range(my_cob_params[r].dmin,dmax+1)]));
  }

  // Fixed models
  for (int k=0; k < fixed_p; ++k) { 
    flank_fixed_PWM.push_back(dmatrix(4, 2*Lmax-fixed_w[k]));
  }


  // Overlapping dimer models
  for (int r = 0; r < number_of_cobs; ++r) {
    const int& dmin = my_cob_params[r].dmin;
    const int& dmax = my_cob_params[r].dmax;
    const int& no = my_cob_params[r].number_of_orientations;
    for (int o=0; o < no; ++o) {
      for (int d=dmin; d <= dmax; ++d) {
	  flank_dimer_PWM[r][o][d] = dmatrix(4, 2*Lmax - my_cob_params[r].dimer_w[d]);
      }
    }
  }

  int dirs[2] = {1, -1};
  int maxdir = use_two_strands ? 2 : 1;

  ////////////////////////////
  //
  // Monomer models
  //
  ////////////////////////////
  
  // Signal from monomeric models
  for (int k=0; k < fixed_p; ++k) {
    dmatrix pwm(4, 2*Lmax-fixed_w[k]);
    // This requires at least gcc 4.9.
    // clang (at least not 3.8) does not support 'declare reduction' even
    // though it defines _OPENMP to 201307 (that is openmp 4.0).
#if defined(_OPENMP) && (_OPENMP >= 201307) && !defined(__clang__)
#pragma omp declare reduction( + : dmatrix : omp_out+=omp_in ) initializer(omp_priv(omp_orig.dim()))
#pragma omp declare reduction( + : std::vector<double, std::allocator<double> > : omp_out+=omp_in ) initializer(omp_priv(std::vector<double>(omp_orig.size())))
#pragma omp parallel for shared(lines,sequences,use_two_strands) reduction(+:pwm) schedule(static)
#endif
    for (int i = 0; i < lines; ++i) {
      const std::string& line = sequences[i];
      const std::string& line_rev = sequences_rev[i];
      for (int j1=0; j1 < fixed_m[i][k]; ++j1) {  // iterates through start positions
	for (int dir=0; dir < maxdir; ++dir) {
	  double z = fixed_Z[i][k][dir][j1];
	  get_new_weights_with_flanks(j1, dirs[dir], z, fixed_w[k], Lmax, fixed_seed[k], line, line_rev, pwm, 
			use_multinomial, true);
	}
      }
    }  // end for lines
    flank_fixed_PWM[k] = pwm;
  } // for k, fixed PWMs


	
  ////////////////////////////
  //
  // Overlapping dimer
  //
  ////////////////////////////
	
  // #pragma omp parallel for schedule(static)
  for (int r = 0; r < number_of_cobs; ++r) {
    for (int o=0; o < my_cob_params[r].number_of_orientations; ++o) {
      for (int d=my_cob_params[r].dmin; d < std::min(0, my_cob_params[r].dmax+1); ++d) {
	if (my_cob_params[r].dimer_lambdas[o][d] == 0.0)
	  continue;

	dmatrix pwm(4, 2*Lmax - my_cob_params[r].dimer_w[d]);
	// This requires gcc 4.9
	// clang (at least not 3.8) does not support 'declare reduction' even
	// though it defines _OPENMP to 201307 (that is openmp 4.0).
#if defined(_OPENMP) && (_OPENMP >= 201307) && !defined(__clang__)
#pragma omp declare reduction( + : dmatrix : omp_out+=omp_in ) initializer(omp_priv(omp_orig.dim()))
#pragma omp declare reduction( + : std::vector<double, std::allocator<double> > : omp_out+=omp_in ) initializer(omp_priv(std::vector<double>(omp_orig.size())))	
#pragma omp parallel for reduction(+:pwm) schedule(static)
#endif
	for (int i = 0; i < lines; ++i) {
	  const std::string& line = sequences[i];
	  const std::string& line_rev = sequences_rev[i];

	  for (int j1=0; j1 < my_cob_params[r].dimer_m[i][d]; ++j1) {  // iterates through start positions
	    for (int dir=0; dir < maxdir; ++dir) {
	      double z = my_cob_params[r].overlapping_dimer_Z[i][o][d][dir][j1];
	      get_new_weights_with_flanks(j1, dirs[dir], z, my_cob_params[r].dimer_w[d], Lmax, my_cob_params[r].dimer_seeds[o][d], 
					  line, line_rev, pwm, use_multinomial, true);
	    }
	  } // for j1
	}  // for i
	flank_dimer_PWM[r][o][d] = pwm;
      } // for d, overlapping dimer PWMs
    } // for o, overlapping dimer PWMs

  } // end for r


      





      
  ////////////////////////////
  //
  // spaced dimer
  //
  ////////////////////////////


  for (int r = 0; r < number_of_cobs; ++r) {
    int tf1 = my_cob_params[r].tf1;
    int tf2 = my_cob_params[r].tf2;
    // temps for spaced
    int max_dist_for_deviation = my_cob_params[r].max_dist_for_deviation;
    for (int o=0; o < my_cob_params[r].number_of_orientations; ++o) {
      //      for (int d=0; d <= dmax; ++d) {
      for (int d=0; d <= max_dist_for_deviation; ++d) {
	int dimer_len = fixed_w[tf1] + d + fixed_w[tf2];
	dmatrix pwm(4, 2*Lmax - dimer_len);

	pwm.fill_with(0.0);

	// This requires gcc 4.9
    // clang (at least not 3.8) does not support 'declare reduction' even
    // though it defines _OPENMP to 201307 (that is openmp 4.0).
#if defined(_OPENMP) && (_OPENMP >= 201307) && !defined(__clang__)
#pragma omp declare reduction( + : dmatrix : omp_out+=omp_in ) initializer(omp_priv(omp_orig.dim()))
#pragma omp declare reduction( + : std::vector<double> : omp_out+=omp_in ) initializer(omp_priv(std::vector<double>(omp_orig.size())))
#pragma omp parallel for reduction(+:pwm) schedule(static)
#endif	    
	for (int i = 0; i < lines; ++i) {
	  const std::string& line = sequences[i];
	  const std::string& line_rev = sequences_rev[i];

	  for (int dir=0; dir < maxdir; ++dir) {

	    for (int j1=0; j1 < my_cob_params[r].dimer_m[i][d]; ++j1) {  // iterates through start positions
	      double z = my_cob_params[r].overlapping_dimer_Z[i][o][d][dir][j1];
	      get_new_spaced_dimer_weights_with_flanks(j1, dirs[dir], z, d, fixed_w[tf1], fixed_w[tf2], Lmax,
						       my_cob_params[r].oriented_dimer_seeds[o].get<0>(), 
						       my_cob_params[r].oriented_dimer_seeds[o].get<1>(),
						       line, line_rev, 
						       pwm,
						       use_multinomial, true);
	    } // for j1
	  }  // for dir
	  
	} // for i in lines
	flank_dimer_PWM[r][o][d] = pwm;
      } // for d, spaced dimer PWMs
    } // for o, spaced dimer PWMs

  } // for end r

	    
  for (int k=0; k < fixed_p; ++k) {
    if (use_pseudo_counts)
      pseudo_counts.add(flank_fixed_PWM[k]);
    normalize_matrix_columns(flank_fixed_PWM[k]);
    //    write_matrix(stdout, flank_fixed_PWM[k], to_string("Flank fixed matrix %i:\n", k), "%.6f");
  }
  
  for (int r = 0; r < number_of_cobs; ++r) {
    int dmin = my_cob_params[r].dmin;
    int dmax = my_cob_params[r].dmax;
    int no = my_cob_params[r].number_of_orientations;
    for (int o=0; o < no; ++o) {
      for (int d=dmin; d <= dmax; ++d) {
	if (my_cob_params[r].dimer_lambdas[o][d] == 0.0)
	  continue;
	if (use_pseudo_counts)
	  pseudo_counts.add(flank_dimer_PWM[r][o][d]);
	normalize_matrix_columns(flank_dimer_PWM[r][o][d]);
      }
    }
    
  }

  return boost::make_tuple(flank_fixed_PWM, flank_dimer_PWM);

} // end of get_models_with_flanks





double
reestimate_dimer_lambdas(const int lines,
			 const int number_of_orientations,
			 const int dmin, const int dmax,
			 const boost::multi_array<double, 2>& dimer_m,
			 const array_5d_type& dimer_Z,
			 boost::multi_array<double, 2>& dimer_lambdas
)
{
  double total_sum = 0.0;
  for (int o=0; o < number_of_orientations; ++o) {
    for (int d=dmin; d <= dmax; ++d) {
      if (dimer_lambdas[o][d] == 0.0)
	continue;
      double k_sum = 0.0;
#pragma omp parallel for shared(use_two_strands) reduction(+:k_sum) schedule(static)
      for (int i = 0; i < lines; ++i) {
	for (int j=0; j < dimer_m[i][d]; ++j) {
	  k_sum += dimer_Z[i][o][d][0][j];
	  if (use_two_strands)
	    k_sum += dimer_Z[i][o][d][1][j];
	}
      }  // end for i
      double temp = k_sum / lines;
      dimer_lambdas[o][d] = temp;
      total_sum += temp;
    } // end for d
  } // end for o

  return total_sum;
}

FloatType
reestimate_PWM_lambdas(int lines, int p, const boost::multi_array<int, 2>& m, const array_4d_type& Z, std::vector<double>& lambda)
{
  FloatType total_sum = 0.0;
  for (int k=0; k < p; ++k) {
    FloatType k_sum = 0.0;
#pragma omp parallel for shared(lines,m,use_two_strands,Z,k) reduction(+:k_sum) schedule(static)
    for (int i = 0; i < lines; ++i) {
      for (int j=0; j < m[i][k]; ++j) 
	k_sum += Z[i][k][0][j];
      if (use_two_strands)
	for (int j=0; j < m[i][k]; ++j) 
	  k_sum += Z[i][k][1][j];
    }  // end for i
    lambda[k] = k_sum / lines;
    total_sum += lambda[k];
  } // end for k

  return total_sum;
}

void
add_columns(std::vector<double>& v, const dmatrix& m)
{
  assert(v.size() == 4);
  for (int j=0; j < 4; ++j) { // for every row
    v[j] += sum(m.row(j));
    //    if (use_two_strands)
    //      v[4-j-1] += sum(m.row(j));
  }
}

std::string
create_overlapping_seed(const std::string& seed1, const std::string& seed2, int d, const gapped_kmer_context& my_gapped_kmer_context)
{
  int k1 = seed1.length();
  int k2 = seed2.length();

  //int my_min_flank = std::min(k1+d,k2+d);
  std::string pattern;
  int hd;
  // if (my_min_flank <= 3) {
  //   hd = hamming_distance_overlapping_seeds_OR;
  //   pattern = combine_seeds_OR(seed1, seed2, d);
  // }
  // else {
  hd = hamming_distance_overlapping_seeds_N;   // this is zero
  pattern = combine_seeds_func(seed1, seed2, d);
    //  }
  assert(pattern.size() == k1 + k2 + d);
  std::string overlapping_seed;
  if (maximize_overlapping_seeds) {
    int begin = k1 + d;
    int len = -d;
    overlapping_seed = my_gapped_kmer_context.max(pattern, hd, begin, len).get<0>();
  }
  else
    overlapping_seed = pattern;
  /*
    for (int i=0; i < k1+d; ++i)                       // Flanks should be the same
    assert(overlapping_seed[i] == pattern[i]);
    for (int i=k1; i < k1+k2+d; ++i)
    assert(overlapping_seed[i] == pattern[i]);
  */
  assert(overlapping_seed.size() == k1 + k2 + d);
  return  overlapping_seed;
}


// Here a sequence is considered as background if the posterior probability of background model
// is higher than the posterior probability of any other model
void
print_background_sequences(const std::vector<std::string>& sequences, 
			   const array_4d_type& fixed_Z, 
			   const std::vector<cob_params_t>& my_cob_params,
			   const std::vector<dmatrix>& fixed_PWM,
			   const std::vector<double>& bg_model,
			   const dmatrix& bg_model_markov)
{
  int lines = sequences.size();
  int L=sequences[0].length();
  int fixed_p = fixed_PWM.size();
  int number_of_cobs = my_cob_params.size();
  std::vector<int> fixed_w(fixed_p);
  std::vector<int> fixed_m(fixed_p);
  for (int k=0; k < fixed_p; ++k) {
    fixed_w[k] = fixed_PWM[k].get_columns();
    fixed_m[k] = L-fixed_w[k]+1;
  }
  FILE* fp = fopen(unbound.c_str(), "w");
  assert(fp != NULL);
  for (int i=0; i < lines; ++i) {
    const std::string& line = sequences[i];

    // First compute the proportion of background as the complement of the sum of other models
    double total_sum = 0;
    for (int k=0; k < fixed_p; ++k) {
      double p = 0.0;
      for (int dir=0; dir < 2; ++dir) {
	for (int j=0; j < fixed_m[k]; ++j) {
	  p += fixed_Z[i][k][dir][j];
	}
      }
      total_sum += p;
    }
    for (int r=0; r < number_of_cobs; ++r) {
      int number_of_orientations = my_cob_params[r].number_of_orientations;
      int dmin = my_cob_params[r].dmin;
      int dmax = my_cob_params[r].dmax;
      int max_dist_for_deviation = my_cob_params[r].max_dist_for_deviation;
      for (int o=0; o < number_of_orientations; ++o) {

	for (int d=dmin; d <= max_dist_for_deviation; ++d) {
	  double p = 0.0;
	  for (int dir=0; dir < 2; ++dir) {
	    for (int j=0; j < my_cob_params[r].dimer_m[i][d]; ++j) {
	      p += my_cob_params[r].overlapping_dimer_Z[i][o][d][dir][j];
	    }
	  }
	  total_sum += p;
	}  // end for d

	for (int d=max_dist_for_deviation+1; d <= dmax; ++d) {
	  double p = 0.0;
	  for (int dir=0; dir < 2; ++dir) {
	    for (int j=0; j < my_cob_params[r].dimer_m[i][d]; ++j) {
	      p += my_cob_params[r].spaced_dimer_Z[i][o][d][dir][j];
	    }
	  }
	  total_sum += p;
	}  // end for d

      }  // end for o
    }  // end for r



    double background = 1.0 - total_sum;
    // Is the posterior probability of any monomeric model higher than the posterior probability of the bg model
    for (int k=0; k < fixed_p; ++k) {
      double p = 0.0;
      for (int dir=0; dir < 2; ++dir) {
	for (int j=0; j < fixed_m[k]; ++j) {
	  p += fixed_Z[i][k][dir][j];
	}
      }
      if (p >= background)
	goto my_exit;
    }
    // Is the posterior probability of any dimeric model higher than the posterior probability of the bg model
    for (int r=0; r < number_of_cobs; ++r) {
      int number_of_orientations = my_cob_params[r].number_of_orientations;
      int dmin = my_cob_params[r].dmin;
      int dmax = my_cob_params[r].dmax;
      int max_dist_for_deviation = my_cob_params[r].max_dist_for_deviation;
      for (int o=0; o < number_of_orientations; ++o) {

	// dimeric models with deviation
	for (int d=dmin; d <= max_dist_for_deviation; ++d) {
	  double p = 0.0;
	  for (int dir=0; dir < 2; ++dir) {                                      // KORJAA TÄMÄ, USE TWO STRANDS
	    for (int j=0; j < my_cob_params[r].dimer_m[i][d]; ++j) {
	      p += my_cob_params[r].overlapping_dimer_Z[i][o][d][dir][j];
	    }
	  }
	  if (p >= background)
	    goto my_exit;
	}  // end for d

	// dimeric models without deviation
	for (int d=max_dist_for_deviation+1; d <= dmax; ++d) {
	  double p = 0.0;
	  for (int dir=0; dir < 2; ++dir) {
	    for (int j=0; j < my_cob_params[r].dimer_m[i][d]; ++j) {
	      p += my_cob_params[r].spaced_dimer_Z[i][o][d][dir][j];
	    }
	  }
	  if (p >= background)
	    goto my_exit;
	}  // end for d

      }  // end for o
    }  // end for r

    fprintf(fp, "%s\n", line.c_str());
    my_exit:
    ;
  } // end for i in lines
  fclose(fp);
}


int
get_number_of_parameters(std::vector<cob_params_t>& my_cob_params, std::vector<dmatrix> fixed_PWM)
{
  int fixed_p=fixed_PWM.size();                // Number of fixed models
  int number_of_parameters = 0;
  number_of_parameters += fixed_p;              // pwm lambdas
  for (int k=0; k < fixed_p; ++k) {
    int w = fixed_PWM[k].get_columns();
    number_of_parameters += 3*w;       // pwm parameters
  }
  number_of_parameters += 3;                    // background
  for (int r=0; r < my_cob_params.size(); ++r) {
    int number_of_orientations = my_cob_params[r].number_of_orientations;
    int dmin = my_cob_params[r].dmin;
    int dmax = my_cob_params[r].dmax;
    int max_dist_for_deviation = my_cob_params[r].max_dist_for_deviation;
    number_of_parameters += number_of_orientations * (dmax-dmin+1);   // lambda parameters
    for (int d=dmin; d <= max_dist_for_deviation; ++d) {
      for (int o=0; o < number_of_orientations; ++o) {
	if (my_cob_params[r].dimer_lambdas[o][d] != 0.0)
	  number_of_parameters += 4 * (2+abs(d));     // correction matrices
      }
    }
  }
  // Background lambda is not counted since it is the complement of the rest of the lambdas
  return number_of_parameters;
}

dmatrix
reverse_complement_markov_model(const dmatrix& m)
{
  assert(m.get_columns() == 4);
  assert(m.get_rows() == 4);
  dmatrix result(m.dim());   // same shape
  for (int i=0; i < 4; ++i) {
    for (int j=0; j < 4; ++j) {
      result(i, j) = m(3-j, 3-i);
    }
  }
  return result;
}

bool
is_almost_palindrome(const dmatrix& m)
{
  dmatrix m_rev = reverse_complement(m);
  return distance(m, m_rev) < 0.001;
}

// uses the zoops model (zero or one occurrence per sequence)
// the error rate lambda[p] equals the quantity 1-gamma of the paper
std::vector<dmatrix>
multi_profile_em_algorithm(const std::vector<std::string>& sequences,
			   std::vector<dmatrix> fixed_PWM,
			   const std::vector<bool>& keep_monomer_fixed,
			   dmatrix& bg_model_markov, 
			   std::vector<double> bg_model,
			   double background_lambda,
			   std::vector<double> fixed_lambda, std::vector<std::string> fixed_seed, 
			   std::vector<cob_params_t>& my_cob_params,
			   double epsilon, double extension_ic_threshold,
			   const gapped_kmer_context& my_gapped_kmer_context)
{
  printf("\nIn multi_profile_em_algorithm:\n");

  int number_of_cobs = my_cob_params.size();
  int directions = use_two_strands ? 2 : 1;

  int fixed_p=fixed_PWM.size();                // Number of fixed models

  assert(fixed_p == fixed_lambda.size());
  assert(fixed_p == keep_monomer_fixed.size());


  typedef BinOp<int>::type func_ptr;
  //typedef const int& (*func_ptr)(const int&, const int&);

  typedef boost::multi_array<dmatrix, 2> cob_of_matrices;

  typedef boost::multi_array<FloatType, 4> array_4d_type;  // for Z variable, indices are (i,k,dir,j)

  std::vector<std::string> sequences_rev(sequences.size());
  for (int i=0; i < sequences.size(); ++i)
    sequences_rev[i] = reverse_complement(sequences[i]);
  
  
  double first_maximum_log_likelihood = 0; //maximum_log_likelihood(fixed_PWM, bg_model, fixed_lambda, 
  //		       background_lambda, my_cob_params, sequences);
  double mll = first_maximum_log_likelihood;

  int lines = sequences.size();
  std::vector<int> fixed_w(fixed_p);
  boost::multi_array<int, 2> fixed_m(boost::extents[lines][fixed_p]);
  std::vector<int> L(lines);
  
  //  std::vector<int> fixed_m(fixed_p);

  
  std::vector<std::vector<double> > pred_flank(fixed_p);
  std::vector<std::vector<double> > succ_flank(fixed_p);
  std::vector<double> pred_ic(fixed_p);
  std::vector<double> succ_ic(fixed_p);


  std::vector<dmatrix> flank_fixed_PWM;
  std::vector<cob_of_matrices> flank_dimer_PWM;


  typedef boost::multi_array<double, 2>::extent_range range;

  std::vector<double> fixed_av_ic(fixed_p);
  const int max_extension_round = 5;
  int extension_round=0;
  int round;
  std::vector<int> iterations;
  while (true) {
    printf("===================================================\n");
    if (local_debug) {
      printf("Extension round %i\n", extension_round);
    }

    

    //////////////
    //
    // Fixed models
    //
    //////////////

    int Lmin = std::numeric_limits<int>::max();
    int Lmax = std::numeric_limits<int>::min();
    for (int i=0; i < lines; ++i) {
      L[i] = sequences[i].length();
      if (L[i] < Lmin)
	Lmin = L[i];
      if (L[i] > Lmax)
	Lmax = L[i];
    }

    myaccumulate<int> acc(0, static_cast<func_ptr>(std::max));   // The cast is needed because std::max is overloaded

    // Set lengths
    acc.reset(0);
    for (int k=0; k < fixed_p; ++k) {
      fixed_w[k] = fixed_PWM[k].get_columns();
      for (int i=0; i < lines; ++i) {
	fixed_m[i][k] = L[i]-fixed_w[k]+1;
	acc(fixed_m[i][k]);
      }
    }

    int fixed_m_max=acc.get();                     // maximum number of motif starting positions in a sequence



    array_4d_type fixed_Z(boost::extents[lines][fixed_p][directions][fixed_m_max]);
    
    ///////////////////////////
    //
    // Overlapping dimer models
    //
    ///////////////////////////
    
    for (int r=0; r < my_cob_params.size(); ++r) {   
      my_cob_params[r].set_lengths();
      my_cob_params[r].initialize_Z(lines);
    }


    std::vector<double> fixed_dist(fixed_p, DBL_MAX);   // Distances between matrices in consecutive rounds
                                                  // is the convergence criterion

    std::vector<boost::multi_array<double, 2> > deviation_dist;
    for (int r=0; r < my_cob_params.size(); ++r) {
      int number_of_orientations = my_cob_params[r].number_of_orientations;
      int dmin = my_cob_params[r].dmin;
      int max_dist_for_deviation = my_cob_params[r].max_dist_for_deviation;
      boost::multi_array<double, 2> temp(boost::extents[number_of_orientations][range(dmin, max_dist_for_deviation+1)]);
      //fill_with(temp, DBL_MAX); 
      deviation_dist.push_back(temp);
    }

    int number_of_parameters = get_number_of_parameters(my_cob_params, fixed_PWM);
    if (local_debug)
      printf("Total number of parameters is %i\n", number_of_parameters);


    for (int r=0; r < my_cob_params.size(); ++r) {
      my_cob_params[r].update_oriented_matrices(fixed_PWM, fixed_seed);
      my_cob_params[r].compute_expected_matrices(fixed_PWM);
    }


    //    std::vector<double> spaced_dimer_dist(spaced_dimer_p, DBL_MAX);   
    bool convergence_criterion_reached = false;
    for (round = 0; round < max_iter and not convergence_criterion_reached; ++round) {

      if (local_debug)
	printf("Round %i\n", round);


      // Print seeds
      printf("Fixed seeds are %s\n", print_vector(fixed_seed).c_str());
      if (use_multinomial and local_debug) {

	
	for (int r=0; r < my_cob_params.size(); ++r) {
	  const cob_params_t& cp = my_cob_params[r];
	  if (local_debug) 
	    print_cob(stdout, cp.dimer_seeds, 
		      to_string("Seeds of overlapping dimer cases %s:\n", cp.name().c_str()), "%s");

	  boost::multi_array<int, 2> seed_counts(get_shape(cp.dimer_seeds));
	  boost::multi_array<double, 2> exp_counts(get_shape(cp.dimer_seeds));

	  std::map<int, int> sites;
	  MY_FOREACH(d, cp.dimer_seeds[0]) {
	    int dimer_len = cp.dimer_w[d];
	    for (int i = 0; i < lines; ++i)
		sites[d] += L[i] - dimer_len + 1;
	  }
	  MY_FOREACH(o, cp.dimer_seeds) {
	    MY_FOREACH(d, cp.dimer_seeds[o]) {
	      int dimer_len = cp.dimer_w[d];
	      double p = pow(4, -dimer_len);
	      seed_counts[o][d] = my_gapped_kmer_context.count(cp.dimer_seeds[o][d]);
	      exp_counts[o][d] = directions * sites[d] * p;
	    }
	  }

	  if (local_debug) {
	    print_cob(stdout, seed_counts, to_string("Seed counts %s:\n", cp.name().c_str()), "%i");
	    print_cob(stdout, exp_counts, to_string("Exp seed counts %s:\n", cp.name().c_str()), "%.2f");
	    printf("\n");
	  }
	}
      }


      /////////////////////////////////
      //
      // compute expectations
      //
      /////////////////////////////////


      std::vector<double> bg_model_rev(bg_model.rbegin(), bg_model.rend()); // use this model when considering the reverse strand
      dmatrix bg_model_markov_rev = reverse_complement_markov_model(bg_model_markov);
      std::vector<FloatType> log_bg_model(4);
      std::vector<FloatType> log_bg_model_rev(4);
      
      FloatType log_background_lambda = log2l(background_lambda);

      std::vector<FloatType> log_fixed_lambda(fixed_p);
      std::vector<matrix<FloatType> > log_fixed_PWM;
      for (int i=0; i < 4; ++i) {
	log_bg_model[i] = log2l(bg_model[i]);
	log_bg_model_rev[i] = log2l(bg_model_rev[i]);
      }
      for (int i=0; i < fixed_p; ++i) {
	log_fixed_lambda[i] = log2l(fixed_lambda[i]);
	log_fixed_PWM.push_back(log2<FloatType>(fixed_PWM[i]));
      }

      int mindmin = std::numeric_limits<int>::max();
      int maxdmax = std::numeric_limits<int>::min();
      for (int r=0; r < number_of_cobs; ++r) {
	mindmin = std::min(mindmin, my_cob_params[r].dmin);
	maxdmax = std::max(maxdmax, my_cob_params[r].dmax);
      }
      
      typedef boost::multi_array<matrix<FloatType>, 2> cob_of_matrices_t;
      typedef cob_of_matrices_t::extent_range range;
      std::vector<cob_of_matrices_t> log_overlapping_models;
      boost::multi_array<matrix<FloatType>, 3> log_oriented_dimer_matrices(boost::extents[number_of_cobs][4][2]);
      boost::multi_array<FloatType, 3> log_dimer_lambda(boost::extents[number_of_cobs][4][range(mindmin,maxdmax+1)]);
      for (int r=0; r < number_of_cobs; ++r) {
	int tf1 = my_cob_params[r].tf1;
	int tf2 = my_cob_params[r].tf2;
	int number_of_orientations = tf1 == tf2 ? 3 : 4;
	int dmin = my_cob_params[r].dmin;
	int dmax = my_cob_params[r].dmax;
	int max_dist_for_deviation = my_cob_params[r].max_dist_for_deviation;
    
	log_overlapping_models.push_back(cob_of_matrices_t(boost::extents[number_of_orientations][range(dmin, max_dist_for_deviation+1)]));
	for (int o=0; o < number_of_orientations; ++o) {
	  log_oriented_dimer_matrices[r][o][0] = log2<FloatType>(my_cob_params[r].oriented_dimer_matrices[o].get<0>());
	  log_oriented_dimer_matrices[r][o][1] = log2<FloatType>(my_cob_params[r].oriented_dimer_matrices[o].get<1>());
	  for (int d=dmin; d <= max_dist_for_deviation; ++d) {
	    log_overlapping_models[r][o][d] = log2<FloatType>(my_cob_params[r].expected_overlapping_dimer_PWMs[o][d] + my_cob_params[r].deviation[o][d]);
	  }
	  for (int d=dmin; d <= dmax; ++d)
	    if (my_cob_params[r].dimer_lambdas[o][d] != 0.0)
	      log_dimer_lambda[r][o][d] = log2l(my_cob_params[r].dimer_lambdas[o][d]);
	    else
	      log_dimer_lambda[r][o][d] = 0.0;
	}
      }

#pragma omp parallel for shared(lines,sequences,bg_model,bg_model_markov,use_two_strands) schedule(static)
      for (int i=0; i < lines; ++i) {
	const std::string& line = sequences[i];
	const std::string& line_rev = sequences_rev[i];
	//printf("Processing sequence %i\n", i);

	// f_0
	FloatType log_background = compute_log_background_probability<FloatType>(line, log_bg_model, bg_model_markov);
	log_background += log_background_lambda;

	// Monomer models
	for (int k=0; k < fixed_p; ++k)
	  expectation_Z_dir_j(fixed_Z, i, k, line, line_rev, fixed_m[i][k], log_fixed_lambda[k], log_fixed_PWM[k], 
			      log_bg_model, log_bg_model_rev, bg_model_markov, bg_model_markov_rev);

	// Dimer models
	for (int r=0; r < my_cob_params.size(); ++r) {
	  int max_dist_for_deviation = my_cob_params[r].max_dist_for_deviation;
	  // Overlapping dimer models
	  for (int o=0; o < my_cob_params[r].number_of_orientations; ++o) {
	    for (int d=my_cob_params[r].dmin; d <= max_dist_for_deviation; ++d) {
	      //const dmatrix& model = my_cob_params[r].overlapping_dimer_PWM[o][d];
	      //const dmatrix& expected = my_cob_params[r].expected_overlapping_dimer_PWMs[o][d];
	      //const dmatrix& deviation = my_cob_params[r].deviation[o][d];
	      //const dmatrix& model = expected + deviation;
	      expectation_Z_dir_j_overlapping(my_cob_params[r].overlapping_dimer_Z, i, o, d, line, line_rev, 
					      my_cob_params[r].dimer_m[i][d], log_dimer_lambda[r][o][d],
					      log_overlapping_models[r][o][d],
					      log_bg_model, log_bg_model_rev, 
					      bg_model_markov, bg_model_markov_rev);
	    }
	  }
	  // Spaced dimer models
	  for (int o=0; o < my_cob_params[r].number_of_orientations; ++o)
	    for (int d=max_dist_for_deviation+1; d <= my_cob_params[r].dmax; ++d){
	      expectation_Z_dir_j_spaced(my_cob_params[r].spaced_dimer_Z, i, o, d, line, line_rev,
					 my_cob_params[r].dimer_m[i][d], log_dimer_lambda[r][o][d],
					 log_oriented_dimer_matrices[r][o][d],
					 log_oriented_dimer_matrices[r][o][d],
					 log_bg_model, log_bg_model_rev, 
					 bg_model_markov, bg_model_markov_rev);
	    }
	} // end for r

	// compute the sum
	
	std::priority_queue<FloatType, std::vector<FloatType>, std::greater<FloatType> > queue;
	//std::vector<FloatType> queue;
	
	queue.push(log_background);
	//queue.push_back(log_background);
	
	// Fixed models
	for (int k=0; k < fixed_p; ++k) 
	  sum_Z_dir_j(fixed_Z, i, k, fixed_m[i][k], queue);

	for (int r=0; r < my_cob_params.size(); ++r) {
	  int max_dist_for_deviation = my_cob_params[r].max_dist_for_deviation;
	  // Overlapping dimer models
	  for (int o=0; o < my_cob_params[r].number_of_orientations; ++o)
	    for (int d=my_cob_params[r].dmin; d <= max_dist_for_deviation; ++d)
	      sum_Z_dir_j(my_cob_params[r].overlapping_dimer_Z, i, o, d, my_cob_params[r].dimer_m[i][d], queue);

	  // Spaced dimer models
	  for (int o=0; o < my_cob_params[r].number_of_orientations; ++o)
	    for (int d=max_dist_for_deviation+1; d <= my_cob_params[r].dmax; ++d)
	      sum_Z_dir_j(my_cob_params[r].spaced_dimer_Z, i, o, d, my_cob_params[r].dimer_m[i][d], queue);
	}

	FloatType normalizing_constant = log_sum(queue);
	//	FloatType normalizing_constant = fixed_sum + overlapping_sum + spaced_sum + background;
	//	printf("Dividing term is %e\n", sum);
	

	// normalize
	
	// Fixed models
	for (int k=0; k < fixed_p; ++k) 
	  normalize_Z_dir_j(fixed_Z, i, k, fixed_m[i][k], normalizing_constant);

	for (int r=0; r < my_cob_params.size(); ++r) {
	  int max_dist_for_deviation = my_cob_params[r].max_dist_for_deviation;
	  // Overlapping dimer models
	  for (int o=0; o < my_cob_params[r].number_of_orientations; ++o)
	    for (int d=my_cob_params[r].dmin; d <= max_dist_for_deviation; ++d)
	      normalize_Z_dir_j(my_cob_params[r].overlapping_dimer_Z, i, o, d, my_cob_params[r].dimer_m[i][d], normalizing_constant);

	  // Spaced dimer models
	  for (int o=0; o < my_cob_params[r].number_of_orientations; ++o)
	    for (int d=1+max_dist_for_deviation; d <= my_cob_params[r].dmax; ++d)
	      normalize_Z_dir_j(my_cob_params[r].spaced_dimer_Z, i, o, d, my_cob_params[r].dimer_m[i][d], normalizing_constant);
	}

      } // for i in lines




      ////////////////////////
      //
      // do the maximization
      //
      ////////////////////////

    
      // re-estimate PWM models

      // Initialize weights, pred_flank, succ_flank
      
      // These are used just to compute the background by subtracting signal from the whole data
      std::vector<double> fixed_signal_sum(4, 0.0);   // These are used to learn new PWM models
      std::vector<double> fixed2_signal_sum(4, 0.0);  // These are not used to learn new PWM models
      std::vector<double> overlapping_signal_sum(4, 0.0);
      std::vector<double> gap_signal_sum(4, 0.0);     // For gap in spaced dimers with d\in[0,max_dist_for_deviation]
      
      dmatrix dinucleotide_signal(4, 4);

      std::vector<dmatrix> new_fixed_PWM;
      std::vector<dmatrix> new_fixed_PWM2;   // These are not used to learn new PWMs

      
      std::vector<cob_of_matrices> overlapping_dimer_weights;
      std::vector<cob_of_matrices> gap_weights;  // covers the area of [0,max_dist_for_deviation] + 1 flank on each side
      std::vector<cob_of_matrices> new_deviation;
      std::vector<boost::multi_array<bool, 2> > empty_gap;
      
      std::vector<boost::tuple<dmatrix,dmatrix> > spaced_dimer_weights_sum;
      std::vector<boost::tuple<dmatrix,dmatrix> > spaced_dimer_weights2_sum;



      // Initialize the models to zero


      // Fixed models
      for (int k=0; k < fixed_p; ++k) { 
	new_fixed_PWM.push_back(dmatrix(4, fixed_w[k]));
	new_fixed_PWM2.push_back(dmatrix(4, fixed_w[k]));
	if (allow_extension) {
	  pred_flank[k].assign(4, 0.0);
	  succ_flank[k].assign(4, 0.0);
	}
      }

      // spaced models
      for (int r = 0; r < number_of_cobs; ++r) {
	int width1 = fixed_w[my_cob_params[r].tf1];
	int width2 = fixed_w[my_cob_params[r].tf2];
	spaced_dimer_weights_sum.push_back(boost::make_tuple(dmatrix(4, width1), 
							     dmatrix(4, width2)));
	spaced_dimer_weights2_sum.push_back(boost::make_tuple(dmatrix(4, width1), 
							      dmatrix(4, width2)));
      }

      // Overlapping dimer models
      for (int r = 0; r < number_of_cobs; ++r) {
	const int& dmin = my_cob_params[r].dmin;
	int max_dist_for_deviation = my_cob_params[r].max_dist_for_deviation;
	const int& no = my_cob_params[r].number_of_orientations;
	const int& k1 = my_cob_params[r].k1;
	const int& k2 = my_cob_params[r].k2;
	overlapping_dimer_weights.push_back(cob_of_matrices(boost::extents[no][range(my_cob_params[r].dmin,max_dist_for_deviation+1)]));
	gap_weights.push_back(cob_of_matrices(boost::extents[no][range(0, max_dist_for_deviation+1)]));
	new_deviation.push_back(cob_of_matrices(boost::extents[no][range(dmin,max_dist_for_deviation+1)]));
	
	empty_gap.push_back(boost::multi_array<bool, 2>(boost::extents[no][range(dmin, max_dist_for_deviation+1)]));
	for (int o=0; o < no; ++o) {
	  for (int d=dmin; d <= max_dist_for_deviation; ++d) {
	    empty_gap[r][o][d] = false;
	    new_deviation[r][o][d] = dmatrix(4, k1 + k2 + d);  // This contains redundant flanks. Just to ease things
	    if (d < 0)
	      overlapping_dimer_weights[r][o][d] = dmatrix(4, my_cob_params[r].dimer_w[d]);
	    else
	      gap_weights[r][o][d] = dmatrix(4, my_cob_params[r].dimer_w[d]);
	  }
	}
      }




      ////////////////////////////////////
      //
      // re-estimate mixing parameters
      //
      ////////////////////////////////////

      double total_sum=0.0;

      // Fixed models
      total_sum += reestimate_PWM_lambdas(lines, fixed_p, fixed_m, fixed_Z, fixed_lambda);  // modifies fixed_lambda


      for (int r = 0; r < number_of_cobs; ++r) {	
	// Overlapping dimers
	total_sum += reestimate_dimer_lambdas(lines,
					      my_cob_params[r].number_of_orientations,
					      my_cob_params[r].dmin, my_cob_params[r].max_dist_for_deviation,
					      my_cob_params[r].dimer_m,
					      my_cob_params[r].overlapping_dimer_Z,
					      my_cob_params[r].dimer_lambdas);              // Modifies dimer_lambdas


	// Spaced dimers
	total_sum += reestimate_dimer_lambdas(lines,
					      my_cob_params[r].number_of_orientations,
					      my_cob_params[r].max_dist_for_deviation+1, my_cob_params[r].dmax,
					      my_cob_params[r].dimer_m,
					      my_cob_params[r].spaced_dimer_Z,
					      my_cob_params[r].dimer_lambdas);                             // Modifies dimer_lambdas
      } // end for r

      assert(total_sum - 1.0 < 0.000001);
      
      background_lambda = 1.0 - total_sum;    // weight for background model

      assert(background_lambda <= 1.0);

      if (background_lambda < 0) {   // to correct rounding errors
	assert(fabs(background_lambda - 0.0) < 0.001);
	background_lambda=0.0;
      }
	

      //minimum_distance_for_learning 

      std::vector<bool> is_fixed_pwm_part_of_cob(fixed_p, false);
      for (int r=0; r < my_cob_params.size(); ++r) {
	using boost::multi_array_types::index_range;
	int dmax = my_cob_params[r].dmax;
	if (dmax >= minimum_distance_for_learning) {
	  boost::multi_array<double, 2> subarray =
	    my_cob_params[r].dimer_lambdas[ boost::indices[index_range()][(long int)minimum_distance_for_learning <= index_range()] ];
	  if (sum(subarray) >= learning_fraction) {
	    is_fixed_pwm_part_of_cob[my_cob_params[r].tf1] = true;
	    is_fixed_pwm_part_of_cob[my_cob_params[r].tf2] = true;
	  }
	}
      }
      printf("Is monomer pwm learnt purely modularly: %s\n", print_vector(is_fixed_pwm_part_of_cob).c_str());




      ///////////////////////////////////////////////////////////////////////////////////
      //
      // Reestimate models
      //
      ///////////////////////////////////////////////////////////////////////////////////


      std::vector<double> dummy(4, 0.0);   // For flanks that are not really computed

      
      ////////////////////////////
      //
      // Monomeric models
      //
      ////////////////////////////
      int dirs[2] = {1, -1};
      int maxdir = use_two_strands ? 2 : 1;
      
      // Signal from monomeric models
      for (int k=0; k < fixed_p; ++k) {
	dmatrix pwm(4, fixed_w[k]);
	std::vector<double>& signal = is_fixed_pwm_part_of_cob[k] ? fixed2_signal_sum : fixed_signal_sum;

	dmatrix temp_dinucleotide_signal(4, 4);
	std::vector<double> temp_signal(4, 0.0);

	int w = fixed_w[k];
	// This requires at least gcc 4.9
	// clang (at least not 3.8) does not support 'declare reduction' even
	// though it defines _OPENMP to 201307 (that is openmp 4.0).
#if defined(_OPENMP) && (_OPENMP >= 201307) && !defined(__clang__)
#pragma omp declare reduction( + : dmatrix : omp_out+=omp_in ) initializer(omp_priv(omp_orig.dim()))
#pragma omp declare reduction( + : std::vector<double, std::allocator<double> > : omp_out+=omp_in ) initializer(omp_priv(std::vector<double>(omp_orig.size())))
#pragma omp parallel for shared(lines,sequences,use_two_strands) reduction(+:pwm,temp_dinucleotide_signal,temp_signal) schedule(static)
#endif
	for (int i = 0; i < lines; ++i) {
	  const std::string& line = sequences[i];
	  const std::string& line_rev = sequences_rev[i];
	  for (int j1=0; j1 < fixed_m[i][k]; ++j1) {  // iterates through start positions
	    for (int dir=0; dir < maxdir; ++dir) {
	      double z = fixed_Z[i][k][dir][j1];
	      get_new_weights(j1, dirs[dir], z, fixed_w[k], fixed_seed[k], line, line_rev, pwm, 
			      dummy, dummy, use_multinomial, false);

	      for (int pos=0; pos < w; ++pos)
		temp_signal[to_int(line[j1+pos])] += z;
	      
	      if (use_markov_background) {
		for (int pos2=0; pos2 < w-1; ++pos2)
		  temp_dinucleotide_signal(to_int(line[j1+pos2]), to_int(line[j1+pos2+1])) += z;
	      }
	    }
	  }
	}  // end for lines
	signal += temp_signal;
	dinucleotide_signal += temp_dinucleotide_signal;

	if (is_fixed_pwm_part_of_cob[k])
	  new_fixed_PWM2[k] += pwm;
	else
	  new_fixed_PWM[k] += pwm;

      } // for k, fixed PWMs


	
      ////////////////////////////
      //
      // Overlapping dimer
      //
      ////////////////////////////
	
      for (int r = 0; r < number_of_cobs; ++r) {
	std::vector<double>& signal = overlapping_signal_sum;
	
	for (int o=0; o < my_cob_params[r].number_of_orientations; ++o) {
	  for (int d=my_cob_params[r].dmin; d < std::min(0, my_cob_params[r].dmax+1); ++d) {
	    if (my_cob_params[r].dimer_lambdas[o][d] == 0.0)
	      continue;

	    bool local_use_two_strands;
	    if (avoid_palindromes and use_two_strands and (o == HH or o==TT) and
		is_almost_palindrome(my_cob_params[r].deviation[o][d])) {
	      local_use_two_strands = false;
	      printf("Avoiding palindrome %s %d iteration %i\n", orients[o], d, round);
	    }
	    else
	      local_use_two_strands = use_two_strands;
	    
	    int w = my_cob_params[r].dimer_w[d];
	    dmatrix temp_dinucleotide_signal(4, 4);
	    std::vector<double> temp_signal(4, 0.0);

	    dmatrix pwm(overlapping_dimer_weights[r][o][d].dim());
	    // This requires gcc 4.9
	    // clang (at least not 3.8) does not support 'declare reduction' even
	    // though it defines _OPENMP to 201307 (that is openmp 4.0).
#if defined(_OPENMP) && (_OPENMP >= 201307) && !defined(__clang__)
#pragma omp declare reduction( + : dmatrix : omp_out+=omp_in ) initializer(omp_priv(omp_orig.dim()))
#pragma omp declare reduction( + : std::vector<double, std::allocator<double> > : omp_out+=omp_in ) initializer(omp_priv(std::vector<double>(omp_orig.size())))	
#pragma omp parallel for reduction(+:pwm,temp_dinucleotide_signal,temp_signal) schedule(static)
#endif
	    for (int i = 0; i < lines; ++i) {
	      const std::string& line = sequences[i];
	      const std::string& line_rev = sequences_rev[i];

	      for (int j1=0; j1 < my_cob_params[r].dimer_m[i][d]; ++j1) {  // iterates through start positions
		for (int dir=0; dir < maxdir; ++dir) {
		  double z = my_cob_params[r].overlapping_dimer_Z[i][o][d][dir][j1];
		  get_new_weights(j1, dirs[dir], z, my_cob_params[r].dimer_w[d], my_cob_params[r].dimer_seeds[o][d], 
				  line, line_rev, pwm, dummy, dummy, use_multinomial, false);

		  for (int pos=0; pos < w; ++pos)
		    temp_signal[to_int(line[j1+pos])] += z;
	      
		  if (use_markov_background) {
		    for (int pos2=0; pos2 < w-1; ++pos2)
		      temp_dinucleotide_signal(to_int(line[j1+pos2]), to_int(line[j1+pos2+1])) += z;
		  }
		} // for dir
	      } // for j1
	    }  // for i
	    signal += temp_signal;
	    dinucleotide_signal += temp_dinucleotide_signal;
	    overlapping_dimer_weights[r][o][d] = pwm;
	  } // for d, overlapping dimer PWMs
	} // for o, overlapping dimer PWMs

      } // end for r


      
	////////////////////////////
	//
	// Spaced dimer
	//
	////////////////////////////


      for (int r = 0; r < number_of_cobs; ++r) {
	int tf1 = my_cob_params[r].tf1;
	int tf2 = my_cob_params[r].tf2;
	// temps for spaced
	dmatrix m1(4, fixed_w[tf1]);
	dmatrix m2(4, fixed_w[tf2]);

	for (int o=0; o < my_cob_params[r].number_of_orientations; ++o) {
	  for (int d=0; d <= my_cob_params[r].dmax; ++d) {

	    m1.fill_with(0.0);
	    m2.fill_with(0.0);

	    int w1 = fixed_w[tf1];
	    int w2 = fixed_w[tf2];
	    std::vector<double>& signal = d >= minimum_distance_for_learning ? fixed_signal_sum : fixed2_signal_sum;
	    
	    std::vector<double> temp_signal(4, 0.0);
	    dmatrix temp_dinucleotide_signal(4, 4);
	    // This requires gcc 4.9
	    // clang (at least not 3.8) does not support 'declare reduction' even
	    // though it defines _OPENMP to 201307 (that is openmp 4.0).
#if defined(_OPENMP) && (_OPENMP >= 201307) && !defined(__clang__)
#pragma omp declare reduction( + : dmatrix : omp_out+=omp_in ) initializer(omp_priv(omp_orig.dim()))
#pragma omp declare reduction( + : std::vector<double> : omp_out+=omp_in ) initializer(omp_priv(std::vector<double>(omp_orig.size())))
#pragma omp parallel for reduction(+:m1,m2, temp_dinucleotide_signal,temp_signal) schedule(static)
#endif	    
	    for (int i = 0; i < lines; ++i) {
	      const std::string& line = sequences[i];
	      const std::string& line_rev = sequences_rev[i];


	      for (int dir=0; dir < maxdir; ++dir) {
		int first = dir == 0 ? 0 : w2 + d;
		int second = dir == 0 ? w1+d : 0;
		for (int j1=0; j1 < my_cob_params[r].dimer_m[i][d]; ++j1) {  // iterates through start positions

		  double z = d <= my_cob_params[r].max_dist_for_deviation ?
		    my_cob_params[r].overlapping_dimer_Z[i][o][d][dir][j1] : 
		    my_cob_params[r].spaced_dimer_Z[i][o][d][dir][j1];
		  get_new_spaced_dimer_weights(j1, dirs[dir], z, d, fixed_w[tf1], fixed_w[tf2], 
					       my_cob_params[r].oriented_dimer_seeds[o].get<0>(), 
					       my_cob_params[r].oriented_dimer_seeds[o].get<1>(),
					       line, line_rev, 
					       m1, m2,
					       use_multinomial);

		  for (int pos=0; pos < w1; ++pos)
		    temp_signal[to_int(line[j1+first+pos])] += z;

		  for (int pos=0; pos < w2; ++pos)
		    temp_signal[to_int(line[j1+second+pos])] += z;

  
		  if (use_markov_background) {
		    // first part
		    for (int pos2=0; pos2 < w1-1; ++pos2)
		      temp_dinucleotide_signal(to_int(line[j1+first+pos2]), to_int(line[j1+first+pos2+1])) += z;
		    // second part
		    for (int pos2=0; pos2 < w2-1; ++pos2)
		      temp_dinucleotide_signal(to_int(line[j1+second+pos2]), to_int(line[j1+second+pos2+1])) += z;
		  }

		}
	      }

	    } // for i in lines
	    signal += temp_signal;
	    dinucleotide_signal += temp_dinucleotide_signal;
	    boost::tuple<dmatrix,dmatrix>& temp = 
	      d >= minimum_distance_for_learning ? spaced_dimer_weights_sum[r] : spaced_dimer_weights2_sum[r];

	    switch (o) {
	    case HT:
	      temp.get<0>() += m1;
	      temp.get<1>() += m2;
	      break;
	    case HH:
	      temp.get<0>() += m1;
	      temp.get<1>() += reverse_complement(m2);
	      break;
	    case TT:
	      temp.get<0>() += reverse_complement(m1);
	      temp.get<1>() += m2;
	      break;
	    case TH:
	      temp.get<0>() += reverse_complement(m1);
	      temp.get<1>() += reverse_complement(m2);
	      break;
	    }	      

	  } // for d, spaced dimer PWMs
	} // for o, spaced dimer PWMs

      } // for end r


      // Learn the 'fixed' models from spaced dimers
      

      for (int r = 0; r < number_of_cobs; ++r) {      // Combine all the dimer half sites into fixed PWMs
	// from spaced models with minimum_distance_for_learning >= d
	new_fixed_PWM[my_cob_params[r].tf1] += spaced_dimer_weights_sum[r].get<0>();     
	new_fixed_PWM[my_cob_params[r].tf2] += spaced_dimer_weights_sum[r].get<1>();

	// from spaced models with 0 <= d < minimum_distance_for_learning, and from monomer models
	new_fixed_PWM2[my_cob_params[r].tf1] += spaced_dimer_weights2_sum[r].get<0>();   
	new_fixed_PWM2[my_cob_params[r].tf2] += spaced_dimer_weights2_sum[r].get<1>();
      } // end for r



      
      	////////////////////////////
	//
	// gap
	//
	////////////////////////////


      for (int r = 0; r < number_of_cobs; ++r) {
	int tf1 = my_cob_params[r].tf1;
	int tf2 = my_cob_params[r].tf2;
	int max_dist_for_deviation = my_cob_params[r].max_dist_for_deviation;

	for (int o=0; o < my_cob_params[r].number_of_orientations; ++o) {
	  for (int d=0; d <= max_dist_for_deviation; ++d) {
	    if (my_cob_params[r].dimer_lambdas[o][d] == 0.0)
	      continue;
	    
	    int dimer_len = fixed_w[tf1] + d + fixed_w[tf2];
	    dmatrix m1(4, dimer_len);

	    m1.fill_with(0.0);
	    bool local_use_two_strands;
	    if (avoid_palindromes and use_two_strands and (o == HH or o==TT) and
		is_almost_palindrome(my_cob_params[r].deviation[o][d])) {
	      local_use_two_strands = false;
	      printf("Avoiding palindrome %s %d iteration %i\n", orients[o], d, round);
	    }
	    else
	      local_use_two_strands = use_two_strands;


	    int w1 = fixed_w[tf1];
	    int w2 = fixed_w[tf2];
	    std::vector<double> temp_signal(4, 0.0);
	    dmatrix temp_dinucleotide_signal(4, 4);
	    // This requires gcc 4.9
	    // clang (at least not 3.8) does not support 'declare reduction' even
	    // though it defines _OPENMP to 201307 (that is openmp 4.0).
#if defined(_OPENMP) && (_OPENMP >= 201307) && !defined(__clang__)
#pragma omp declare reduction( + : dmatrix : omp_out+=omp_in ) initializer(omp_priv(omp_orig.dim()))
#pragma omp declare reduction( + : std::vector<double> : omp_out+=omp_in ) initializer(omp_priv(std::vector<double>(omp_orig.size())))
#pragma omp parallel for reduction(+:m1, temp_dinucleotide_signal, temp_signal) schedule(static)
#endif	    
	    for (int i = 0; i < lines; ++i) {
	      const std::string& line = sequences[i];
	      const std::string& line_rev = sequences_rev[i];

	      for (int dir=0; dir < maxdir; ++dir) {
		int first = dir == 0 ? w1 : w2;  // in the second strand the binding sites are in different order, first gap pos 
		int last = first+d;              // last gap pos
		for (int j1=0; j1 < my_cob_params[r].dimer_m[i][d]; ++j1) {  // iterates through start positions
		  double z = my_cob_params[r].overlapping_dimer_Z[i][o][d][dir][j1];
		  get_new_gap_weights(j1, dirs[dir], z, d, fixed_w[tf1], fixed_w[tf2], 
				      my_cob_params[r].oriented_dimer_seeds[o].get<0>(), 
				      my_cob_params[r].oriented_dimer_seeds[o].get<1>(),
				      line, line_rev, 
				      m1,
				      use_multinomial);
		  for (int pos=first; pos < last; ++pos)
		    temp_signal[to_int(line[j1+pos])] += z;

		} // for j1
	      } // for dir
	    } // for i in lines
	    gap_signal_sum += temp_signal;
	    gap_weights[r][o][d] = m1;

	  } // for d, spaced dimer PWMs
	} // for o, spaced dimer PWMs

      } // for end r


      
      /////////////////////////
      //
      // Re-estimate background
      //
      /////////////////////////
      
      std::vector<double> total_signal_sum(4, 0.0);

	
      total_signal_sum += fixed_signal_sum;
      if (local_debug)
	printf("fixed signal: %s\n", print_vector(fixed_signal_sum).c_str());

      // These weren't used to learning fixed PWMs.
      // Only to form background model by subtracting signal from full data
      // if (use_two_strands)
      // 	fixed2_signal_sum /= 2.0;
      total_signal_sum += fixed2_signal_sum;
      if (local_debug)
	printf("fixed2 signal: %s\n", print_vector(fixed2_signal_sum).c_str());

      // if (use_two_strands)
      // 	overlapping_signal_sum /= 2.0;

      total_signal_sum += overlapping_signal_sum;
      if (local_debug)
	printf("overlapping signal: %s\n", print_vector(overlapping_signal_sum).c_str());

      // if (use_two_strands)
      // 	gap_signal_sum /= 2.0;
      total_signal_sum += gap_signal_sum;
      if (local_debug)
	printf("gap signal: %s\n", print_vector(gap_signal_sum).c_str());
      
      if (local_debug)
	printf("Total signal is %s\n", print_vector(total_signal_sum).c_str());
      dvector bg_temp = background_frequencies;    // always counted only in single direction!!!!!
      dmatrix bg_markov_temp = background_frequency_matrix;
      if (local_debug)
	printf("Data distribution is %s\n", print_vector(bg_temp).c_str());

      //recompute the background probabilities, without motif occurences
      for (int i=0; i < 4; ++i) {
	bg_model[i] = bg_temp[i] - total_signal_sum[i];
	assert(bg_model[i] >= 0);
      }

      bg_model_markov = bg_markov_temp - dinucleotide_signal;
      if (local_debug) {
	printf("Unnormalized background distribution (intermed): %s\n", print_vector(bg_model).c_str());
	write_matrix(stdout, bg_model_markov, "Unnormalized dinucleotide background:\n", "%.6f");
      }
      if (use_pseudo_counts) {
	pseudo_counts.add(bg_model);
	pseudo_counts.add(bg_model_markov);
      }
      normalize_vector(bg_model);
      normalize_matrix_rows(bg_model_markov);

      if (local_debug) {
	printf("Background distribution (intermed): %s\n", print_vector(bg_model).c_str());
	write_matrix(stdout, bg_model_markov, "Dinucleotide background:\n", "%.6f");
      }


      
      ///////////////////////////////////////////////////////////////////////////
      //
      // Normalize PWMs, compute ICs for PWM positions and for flanking positions
      //
      ///////////////////////////////////////////////////////////////////////////

      printf("\n");

      // Monomer models
      for (int k=0; k < fixed_p; ++k) {
        if (local_debug)
	  write_matrix(stdout, new_fixed_PWM[k], to_string("Unnormalized fixed matrix %i:\n", k).c_str(), "%.6f");

	if (use_pseudo_counts)
	  pseudo_counts.add(new_fixed_PWM[k]);
	normalize_matrix_columns(new_fixed_PWM[k]);

	assert(is_column_stochastic_matrix(new_fixed_PWM[k]));

	if (local_debug) {
	  write_matrix(stdout, new_fixed_PWM[k], to_string("Intermediate fixed matrix %i:\n", k).c_str(), "%.6f");
	  std::vector<double> ic(fixed_w[k]);
	  for (int i=0; i < fixed_w[k]; ++i)
	    ic[i] = information_content(new_fixed_PWM[k].column(i), bg_model);
	  printf("Information content by columns\n");
	  printf("%s\n", print_vector(ic, "\t", 2).c_str());
	}

	
	if (allow_extension) {
	  if (use_pseudo_counts) {
	    pseudo_counts.add(pred_flank[k]);
	    pseudo_counts.add(succ_flank[k]);
	  }
	  normalize_vector(pred_flank[k]);
	  normalize_vector(succ_flank[k]);
	}
	
      }

      // Overlapping dimer models
      for (int r = 0; r < number_of_cobs; ++r) {
	for (int o=0; o < my_cob_params[r].number_of_orientations; ++o) {
	  for (int d=my_cob_params[r].dmin; d < std::min(0, my_cob_params[r].dmax+1); ++d) {
	    if (use_pseudo_counts)
	      pseudo_counts.add(overlapping_dimer_weights[r][o][d]);
	    normalize_matrix_columns(overlapping_dimer_weights[r][o][d]);
	    assert(is_column_stochastic_matrix(overlapping_dimer_weights[r][o][d]));
	  }
	}
      } // end for r

      // gap models
      for (int r = 0; r < number_of_cobs; ++r) {
	int max_dist_for_deviation = my_cob_params[r].max_dist_for_deviation;
	for (int o=0; o < my_cob_params[r].number_of_orientations; ++o) {
	  for (int d=0; d <= max_dist_for_deviation; ++d) {
	    double mysum = sum(gap_weights[r][o][d]);
	    if (extra_debug) {
	      write_matrix(stdout, gap_weights[r][o][d],
			   to_string("Observed unnormalized gap matrix %s %s %i:\n",
				     my_cob_params[r].name().c_str(), orients[o], d).c_str(), "%.6f");
	      printf("Sum is %e\n", mysum);
	    }
	    if (mysum == 0.0)
	      empty_gap[r][o][d]=true;
	    if (use_pseudo_counts)
	      pseudo_counts.add(gap_weights[r][o][d]);
	    normalize_matrix_columns(gap_weights[r][o][d]);
	    assert(is_column_stochastic_matrix(gap_weights[r][o][d]));
	  }
	}
      } // end for r


      
      for (int k=0; k < fixed_p; ++k) {
	if (keep_monomer_fixed[k]) {
	  new_fixed_PWM[k] = fixed_PWM[k];
	}	  
      }

      
      ///////////////
      //
      // Adjust seeds
      //
      ///////////////
      
      for (int k=0; k < fixed_p; ++k) {
	fixed_seed[k] = string_giving_max_score(new_fixed_PWM[k]);
      }
      
      if (use_multinomial and adjust_seeds) {

	for (int r = 0; r < number_of_cobs; ++r) {
	  int number_of_orientations = my_cob_params[r].number_of_orientations;
	  int dmin = my_cob_params[r].dmin;
	  int dmax = my_cob_params[r].dmax;
	  int tf1 = my_cob_params[r].tf1;
	  int tf2 = my_cob_params[r].tf2;
	  for (int o=0; o < number_of_orientations; ++o) {
	    std::string seed1, seed2;
	    boost::tie(seed1, seed2) = //my_cob_params[r].oriented_dimer_seeds[0]; 
	      get_seeds_according_to_hetero_orientation(o, fixed_seed[tf1], fixed_seed[tf2]);
	    for (int d=dmin; d < std::min(0, dmax+1); ++d) {
	      my_cob_params[r].dimer_seeds[o][d] = create_overlapping_seed(seed1, seed2, d, my_gapped_kmer_context);
	    }  // end for d

	  }  // end for o
	} // end for r

      } // end if use_multinomial


      for (int r=0; r < my_cob_params.size(); ++r) {
	my_cob_params[r].update_oriented_matrices(new_fixed_PWM, fixed_seed);
	my_cob_params[r].compute_expected_matrices(new_fixed_PWM);
      }
      
      //////////////////////
      //
      // Prune dimeric cases
      //
      //////////////////////
      
      for (int r = 0; r < number_of_cobs; ++r) {
	const int& w1 = my_cob_params[r].k1;
	//	const int& w2 = my_cob_params[r].k2;
	for (int o=0; o < my_cob_params[r].number_of_orientations; ++o) {
	  for (int d=my_cob_params[r].dmin; d < std::min(0, my_cob_params[r].dmax+1); ++d) {
	    if (my_cob_params[r].dimer_lambdas[o][d] != 0.0) {
	      double ic = average_information_content(overlapping_dimer_weights[r][o][d].cut(0, w1+d-1, 4, 2-d));
	      double lambda = my_cob_params[r].dimer_lambdas[o][d];
	      bool too_weak = ic < ic_threshold;
	      // bool too_faraway = hamming_distance(my_cob_params[r].dimer_seeds[o][d], 
	      // 					  string_giving_max_score(overlapping_dimer_weights[r][o][d])) >= hamming_threshold;
	      bool too_small = lambda < cob_cutoff;
	      //	      if (too_weak or (use_multinomial and too_faraway) or too_small) {
	      if (too_weak or too_small) {
		my_cob_params[r].dimer_lambdas[o][d] = 0.0;
		if (local_debug) {
		  printf("Excluded dimer case %s %s %i:", my_cob_params[r].name().c_str(), orients[o], d);
		  if (too_weak)
		    printf(" Too weak (ic=%f)", ic);
		  if (too_small)
		    printf(" Too small (lambda=%f)", lambda);
		  printf("\n");
		  //			 too_faraway ? "Too far away" : "");
		}
	      }
	    }
	  }  // end for d
	 for (int d=0; d <= my_cob_params[r].dmax; ++d) {
	   if (my_cob_params[r].dimer_lambdas[o][d] != 0.0) {
	     double lambda = my_cob_params[r].dimer_lambdas[o][d];
	     bool too_small = lambda < cob_cutoff;
	     if (too_small) {
		my_cob_params[r].dimer_lambdas[o][d] = 0.0;
		if (local_debug) {
		  printf("Excluded dimer case %s %s %i: Too small (%f)\n", my_cob_params[r].name().c_str(), orients[o], d, lambda);
		}
	     }
	   }
	 } // end for d
	 
	}  // end for o
      }  // end for r

      
      /////////////////////////////////
      //
      // print overlapping dimer models
      //
      /////////////////////////////////
      
      if (extra_debug) {
	for (int r = 0; r < number_of_cobs; ++r) {
	  for (int o=0; o < my_cob_params[r].number_of_orientations; ++o) {
	    for (int d=my_cob_params[r].dmin; d < std::min(0, my_cob_params[r].dmax+1); ++d) {
	      if (my_cob_params[r].dimer_lambdas[o][d] != 0.0) {
		std::vector<double> ic(my_cob_params[r].dimer_w[d]);
		for (int i=0; i < my_cob_params[r].dimer_w[d]; ++i) {
		  ic[i] = information_content(overlapping_dimer_weights[r][o][d].column(i), bg_model);
		}
		printf("\n");
		write_matrix(stdout, overlapping_dimer_weights[r][o][d],
			     to_string("Observed overlapping dimer case matrix %s %s %i:\n", 
				       my_cob_params[r].name().c_str(), orients[o], d).c_str(), "%.6f");
		printf("Information content of overlapping dimer by columns\n");
		printf("%s\n", print_vector(ic, "\t", 2).c_str());
	      }
	    }
	  }
	} // end for r
      }
      
      /////////////////////////////////
      //
      // print gap dimer models
      //
      /////////////////////////////////

      if (extra_debug) {
	for (int r = 0; r < number_of_cobs; ++r) {
	  int max_dist_for_deviation = my_cob_params[r].max_dist_for_deviation;
	  for (int o=0; o < my_cob_params[r].number_of_orientations; ++o) {
	    for (int d=0; d <= max_dist_for_deviation; ++d) {
	      if (my_cob_params[r].dimer_lambdas[o][d] != 0.0) {
		std::vector<double> ic(my_cob_params[r].dimer_w[d]);
		for (int i=0; i < my_cob_params[r].dimer_w[d]; ++i) {
		  ic[i] = information_content(gap_weights[r][o][d].column(i), bg_model);
		}
		printf("\n");
		write_matrix(stdout, gap_weights[r][o][d],
			     to_string("Observed gap matrix %s %s %i:\n", 
				       my_cob_params[r].name().c_str(), orients[o], d).c_str(), "%.6f");
		printf("Information content of gap by columns\n");
		printf("%s\n", print_vector(ic, "\t", 2).c_str());
	      }
	    }
	  }
	} // end for r
      }

      
      ////////////////////////////
      //
      // Reestimate the deviations
      //
      ////////////////////////////
      
      for (int r = 0; r < number_of_cobs; ++r) {
	const int& w1 = my_cob_params[r].k1;
	//	const int& w2 = my_cob_params[r].k2;
	for (int o=0; o < my_cob_params[r].number_of_orientations; ++o) {
	  for (int d=my_cob_params[r].dmin; d < std::min(0, my_cob_params[r].dmax+1); ++d) {
	    if (my_cob_params[r].dimer_lambdas[o][d] != 0.0) {
	      for (int row = 0; row < 4; ++row) {
		for (int column = w1 + d - 1; column < w1 + 1; ++column) // The overlapping part with flanks of 1bp on both sides
		  new_deviation[r][o][d](row,column) = 
		    overlapping_dimer_weights[r][o][d](row,column) - my_cob_params[r].expected_overlapping_dimer_PWMs[o][d](row,column);
	      }
	    }
	  }  // end for d
	 for (int d=0; d <= my_cob_params[r].max_dist_for_deviation; ++d) {
	   if (my_cob_params[r].dimer_lambdas[o][d] != 0.0) {
	     if (empty_gap[r][o][d])
	       continue;   // The deviation array stays empty
	     int first, last;
	     // The interval [first,last] is the gap.
	     first = w1;
	     last = w1 + d - 1;
	     for (int row = 0; row < 4; ++row) {
	       for (int column = first - 1; column <= last + 1; ++column) // The gap part with flanks of 1bp on both sides
		 new_deviation[r][o][d](row,column) = 
		   gap_weights[r][o][d](row,column) - my_cob_params[r].expected_overlapping_dimer_PWMs[o][d](row,column);
	     } 
	   }
	 } // end for d
	}  // end for o
      }  // end for r

      
      int number_of_parameters = get_number_of_parameters(my_cob_params, fixed_PWM);
      if (local_debug)
	printf("Total number of parameters is %i\n", number_of_parameters);

      

      // Print reestimated mixing parameters
      if (local_debug ) {
	printf("\n");
	printf("Intermediate background lambda is %f\n", background_lambda);
	printf("Intermediate fixed lambdas are %s\n", print_vector(fixed_lambda).c_str());
	if (use_dimers) {
	  for (int r = 0; r < number_of_cobs; ++r) {
	    print_cob(stdout, my_cob_params[r].dimer_lambdas, to_string("Intermediate dimer lambdas %s:\n", my_cob_params[r].name().c_str()), "%e");
	    printf("Intermediate sum of dimer lambdas %s is %.4f\n", my_cob_params[r].name().c_str(),
		   sum(my_cob_params[r].dimer_lambdas));
	  }
	}
      }




      assert(background_lambda >= 0);
      assert(background_lambda <= 1);

      if (allow_extension) {
	for (int k=0; k < fixed_p; ++k) {
	  pred_ic[k] = information_content(pred_flank[k], bg_model);
	  succ_ic[k] = information_content(succ_flank[k], bg_model);
	  printf("Preceeding flank for matrix %i was %s with ic %.2f\n", k,
		 print_vector(pred_flank[k]).c_str(), 
		 pred_ic[k]);
	  printf("Succeeding flank for matrix %i was %s with ic %.2f\n", k,
		 print_vector(succ_flank[k]).c_str(), 
		 succ_ic[k]);
	}
      }
      
     // Print deviation tables
      if (extra_debug) {
	for (int r = 0; r < number_of_cobs; ++r) {
	  int max_dist_for_deviation = my_cob_params[r].max_dist_for_deviation;
	  for (int o=0; o < my_cob_params[r].number_of_orientations; ++o) {
	    for (int d=my_cob_params[r].dmin; d <= max_dist_for_deviation; ++d) {
	      if (my_cob_params[r].dimer_lambdas[o][d] != 0.0) {
		write_matrix(stdout, new_deviation[r][o][d], 
			     to_string("Deviation matrix %s %s %i:\n", my_cob_params[r].name().c_str(),
				       orients[o], d).c_str(), "%.6f");
	      }
	    }
	  }
	}  // end for r
      }


      ///////////////////////
      // Compute distances between old and new models

      typedef BinOp<double>::type my_func;
      myaccumulate<double> acc(0, static_cast<my_func>(std::max));

      printf("\n");


      for (int k=0; k < fixed_p; ++k)
	fixed_dist[k]=distance(new_fixed_PWM[k], fixed_PWM[k]);
      if (local_debug)
	printf("Fixed distances are %s\n", print_vector(fixed_dist).c_str());
      acc(max_element(fixed_dist));

      // Compute distances between 'bridging' table with this round and the previous round
      for (int r = 0; r < number_of_cobs; ++r) {      
	int max_dist_for_deviation = my_cob_params[r].max_dist_for_deviation;
	const int& w1 = my_cob_params[r].k1;
	for (int o=0; o < my_cob_params[r].number_of_orientations; ++o) {
	  for (int d=my_cob_params[r].dmin; d < std::min(0, my_cob_params[r].dmax+1); ++d) {
	    int first;
	    //int last;
	    // The interval [first,last] is either the overlap area or the gap.
	    if (d < 0) {
	      first = w1 + d;
	      //last = w1 - 1;
	    } else {
	      first = w1;
	      //last = w1 + d - 1;
	    }
	    if (my_cob_params[r].dimer_lambdas[o][d] != 0.0)
	      //	      deviation_dist[r][o][d]=distance(new_deviation[r][o][d], my_cob_params[r].deviation[o][d]);
	      deviation_dist[r][o][d]=distance(overlapping_dimer_weights[r][o][d].cut(0, first-1, 4, 2-d),
					       my_cob_params[r].overlapping_dimer_PWM[o][d].cut(0, first-1, 4, 2-d));
	    else
	      deviation_dist[r][o][d]=0.0;
	  } // for d
	  for (int d=0; d <= max_dist_for_deviation; ++d) {
	    int first;
	    //int last;
	    // The interval [first,last] is either the overlap area or the gap.
	    if (d < 0) {
	      first = w1 + d;
	      //last = w1 - 1;
	    } else {
	      first = w1;
	      //last = w1 + d - 1;
	    }
	    if (my_cob_params[r].dimer_lambdas[o][d] != 0.0) {
	      deviation_dist[r][o][d]=distance(gap_weights[r][o][d].cut(0, first-1, 4, 2+d),
					       my_cob_params[r].overlapping_dimer_PWM[o][d].cut(0, first-1, 4, 2+d));
	      overlapping_dimer_weights[r][o][d] = my_cob_params[r].expected_overlapping_dimer_PWMs[o][d];

	      // The next line is doing a little too much, since the flanks of overlapping_dimer_PWM aren't used anywhere
	      overlapping_dimer_weights[r][o][d].inject(gap_weights[r][o][d].cut(0, first-1, 4, 2+d), 0, first-1);
	    }
	    else
	      deviation_dist[r][o][d]=0.0;
	  } // for d
	}
	if (local_debug and my_cob_params[r].dmin < 0) {
	  print_cob(stdout, deviation_dist[r], to_string("Deviation %s distances are\n", 
		    my_cob_params[r].name().c_str()), "%f");
	  printf("Max distance in deviation %s is %f\n", my_cob_params[r].name().c_str(),
		 max_element(deviation_dist[r]));
	}
	acc(max_element(deviation_dist[r]));
      }

      double total_maximum_distance = acc.get();
      if (local_debug)
        printf("Total maximum distance is %f\n\n", total_maximum_distance);


      ///////////////////////
      // Replace old models

      
      fixed_PWM = new_fixed_PWM;              // Replace old fixed models
      
      for (int r = 0; r < number_of_cobs; ++r) {
	my_cob_params[r].deviation = new_deviation[r];
	my_cob_params[r].overlapping_dimer_PWM = overlapping_dimer_weights[r];  // This is needed only to compute the distance for convergence
      }

      if (local_debug) {
	for (int k=0; k < fixed_p; ++k)
	  fixed_av_ic[k] = average_information_content(fixed_PWM[k], bg_model);
	printf("Intermediate average information content of fixed models: %s\n", print_vector(fixed_av_ic).c_str());
      }

      bg_model_rev = std::vector<double>(bg_model.rbegin(), bg_model.rend()); // use this model when considering the reverse strand
      bg_model_markov_rev = reverse_complement_markov_model(bg_model_markov);

      mll = complete_data_log_likelihood(fixed_PWM,
					 bg_model, bg_model_markov,
					 bg_model_rev, bg_model_markov_rev,
					 fixed_lambda, background_lambda, my_cob_params, fixed_Z, sequences, sequences_rev, fixed_w, fixed_m);
      if (local_debug)
	printf("Log likelihood is %f\n", mll);

      bool deviation_converged = true;
      for (int r=0; r < my_cob_params.size(); ++r) {
	if (max_element(deviation_dist[r]) >= epsilon) {
	  deviation_converged = false;
	  break;
	}
      }

      convergence_criterion_reached = 
	max_element(fixed_dist) < epsilon and deviation_converged;

      if ((convergence_criterion_reached or round+1 >= max_iter) and get_full_flanks)
	boost::tie(flank_fixed_PWM, flank_dimer_PWM) =
	  get_models_with_flanks(sequences,
				 sequences_rev,
				 fixed_seed,
				 fixed_w,
				 Lmax,
				 fixed_m,
				 my_cob_params,
				 fixed_Z);

      if (local_debug) {
 	printf("---------------------------------------------------\n");
        fflush(stdout);
      }
    } // for round
    iterations.push_back(round);

    if (unbound.length() != 0)
      print_background_sequences(sequences, fixed_Z, my_cob_params, fixed_PWM, bg_model, bg_model_markov);
    
    ++extension_round;
      

    if (not allow_extension or extension_round == max_extension_round)
      break;

    bool all_below = true; 
    for (int k=0; k < fixed_p; ++k) {
      if (pred_ic[k] >= extension_ic_threshold or succ_ic[k] >= extension_ic_threshold) {
	all_below = false;
	break;
      }
    }
    if (all_below)
      break;
    for (int k=0; k < fixed_p; ++k) {
      char nucs[] = "ACGT";
      if (pred_ic[k] >= extension_ic_threshold) {
	fixed_PWM[k].insert_column(0, pred_flank[k]);
	fixed_seed[k].insert(0, 1, nucs[arg_max(pred_flank[k])]);  // THIS SHOULD BE CORRECTED. HOW IS NEW BASE DECIDED, NOW MAX
      }
      if (succ_ic[k] >= extension_ic_threshold) {
	fixed_PWM[k].insert_column(fixed_PWM[k].get_columns(), succ_flank[k]);
	fixed_seed[k].insert(fixed_seed[k].length(), 1, nucs[arg_max(succ_flank[k])]);  // THIS SHOULD BE CORRECTED
      }
      if (pred_ic[k] >= extension_ic_threshold or succ_ic[k] >= extension_ic_threshold) {
	write_matrix(stdout, fixed_PWM[k], to_string("Extended matrix %i:\n", k).c_str(), "%.6f");
      }
    }    
  } // while extension_round  

  /******************************************************************************************/

  printf("%s\n", std::string(40, '*').c_str());
  printf("Results:\n");

  
  // Print deviations
  for (int r = 0; r < number_of_cobs; ++r) {
    int max_dist_for_deviation = my_cob_params[r].max_dist_for_deviation;
    for (int o=0; o < my_cob_params[r].number_of_orientations; ++o) {
      for (int d=my_cob_params[r].dmin; d <= max_dist_for_deviation; ++d) {
	if (my_cob_params[r].dimer_lambdas[o][d] != 0.0) {
	  write_matrix(stdout, my_cob_params[r].deviation[o][d], 
		       to_string("Deviation matrix %s %s %i:\n", my_cob_params[r].name().c_str(),
				 orients[o], d).c_str(), "%.6f");
	}
      }
    }
  }  // end for r

  // Print flanks

  if (get_full_flanks) {
    for (int k=0; k < fixed_p; ++k) {
      write_matrix(stdout, flank_fixed_PWM[k], to_string("Flank fixed matrix %i:\n", k), "%.6f");
    }
  
    for (int r = 0; r < number_of_cobs; ++r) {
      int dmin = my_cob_params[r].dmin;
      int dmax = my_cob_params[r].dmax;
      int no = my_cob_params[r].number_of_orientations;
      int tf1 = my_cob_params[r].tf1;
      int tf2 = my_cob_params[r].tf2;
      for (int o=0; o < no; ++o) {
	for (int d=dmin; d <= dmax; ++d) {
	  if (my_cob_params[r].dimer_lambdas[o][d] == 0.0)
	    continue;
	  write_matrix(stdout, flank_dimer_PWM[r][o][d], to_string("Flank dimer case matrix %i-%i %s %i:\n", tf1,tf2, orients[o], d), "%.6f");
	}
      }
    
    }
  }
  
  printf("Maximum log likelihood in the beginning: %i\n", 
	 (int)first_maximum_log_likelihood);
  printf("Maximum log likelihood in the end: %i\n", (int)mll);

  //printf("AIC score: %lg\n", aic_score(maximum_log_likelihood, lines, w[0]));
  printf("EM-algorithm took %s = %d iterations \n", print_vector(iterations, "+", 0).c_str(), sum(iterations));
  printf("\n");
  printf("Background lambda is %f\n", background_lambda);
  printf("Fixed lambdas are %s\n", print_vector(fixed_lambda).c_str());
  if (use_output) {
    std::string file = to_string("%s/monomer_weights.txt", outputdir.c_str());
    FILE* fp = fopen(file.c_str(), "w");
    if (fp == NULL) {
      fprintf(stderr, "Couldn't open file %s for writing\n", file.c_str());
      exit(1);
    }
    fprintf(fp, "%s\n", print_vector(names, ",", 0).c_str());
    fprintf(fp, "%s\n", print_vector(fixed_lambda, ",", 6).c_str());
    fclose(fp);
  }


  double total_dimer_lambda=0.0;
  if (use_dimers) {
    for (int r = 0; r < number_of_cobs; ++r) {
      int max_dist_for_deviation = my_cob_params[r].max_dist_for_deviation;
      const cob_params_t& cp = my_cob_params[r];
      print_cob(stdout, cp.dimer_lambdas, 
		to_string("Dimer lambdas %s:\n", cp.name().c_str()), 
		"%.5f");
      if (use_output)
	write_cob_file(to_string("%s/%s-%s.cob", outputdir.c_str(), names[cp.tf1].c_str(), names[cp.tf2].c_str()), cp.dimer_lambdas);
      printf("Sum of dimer lambdas of cob table %s is %.4f\n", 
	     cp.name().c_str(),
	     sum(cp.dimer_lambdas));
      total_dimer_lambda += sum(cp.dimer_lambdas);
      std::vector<std::string> excluded;
      MY_FOREACH(o, cp.dimer_lambdas) {
	MY_FOREACH(d, cp.dimer_lambdas[o]) {
	  if (d <= max_dist_for_deviation) {
	    if (cp.dimer_lambdas[o][d] == 0.0) {
	      excluded.push_back(to_string("%s %i", orients[o], d));
	    } else {
	      if (use_output) {
		write_matrix_file(to_string("%s/%s-%s.%s.%i.dev", outputdir.c_str(), names[cp.tf1].c_str(), names[cp.tf2].c_str(),
					    orients[o], d),
				  cp.deviation[o][d]);
	      }
	    }
	  }
	  else { // d>=0
 	    if (cp.dimer_lambdas[o][d] == 0.0) {
              excluded.push_back(to_string("%s %i", orients[o], d));
	    }
	  }
	}
      }
      printf("Cob %s excluded cases: %s\n", cp.name().c_str(), print_vector(excluded).c_str());
    }
  }

  assert(fabs(background_lambda + sum(fixed_lambda) + total_dimer_lambda) - 1.0 < 0.001);

  printf("\n");

  printf("Background distribution: %s\n", print_vector(bg_model).c_str());
  printf("Fixed seeds are %s\n", print_vector(fixed_seed).c_str());

  for (int k=0; k < fixed_p; ++k) {
    write_matrix(stdout, fixed_PWM[k], to_string("Fixed matrix %i:\n", k), "%.6f");
    if (use_output)
      write_matrix_file(to_string("%s/%s.pfm", outputdir.c_str(), names[k].c_str()), fixed_PWM[k]);
    fixed_av_ic[k] = average_information_content(fixed_PWM[k], bg_model);
  }
  printf("Average information content in fixed PWMs: %s\n", print_vector(fixed_av_ic).c_str());

  printf("%s\n", std::string(40, '*').c_str());

  
  return fixed_PWM;

} // end multi_profile_em_algorithm





bool
ends_with(const std::string& s, const std::string& ending)
{
  int k=ending.length();
  if (s.length() < k)
    return false;

  return s.substr(s.length()-k, k) == ending;
}

std::string
replace_ending(std::string s, const std::string& old_ending, const std::string& new_ending)
{
  int k = old_ending.length();
  assert(s.length() >= k);

  s.replace(s.length()- k, k, new_ending);
  return s;
}


// print a help message of positional parameters,
// using boost's program options library
void
print_positional_options(const po::positional_options_description& popt,
			 const po::options_description& o)
{
  printf("Positional parameters:\n");
  int count = popt.max_total_count();
  for (int i=0; i < count; ++i) {
    const std::string& name = popt.name_for_position(i);
    const std::string& desc =  o.find(name, false).description();
    printf("  %-22s%s\n", name.c_str(), desc.c_str());
  }
}

// Not used currently
class mytempfile
{
public:

  mytempfile() 
  {
    char templ[] = "/tmp/mytempXXXXXX";
    fd = mkstemp(templ);
    //    f = fdopen(fd);
    name = templ;
    reference_count = new int(1);
    //    *reference_count = 1;
  }

  mytempfile(const mytempfile& orig) {
    fd = orig.fd;
    reference_count = orig.reference_count;
    name = orig.name;
    ++*reference_count;
  }

  ~mytempfile() 
  {
    --*reference_count;
    if (*reference_count == 0) {
      //    fclose(f);
      close(fd);
      remove(name.c_str());
      delete reference_count;
    }
  }

  std::string name;
  //  FILE* f;
private:


  // NOT IMPLEMENTED
  mytempfile&
  operator==(const mytempfile& rhs)
  {
    return *this;
  }

  int fd;
  int* reference_count;
};


// Not used!!!
void
pwm_list_to_logos(const std::string& result_filename, const std::vector<dmatrix>& matrices)
{
  int p = matrices.size();
  std::string append_command = "convert ";
  std::vector<mytempfile> output;
  bool use_spacek = true;  // Use Jussi's program to create logos

  std::string logo_program;
  if (use_spacek)
    logo_program = "spacek40 --logo";
  else
    logo_program = "to_logo.sh";

  for (int i=0; i < p; ++i) {
    mytempfile input;
    output.push_back(mytempfile());

    write_matrix_file(input.name, matrices[i]);

    std::string logo_command = to_string("%s %s %s > /dev/null", 
					 logo_program.c_str(), input.name.c_str(), output[i].name.c_str());
    printf("%s\n", logo_command.c_str());

    int code = system(logo_command.c_str());
    //error(code != 0, "system returned non-zero result");
    if (code != 0)
      return;
    append_command += " ";
    append_command += output[i].name;
    if (use_spacek)
      append_command += ".png";
  }
  append_command += " ";
  if (use_spacek)
    append_command += " -append ";
  else
    append_command += " +append ";
  append_command += result_filename;
  printf("%s\n", append_command.c_str());
  int code = system(append_command.c_str());
  //error(code != 0, "system returned non-zero result");
  if (code != 0)
    return;

}


dmatrix
multinomial1_multimer_bernoulli_corrected(const std::string& seed, const std::vector<std::string>& sequences)
{
  int lines = sequences.size();
  int L = sequences[0].length();
  int k = seed.length();
  int seed_count;
  int multinomial_count;
  dmatrix multinomial1_motif;
  boost::tie(multinomial1_motif, seed_count, multinomial_count) = find_snips_multimer_helper(seed, sequences);
  double mean = lines * (use_two_strands ? 2 : 1) * (L-k+1) * pow(4, -k);
      
  dmatrix mean_matrix(4, k);
  mean_matrix.fill_with(mean);
  dmatrix result;
  result = multinomial1_motif - mean_matrix;
  
  result.apply(cut);

  return result;
}

cob_params_t
create_cob(cob_combination_t cob_combination,
	   const std::vector<std::string>& fixed_seeds,
	   const std::vector<dmatrix>& fixed_M,
	   double dimer_lambda_fraction,
	   const overlapping_dimer_cases_t& overlapping_dimer_cases, 
	   const spaced_dimer_cases_t& spaced_dimer_cases,
	   const std::vector<std::string>& sequences,
	   const std::vector<int>& L,
	   const gapped_kmer_context& my_gapped_kmer_context, int dmin, int dmax, int max_dist_for_deviation)
{
  int tf1, tf2;
  boost::tie(tf1, tf2) = cob_combination;
  //  int L = sequences[0].length();
  int number_of_orientations = tf1 == tf2 ? 3 : 4;
  
  const int number_of_spaced_dimer_cases = use_dimers ? (dmax + 1) * number_of_orientations : 0;
  const int number_of_overlapping_dimer_cases = use_dimers ? (-dmin) * number_of_orientations : 0;
    

  //typedef typename double_cob_type::extent_range range;
  typedef boost::multi_array<double, 2>::extent_range range;
  boost::multi_array<double, 2> dimer_lambdas(boost::extents[number_of_orientations][range(dmin,dmax+1)]);
  int temp = std::min(0, dmax+1);
  boost::multi_array<std::string, 2> dimer_seeds(boost::extents[number_of_orientations][range(dmin,temp)]);

  boost::multi_array<dmatrix, 2> overlapping_dimer_PWM(boost::extents[number_of_orientations][range(dmin, max_dist_for_deviation+1)]);


  if (use_dimers) {
    // These were given on the command line
    double uniform_lambda = dimer_lambda_fraction / (number_of_overlapping_dimer_cases + number_of_spaced_dimer_cases);
    /*
    for (int i=0; i < overlapping_dimer_cases.size(); ++i) {
      //overlapping_dimer_lambda.push_back(uniform_lambda);
      int o = overlapping_dimer_cases[i].get<2>();
      int d = overlapping_dimer_cases[i].get<3>();
      std::string seed = overlapping_dimer_cases[i].get<4>();
      dimer_seeds[o][d] = seed;
      dimer_lambdas[o][d] = uniform_lambda;
    }
     
    // These were given on the command line
    for (int i=0; i < spaced_dimer_cases.size(); ++i) {
      int o = spaced_dimer_cases[i].get<2>();
      int d = spaced_dimer_cases[i].get<3>();
      dimer_lambdas[o][d] = uniform_lambda;
      //      spaced_dimer_lambda.push_back(uniform_lambda);
    }
    */
    
    // Fill dimer matrix with the rest of overlapping cases
    for (int o=0; o < number_of_orientations; ++o) {
      std::string seed1, seed2;
      boost::tie(seed1, seed2) = get_seeds_according_to_hetero_orientation(o, fixed_seeds[tf1], fixed_seeds[tf2]);
      for (int d=dmin; d < std::min(0, dmax+1); ++d) {
	if (dimer_seeds[o][d] != "")                  // value already set
	  continue;
	dimer_seeds[o][d] = create_overlapping_seed(seed1, seed2, d, my_gapped_kmer_context);
	dimer_lambdas[o][d] = uniform_lambda;
      }
    }
 
    // Fill dimer matrix with the rest of spaced cases
    for (int o=0; o < number_of_orientations; ++o) {
      for (int d=0; d <= dmax; ++d) {
	if (dimer_lambdas[o][d] != 0.0)                  // value already set
	  continue;
	dimer_lambdas[o][d] = uniform_lambda;
	//	spaced_dimer_lambda.push_back(uniform_lambda);
      }
    }

  } // end if use_dimers


  //  if (seeds_given) {
  std::vector<boost::tuple<dmatrix,dmatrix> > oriented_dimer_matrices(4);
  for (int o=0; o < number_of_orientations; ++o) {
    oriented_dimer_matrices[o] = get_matrices_according_to_hetero_orientation(o, fixed_M[tf1], fixed_M[tf2]);
    
    for (int d=dmin; d < std::min(0, dmax+1); ++d) {
      std::string seed = dimer_seeds[o][d];
      overlapping_dimer_PWM[o][d] = multinomial1_multimer_bernoulli_corrected(seed, sequences);
      //       = find_snips_multimer(seed, sequences);
      if (use_pseudo_counts)
	pseudo_counts.add(overlapping_dimer_PWM[o][d]);
      normalize_matrix_columns(overlapping_dimer_PWM[o][d]);
      /*
	write_matrix(stdout, overlapping_dimer_PWM[o][d], 
	to_string("Initial overlapping_dimer matrix %i-%i %s %i from seed %s:\n", 
	tf1, tf2,
	orients[o], d, seed.c_str()), "%.6f");
      */
    }  // end for d
    for (int d=0; d <= max_dist_for_deviation; ++d) {
      //      int dimer_len = fixed_seeds[tf1].length() + fixed_seeds[tf2].length() + d;
      //      overlapping_dimer_PWM[o][d] = dmatrix(4, dimer_len);
      overlapping_dimer_PWM[o][d] =
	normalize_matrix_columns_copy(matrix_product(oriented_dimer_matrices[o].get<0>(),
						     oriented_dimer_matrices[o].get<1>(),
						     d));


	}
  } // end for o
  //  }

  cob_params_t cb(tf1, tf2, dimer_lambdas, dimer_seeds, overlapping_dimer_PWM, fixed_seeds, L, dmin, dmax, max_dist_for_deviation);

  return cb;
}


dmatrix
get_meme_init_pwm(const std::string& seed)
{
  double x = 0.5;
  int k = seed.length();
  dmatrix result(4, k);
  for (int c=0; c < k; ++c) {
    int seed_row = to_int(seed[c]);
    for (int r=0; r < 4; ++r) {
      result(r,c) = (r == seed_row ? x : x / 3);
    }
  }

  return result;
}

std::vector<cob_combination_t>
get_all_cob_combinations(int fixed_p)
{
  std::vector<cob_combination_t> cob_combinations;
  for (int i=0; i < fixed_p; ++i) {
    for (int j=i; j < fixed_p; ++j) {
      cob_combinations.push_back(boost::make_tuple(i,j));
    }
  }
  return cob_combinations;
}

int main(int argc, char* argv[])
{
  TIME_START(t);
  WALL_TIME_START(t2);
#ifdef FE_NOMASK_ENV
  feenableexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW | FE_UNDERFLOW);   // These exceptions cause trap to occur
#endif

  if (argc > 1) {
    printf("Command line was: ");
    print_command_line(argc, argv);
  }
  
  use_two_strands = true;
  bool dmin_given = false;
  bool dmax_given = false;
  bool max_dist_for_deviation_given = false;
  int unique_param = -1; // either -1, 0, 1, ...
  std::vector<int> dmina;
  std::vector<int> dmaxa;
  std::vector<int> max_dist_for_deviationa;
  using namespace boost;
  using std::string;
  using std::cout;
  double epsilon = 0.01;
  double extension_threshold = 2.00;
  std::vector<cob_combination_t> cob_combinations;
  string prior_parameter = "dirichlet";

  std::vector<std::string> sequences;
  std::vector<bool> keep_monomer_fixed;
  
  combine_seeds_func = combine_seeds_methods[default_combine_method].function;

  int number_of_threads=1;
#if defined(_OPENMP)
  {
    char* p = getenv("OMP_NUM_THREADS"); 
    if (p)
      number_of_threads = atoi(p);
  }
#endif

  
  // Declare the supported options.
  po::options_description hidden("Allowed options");
  hidden.add_options()
    ("sequencefile",                           "Name of the sequence file containing sequences")
    ("parameterlist", po::value<string>(), "Seeds, or names of the initial matrix files, separated by commas")
    ;

  po::options_description nonhidden("Allowed options");
  nonhidden.add_options()
    ("help",                           "Produce help message")
    ("matrices",                          "Matrix filenames are given as parameters instead of seeds")
    ("keep-monomer-fixed",            po::value<std::string>(), "A list of indices of monomers to keep fixed during EM algorithms. You can also specify 'all' as parameter, default: no monomers are fixed.")
    ("disable-multinomial",                    m("Use plain alignment of sequences instead of the multinomial model", not use_multinomial).c_str())
    ("directional-seed",                    m("Are overlapping seeds required to be non-palindromic", require_directional_seed).c_str())
    ("hamming-radius", po::value<int>(),     m("Hamming radius", hamming_radius).c_str())
    ("meme-init",                    m("Derive the initial PWM from the seed using MEME style initialization instead of multinomial style", use_meme_init).c_str())
    ("flanks",  m("Get full flanks for each model", get_full_flanks).c_str())

    ("max-iter", po::value<int>(),   m("Maximum number of iterations", max_iter).c_str())
    ("epsilon", po::value<double>(),        m("Epsilon defines convergence of EM (double). Elementwise maximum distance between two consequent matrices", epsilon).c_str())

    ("cob", po::value<string>(),     "List of cob tables wanted. Example format: --cob '0-0,1-1,0-1', "
                                     "where the numbers refer to the indices of the list of monomers. "
                                     "Note, that the indexing starts from zero. "
                                     "Use option '--cob all' to create all combinations.")
    ("cob-cutoff", po::value<double>(),        m("Cob cutoff constant", cob_cutoff).c_str())

    ("dmin", po::value<std::string>(),   "Smallest negative distance in dimer, comma separated list or a single global value, default: half of the length of the shorter monomer")
    ("dmax", po::value<std::string>(),   m("Maximum positive distance in dimer, comma separated list or a single global value", global_dmax).c_str())
    ("max-gap-learned", po::value<std::string>(),   m("Maximum positive distance in dimer for which the gap is learned, comma separated list or a single global value", global_max_dist_for_deviation).c_str())
    ("minimum-distance-for-learning", po::value<int>(),   m("Minimum distance in dimer cases for a site to be used for monomer model learning", minimum_distance_for_learning).c_str())
    ("ic-threshold", po::value<double>(),   m("Information content threshold for an overlapping dimer to be accepted", ic_threshold).c_str())
    ("min-fraction-for-learning", po::value<double>(),   m("The monomers are learned from the dimeric cases alone, when the fraction of dimeric cases "
							   "with distance >=minimum-distance-for-learning is at least this threshold", learning_fraction).c_str())
    
    // positional-background is not implemented currently
    //    ("positional-background",        m("Use positional model for background", use_positional_background).c_str())

    // not implemented
    //("markov-model",                 m("Use markov model for background", use_markov_background).c_str())

    // not used currently ("extension-threshold", po::value<double>(),   m("Information content threshold for extending models. Default: no extension", extension_threshold).c_str())
    //("fixed-lambdalist", po::value<string>(),     "Comma separated list of mixing parameters, one for each matrix")


    ("unbound", po::value<string>(&unbound), "File to store the unbound sequences")
    ("output", m("Write model parameters to files", use_output).c_str())
    ("names", po::value<string>(),     "Names for the monomer models. Comma separated list. Default: TF0,TF1, ...")
    ("outputdir", po::value<string>(), m("Directory where to put the learned matrices", outputdir).c_str())
    ("number-of-threads", po::value<int>(), m("Number of parallel threads", number_of_threads).c_str())
    ("prior", po::value<std::string>(), m("Choose either addone or dirichlet prior, or none for no pseudo-counts", prior_parameter).c_str())
    ("single-strand", m("Only consider binding sites in forward strand", not use_two_strands).c_str())
    ("unique", po::value<std::string>(), "Uniqueness of sequences. Either off, unique, 1, 2, 3, ..., default: off")
    ("quiet", m("Don't print intermediate results", not local_debug).c_str())
    ;

  po::options_description desc("Allowed options");
  desc.add(hidden).add(nonhidden);

  po::positional_options_description popt;
  popt.add("sequencefile", 1);
  popt.add("parameterlist", 1);

  
  string seqsfile;
  //string freqs;
  std::vector<std::string> fixedmatrixfilelist;
  std::vector<std::string> fixedseedlist;

  double background_lambda;
  std::vector<double> fixed_lambda;
  std::vector<double> overlapping_dimer_lambda;
  std::vector<double> spaced_dimer_lambda;


  overlapping_dimer_cases_t overlapping_dimer_cases;
  spaced_dimer_cases_t spaced_dimer_cases;

  po::variables_map vm;
  int fixed_p = 0; // number of fixed models
  // int overlapping_dimer_p = 0; // number of overlapping models
  // int spaced_dimer_p = 0; // number of spaced models
  boost::multi_array<std::string, 2> dimer_seeds;

   ////////////////////////////////
  // parse command line parameters
  ////////////////////////////////

  try {
    po::store(po::command_line_parser(argc, argv).
	      options(desc).positional(popt).run(), vm);
    po::notify(vm);

    // are all positional parameters present
    bool all_positional = vm.count("sequencefile") && vm.count("parameterlist");

    if (vm.count("help") || not all_positional) {
      cout << basename(argv[0]) << " version " << PACKAGE_VERSION << std::endl;
      cout << "Usage:\n" 
	   << basename(argv[0]) << " [options] sequencefile seedlist\n"
	   << basename(argv[0]) << " [options] --matrices sequencefile matrixfilelist\n";
      print_positional_options(popt, hidden);
      cout << nonhidden << "\n";
      return 1;
    }

    if (vm.count("markov-model")) 
      use_markov_background = true;

    if (vm.count("single-strand")) 
      use_two_strands = false;

    if (vm.count("unique")) {
      std::string arg = vm["unique"].as< string >();
      if (arg == "off")
	unique_param = -1;
      else if (arg == "unique")
	unique_param = 0;
      else {
	unique_param = atoi(arg);
	error(unique_param <= 0, "Argument to --unique must be one of the following: off, unique, 1, 2, 3, ...");
      }
    }
    

    if (vm.count("quiet")) 
      local_debug = false;

    if (vm.count("flanks")) 
      get_full_flanks = true;

    if (vm.count("positional-background")) 
      use_positional_background = true;

    if (vm.count("disable-multinomial")) 
      use_multinomial = false;
 
    if (vm.count("directional-seed")) 
      require_directional_seed = true;

    if (vm.count("meme-init")) 
      use_meme_init = true;

    if (vm.count("output")) 
      use_output = true;

    if (vm.count("outputdir")) {
      outputdir = vm["outputdir"].as< string >();
      use_output = true;
    }

    if (vm.count("prior")) {
      prior_parameter = vm["prior"].as< string >();
      if (prior_parameter == "addone" || prior_parameter == "dirichlet")
	use_pseudo_counts=true;
      else if (prior_parameter == "none") {
	use_pseudo_counts=false;
      }
      else {
	error(true, to_string("Invalid prior: %s.\nAllowed options are addone, dirichlet, or none\n", prior_parameter.c_str()));
      }
    }   

#if defined(_OPENMP)
    if (vm.count("number-of-threads")) {
      number_of_threads = vm["number-of-threads"].as< int >();
    }
    omp_set_num_threads(number_of_threads);
#endif

    if (vm.count("hamming-radius"))
      hamming_radius = vm["hamming-radius"].as< int >();
    error(hamming_radius <= 0, "hamming-radius must be positive integer.");

    if (vm.count("max-iter"))
      max_iter = vm["max-iter"].as< int >();
    error(max_iter <= 0, "max-iter must be positive integer.");

    if (vm.count("cob-cutoff"))
      cob_cutoff = vm["cob-cutoff"].as< double >();
    error(cob_cutoff <= 0.0 or cob_cutoff >= 1.0, "cob-cutoff must be between 0 and 1.");

    if (vm.count("minimum-distance-for-learning"))
      minimum_distance_for_learning = vm["minimum-distance-for-learning"].as< int >();
    error(minimum_distance_for_learning < 0, "minimum-distance-for-learning must be non-negative.");

    if (vm.count("min-fraction-for-learning"))
      learning_fraction = vm["min-fraction-for-learning"].as< double >();
    error(learning_fraction < 0.0 or learning_fraction > 1.0, "min-fraction-for-learning must be between 0 and 1.");

    if (vm.count("ic-threshold"))
      ic_threshold = vm["ic-threshold"].as< double >();
    error(ic_threshold < 0.0 or ic_threshold > 1.0, "ic-threshold must be between 0 and 1.");

    seqsfile   = vm["sequencefile"].as< string >();

    if (vm.count("overlap-maximize")) {
      hamming_distance_overlapping_seeds_N = vm["overlap-maximize"].as< int >();
      maximize_overlapping_seeds = true;
    }
    error(hamming_distance_overlapping_seeds_N < 0 or hamming_distance_overlapping_seeds_N > 2, "overlap maximize must be either 0, 1 or 2.");

    // if (vm.count("min-flank"))
    //   min_flank = vm["min-flank"].as< int >();
    // error(min_flank <= 0, "min-flank must be positive.");

    if (vm.count("overlap-combine")) {
      std::string method = vm["overlap-combine"].as<std::string>();
      if (method == combine_seeds_methods[COMBINE_SEEDS_OR].name)
	combine_seeds_func = combine_seeds_methods[COMBINE_SEEDS_OR].function;
      else if (method == combine_seeds_methods[COMBINE_SEEDS_AND].name)
	combine_seeds_func = combine_seeds_methods[COMBINE_SEEDS_AND].function;
      else if (method == combine_seeds_methods[COMBINE_SEEDS_N].name)
	combine_seeds_func = combine_seeds_methods[COMBINE_SEEDS_N].function;
      else {
	error(true, "The overlap-combine method must be either union, intersection or N.");
      }

    }

    // Command line parameters for PWM models

    if (vm.count("parameterlist")) {
      if (vm.count("matrices")) {         // matrices were given
	fixedmatrixfilelist = split(vm["parameterlist"].as< string >(), ',');
	fixed_p = fixedmatrixfilelist.size();   // number of monomer models
	seeds_given=false;
      } else {                         // Initial seed were given
	fixedseedlist = split(vm["parameterlist"].as< string >(), ',');
	BOOST_FOREACH(std::string& seed, fixedseedlist) {
	  boost::to_upper(seed);
	  error(not is_iupac_string(seed), "Seed must be an IUPAC string");
	}
	seeds_given=true;
	fixed_p = fixedseedlist.size();   // number of models
      }
    }

    if (vm.count("cob")) {
      if (vm["cob"].as<std::string>() == "all")
	cob_combinations = get_all_cob_combinations(fixed_p);
      else {
	std::vector<std::string> cob_tables = split(vm["cob"].as<std::string>(), ',');
	BOOST_FOREACH(std::string s, cob_tables) {
	  std::vector<std::string> temp =split(s, '-');
	  error(temp.size() != 2,
		to_string("Cob combinations given to --cob should be pairs of integers separated by '-'. Error in %s", s.c_str()));

	  int tf1 = atoi(temp[0]);
	  int tf2 = atoi(temp[1]);
	  error(tf1 < 0 or tf1 >= fixed_p, to_string("TF index pair in %s out of range.", s.c_str()));
	  error(tf2 < 0 or tf2 >= fixed_p, to_string("TF index pair in %s out of range.", s.c_str()));
	  //	  assert(0 <= tf1 and tf1 < fixed_p);
	  //	  assert(0 <= tf2 and tf2 < fixed_p);
	  cob_combinations.push_back(make_tuple(tf1, tf2));
	}
      }
    }
    int number_of_cobs = cob_combinations.size();
    std::string cob_combinations_string;
    for (int i=0; i < number_of_cobs; ++i) {
      std::string c = to_string("%i-%i", cob_combinations[i].get<0>(), cob_combinations[i].get<1>());
      if (i == 0)
	cob_combinations_string = c;
      else
	cob_combinations_string.append("," + c);
    }
    printf("Cob combinations are %s\n", cob_combinations_string.c_str());
    if (vm.count("dmin")) {
      std::vector<std::string> dmin_array = split(vm["dmin"].as<std::string>(), ',');

      if (dmin_array.size() == 1 and number_of_cobs > 0) {
	int value = atoi(dmin_array[0]);
	error(value > 0, "dmin values must be non-positive");
	dmina.assign(number_of_cobs, value);
      }
      else {
	error(dmin_array.size() != number_of_cobs, "Different number of cob tables and dmin values");
	BOOST_FOREACH(std::string s, dmin_array) {
	  int value = atoi(s);
	  error(value > 0, "dmin values must be non-positive");
	  dmina.push_back(value);
	}
      }
      dmin_given = true;
    }

    if (vm.count("dmax")) {
      std::vector<std::string> dmax_array = split(vm["dmax"].as<std::string>(), ',');

      if (dmax_array.size() == 1 and number_of_cobs > 0) {
	int value = atoi(dmax_array[0]);
	error(value < 0, "dmax values must be non-negative");
	dmaxa.assign(number_of_cobs, value);
      }
      else {
	error(dmax_array.size() != number_of_cobs, "Different number of cob tables and dmax values");
	BOOST_FOREACH(std::string s, dmax_array) {
	  int value = atoi(s);
	  error(value < 0, "dmax values must be non-negative");
	  dmaxa.push_back(value);
	}
      }
      dmax_given = true;
    }

    if (vm.count("max-gap-learned")) {
      std::vector<std::string> max_dist_for_deviation_array = split(vm["max-gap-learned"].as<std::string>(), ',');

      if (max_dist_for_deviation_array.size() == 1 and number_of_cobs > 0) {
	int value = atoi(max_dist_for_deviation_array[0]);
	//error(value < 0, "dmax values must be non-negative");
	max_dist_for_deviationa.assign(number_of_cobs, value);
      }
      else {
	error(max_dist_for_deviation_array.size() != number_of_cobs,
	      "Different number of cob tables and max-gap-learned values");
	BOOST_FOREACH(std::string s, max_dist_for_deviation_array) {
	  int value = atoi(s);
	  //error(value < 0, "dmax values must be non-negative");
	  max_dist_for_deviationa.push_back(value);
	}
      }
      max_dist_for_deviation_given = true;
    }

    if (vm.count("names")) {         // Names given for the monomer models/seeds
      std::vector<std::string> namelist = split(vm["names"].as< string >(), ',');
      error(namelist.size() != fixed_p, "For each seed/pfm a name should be given");
      names = namelist;
      // BOOST_FOREACH(std::string& seed, fixedseedlist) {
      // }
    }
    else {
      names.clear();
      for (int i=0; i < fixed_p; ++i)
	names.push_back(to_string("TF%i", i));
    }

    // Command line parameters for cob models
    /*
    if (vm.count("overlapping-dimer") and use_dimers) {
      std::vector<std::string> v = vm["overlapping-dimer"].as<std::vector<std::string> >();
      printf("Overlapping dimer cases:\n");
      for (int a=0; a < v.size(); ++a) {
	std::vector<std::string> parts = split(v[a], ',');
	assert(parts.size() == 5);

	int i = atoi(parts[0]);
	int j = atoi(parts[1]);
	int o = string_to_orientation(parts[2]);
	int d = atoi(parts[3]);
	std::string seed = parts[4];

	assert(0 <= i && i < fixed_p);
	assert(0 <= j && j < fixed_p);

	overlapping_dimer_cases.push_back(boost::make_tuple(i, j, o, d, seed));
	printf("%i: %s\n", a, v[a].c_str());
      }
      //      overlapping_dimer_p = overlapping_dimer_cases.size();
    }

    if (vm.count("spaced-dimer") and use_dimers) {
      std::vector<std::string> v = vm["spaced-dimer"].as<std::vector<std::string> >();
      printf("Spaced dimer cases:\n");
      for (int a=0; a < v.size(); ++a) {
	std::vector<std::string> parts = split(v[a], ',');
	assert(parts.size() == 4);

	int i = atoi(parts[0]);
	int j = atoi(parts[1]);
	int o = string_to_orientation(parts[2]);
	int d = atoi(parts[3]);

	assert(0 <= i && i < fixed_p);
	assert(0 <= j && j < fixed_p);

	spaced_dimer_cases.push_back(boost::make_tuple(i, j, o, d));
	printf("%i: %s\n", a, v[a].c_str());
      }
      //spaced_dimer_p = spaced_dimer_cases.size();
    }
    */

    /*
    // Command line parameter disabled
    if (vm.count("fixed-lambdalist")) {
      std::vector<std::string> fixedlambdastringlist  = split(vm["fixed-lambdalist"].as< std::string >(), ',');
      for (int i=0; i < fixedlambdastringlist.size(); ++i) {
	fixed_lambda.push_back(atof(fixedlambdastringlist[i].c_str()));
	assert(fixed_lambda[i] >= 0.0 and fixed_lambda[i] <= 1.0);
      }
    } else {      // uniform distribution
    }
    assert(fixed_lambda.size() == fixed_p);
    */
    fixed_lambda.assign(fixed_p, 0.0);  // Will be assigned later

    if (vm.count("keep-monomer-fixed")) {
      std::string param = vm["keep-monomer-fixed"].as< std::string >();
      if (param == "all") 
	keep_monomer_fixed.assign(fixed_p, true);
      else {
	keep_monomer_fixed.assign(fixed_p, false);
	std::vector<std::string> param_list = split(param, ',');
	BOOST_FOREACH(std::string s, param_list) {
	  int tf = atoi(s);
	  if (tf < 0 || tf >= fixed_p) {
	    fprintf(stderr, "Parameter %s to option --keep-monomer-fixed has to be an integer in the range [0,%i)\n", s.c_str(), fixed_p);
	    exit(1);
	  }
	  keep_monomer_fixed[tf] = true;
	}
      }
    }
    else
      keep_monomer_fixed.assign(fixed_p, false);

    if (vm.count("epsilon")) {
      epsilon    = vm["epsilon"].as< double >();
      error(epsilon <= 0.0, "Epsilon must be positive real number.");
    }

    if (vm.count("extension-threshold")) {
      extension_threshold = vm["extension-threshold"].as< double >();
      allow_extension = true;
      error(extension_threshold <= 0.0, "Epsilon must be positive real number.");
    }

  }
  catch (std::exception& e) {
    std::cerr << e.what() << std::endl;
    std::cerr << desc << "\n";
    exit(1);
  }


  char hostname[200+1];
  hostname[200] = 0;
  gethostname(hostname, 200);
  time_t now;
  time(&now);
  printf("Starting program at %s", asctime(localtime(&now))); 
  printf("MODER version %s\n", PACKAGE_VERSION);
  printf("Running on host: %s\n", hostname);
  printf("Using %i openmp threads.\n", number_of_threads);
  
  if (combine_seeds_func == combine_seeds_OR)
    printf("Using OR combine function.\n");
  else if (combine_seeds_func == combine_seeds_AND)
    printf("Using AND combine function.\n");
  else if (combine_seeds_func == combine_seeds_masked)
    printf("Using masked (N) combine function.\n");
  else {
    error(true, "Error: combine_seeds_func not initialized.\n");
  }

  double background_lambda_fraction = 0.5;
  double PWM_lambda_fraction = cob_combinations.size() == 0 ? 0.5 : 0.3;
  double per_PWM_lambda_fraction = PWM_lambda_fraction / fixed_p;
  double dimer_lambda_fraction = 
    1.0 - (background_lambda_fraction + fixed_p * per_PWM_lambda_fraction);


  background_lambda = background_lambda_fraction;

  for (int k=0; k < fixed_p; ++k) {
    fixed_lambda[k] = per_PWM_lambda_fraction;
  }

  //  if (use_multinomial)
  //assert(seeds_given);
  
  int lines, bad_lines;


  std::vector<std::string> parts = split(seqsfile, '.');
  std::string extension = boost::to_lower_copy(parts.back());

  if (extension == "fasta" or extension == "fa") {
    printf("Reading from fasta file\n");
    boost::tie(lines, bad_lines) = read_fasta_sequences(seqsfile, sequences);
  }
  else {
    printf("Reading from sequence-per-line file\n");
    boost::tie(lines, bad_lines) = read_sequences(seqsfile, sequences);
  }
  //  for (int i=0; i < sequences.size(); ++ i)
  //    printf("%i: %s\n", i, sequences[i].c_str());
  
  check_data(sequences);
  printf("Read %zu good lines from file %s\n", sequences.size(), seqsfile.c_str());
  printf("Discarded %i bad lines\n", bad_lines);
  if (unique_param >= 0) {
    TIME_START(s1);
    sequences = remove_duplicate_reads_faster(sequences, unique_param);
    TIME_PRINT("\tRemoving Hamming duplicates took %.2f seconds.\n", s1);
  }
  printf("Using %zu sequences\n", sequences.size());
  std::vector<int> L(lines);
  int Lmin = std::numeric_limits<int>::max();
  int Lmax = std::numeric_limits<int>::min();
  for (int i=0; i < lines; ++i) {
    L[i] = sequences[i].length();
    if (L[i] < Lmin)
      Lmin = L[i];
    if (L[i] > Lmax)
      Lmax = L[i];
  }
  printf("Minimum sequence length is %i\n", Lmin);
  printf("Maximum sequence length is %i\n", Lmax);
  printf("Use markov model for background: %s\n", 
	 yesno(use_markov_background));
  printf("Use positional model for background: %s\n", 
	 yesno(use_positional_background));
  printf("Use two dna strands: %s\n", yesno(use_two_strands));
  printf("Avoid palindromes: %s\n", yesno(avoid_palindromes));
  printf("Use multinomial model: %s\n", yesno(use_multinomial));
  if (use_multinomial)
    printf("Hamming radius is %i\n", hamming_radius);  
  printf("Use pseudo counts: %s\n", yesno(use_pseudo_counts));
  printf("Get full flanks: %s\n", yesno(get_full_flanks));
  printf("Allow extension: %s\n", yesno(allow_extension));
  if (allow_extension)
    printf("extension threshold: %.2f\n", extension_threshold);
  printf("Minimum distance for learning is %i\n", minimum_distance_for_learning);
  printf("Minimum fraction for learning is %f\n", learning_fraction);
  printf("Information content threshold is %f\n", ic_threshold);
  printf("Cob cutoff constant is %f\n", cob_cutoff);
  
  
  printf("Epsilon is %g\n", epsilon);
  printf("Maximum number of iterations is %i\n", max_iter);
  std::vector<int> character_frequencies;
  
  boost::tie(background_frequencies, background_frequency_matrix, character_frequencies) = count_background(sequences);
  int CG = background_frequencies[1] + background_frequencies[2];
  printf("CG content: %lg\n", (double)CG/sum(background_frequencies));
  printf("Sequences contain %i characters\n", character_count);

  //if (use_markov_background) {
  printf("Background frequency matrix:\n");
  background_frequency_matrix.print(10);
  printf("Di-nucleotide count is %i\n", digram_count);
  background_probability_matrix = background_frequency_matrix;
  normalize_matrix_rows(background_probability_matrix);
  printf("Background probility matrix:\n");
  background_probability_matrix.print(10);
  assert(is_row_stochastic_matrix(background_probability_matrix));
    //}

  // normalize background vector
  background_probabilities = background_frequencies;
  normalize_vector(background_probabilities);
  
  printf("Empirical background probability: ");
  printf("[%g,%g,%g,%g]\n", 
	 background_probabilities[0],background_probabilities[1],
	 background_probabilities[2],background_probabilities[3]);

  printf("Using background probability: ");
  printf("[%g,%g,%g,%g]\n", background_probabilities[0],background_probabilities[1],
	 background_probabilities[2],background_probabilities[3]);

  if (prior_parameter == "addone")
    pseudo_counts.use_add_one(0.000001);
  else if (prior_parameter == "dirichlet")
    pseudo_counts.use_dirichlet(0.01, background_probabilities);

  if (use_pseudo_counts) {
    printf("Use prior: %s\n", prior_parameter.c_str());
    printf("\twith probabilities: ");
    std::vector<double> d = pseudo_counts.get();
    for (int i=0; i < 4; ++i)
      printf("%.4e ", d[i]);
    printf("\n");
  }

  //  positional_background = count_positional_background(sequences);

  //if (use_pseudo_counts)
  //  add_pseudo_counts(positional_background);
  normalize_matrix_columns(positional_background);
  assert(is_column_stochastic_matrix(positional_background));
  // assert(is_palindromic_matrix(positional_background));
  // write_matrix(stdout, positional_background, "Positional background probility matrix:\n", "%10lf");

  printf("Keep monomer fixed: %s\n", print_vector(keep_monomer_fixed).c_str());

  // initial motifs
  std::vector<dmatrix> fixed_M(fixed_p);

  if (seeds_given) {
    printf("Initial fixed seeds are %s\n", print_vector(fixedseedlist).c_str());

    for (int k=0; k < fixed_p; ++k) {
      error(fixedseedlist[k].length() > Lmin, "Seed is longer than the sequences");
      if (not use_meme_init) {
	fixed_M[k] = multinomial1_multimer_bernoulli_corrected(fixedseedlist[k], sequences);
	//fixed_M[k] = find_snips_multimer_helper(fixedseedlist[k], sequences).get<0>();
	write_matrix(stdout, fixed_M[k], to_string("Unnormalized initial fixed matrix %i from seed %s:\n", 
						   k, fixedseedlist[k].c_str()), "%.6f");
	if (use_pseudo_counts)
	  pseudo_counts.add(fixed_M[k]);
	normalize_matrix_columns(fixed_M[k]);
	if (use_multinomial and adjust_seeds)
	  fixedseedlist[k] = string_giving_max_score(fixed_M[k]);
      }
      else
	fixed_M[k] = get_meme_init_pwm(fixedseedlist[k]);
      write_matrix(stdout, fixed_M[k], to_string("Initial fixed matrix %i from seed %s:\n", 
						 k, fixedseedlist[k].c_str()), "%.6f");
    }
  } else {
    fixedseedlist.resize(fixed_p);
    for (int k=0; k < fixed_p; ++k) {
      fixed_M[k] = read_matrix_file(fixedmatrixfilelist[k]);
      error(fixed_M[k].get_columns() > Lmin, "Matrix is wider than the sequences");
      if (use_pseudo_counts)
	pseudo_counts.add(fixed_M[k]);
      normalize_matrix_columns(fixed_M[k]);
      write_matrix(stdout, fixed_M[k], to_string("Initial matrix %i from file %s:\n", k, fixedmatrixfilelist[k].c_str()), "%.6f");
      assert(is_column_stochastic_matrix(fixed_M[k]));
      fixedseedlist[k] = string_giving_max_score(fixed_M[k]);
      //if (use_multinomial)
      //	assert(M[k].get_columns() == seedlist[k].length());
    }
    printf("Initial fixed seeds are %s\n", print_vector(fixedseedlist).c_str());

    printf("Read %i matrices\n", fixed_p);
    
  }



  gapped_kmer_context my_gapped_kmer_context(sequences);

  if (outputdir != ".")
    mkdir(outputdir.c_str(), S_IRWXU);
  //pwm_list_to_logos(to_string("%s/init.png", outputdir.c_str()), M);

  if (cob_combinations.size() > 0)
    dimer_lambda_fraction /= cob_combinations.size();

    
  std::vector<cob_params_t> my_cob_params;
  int i=0;
  BOOST_FOREACH(cob_combination_t cob_combination, cob_combinations) {
    int tf1, tf2;
    boost::tie(tf1, tf2) = cob_combination;
    int dmin, dmax;
    int max_dist_for_deviation;   // The deviation tables will be learned for distances in [dmin,max_dist_for_deviation]
    int k1 = fixedseedlist[tf1].length();
    int k2 = fixedseedlist[tf2].length();
    int mink = std::min(k1, k2);

    if (dmin_given)
      dmin = dmina[i];
    else
      dmin = -mink/2;

    if (dmax_given)
      dmax = dmaxa[i];
    else
      dmax = global_dmax;

    if (max_dist_for_deviation_given)
      max_dist_for_deviation = max_dist_for_deviationa[i];
    else
      max_dist_for_deviation = global_max_dist_for_deviation;
    
    // Make sure these are within reasonable limits
    dmin = std::max(-mink + 1, dmin);
    //    dmax = std::min(dmax, Lmax - k1 - k2);
    dmax = std::min(dmax, Lmin - k1 - k2);
    error(dmax < dmin, "Requested dimeric cases do not fit into the input sequence");
    max_dist_for_deviation = std::min(max_dist_for_deviation, dmax);
    cob_params_t cp = create_cob(cob_combination, fixedseedlist, fixed_M, dimer_lambda_fraction, 
				 overlapping_dimer_cases,
				 spaced_dimer_cases, sequences, L, my_gapped_kmer_context, dmin, dmax, max_dist_for_deviation);
    cp.update_oriented_matrices(fixed_M, fixedseedlist);
    cp.compute_expected_matrices(fixed_M);
    cp.compute_deviation_matrices();
    printf("Maximum distance for gap learning for cob case %s is %i\n", cp.name().c_str(), max_dist_for_deviation);

    for (int o=0; o < cp.number_of_orientations; ++o) {
      for (int d=cp.dmin; d < std::min(0, cp.dmax+1); ++d) {
	write_matrix(stdout, cp.deviation[o][d], 
		       to_string("Initial deviation matrix %s %s %i:\n", cp.name().c_str(),
				 orients[o], d).c_str(), "%.6f");
      }
    }

    my_cob_params.push_back(cp);
    ++i;
  }

  double dimer_lambdas_sum = 0.0;
  for (int r=0; r < my_cob_params.size(); ++r)
    dimer_lambdas_sum += sum(my_cob_params[r].dimer_lambdas);
  double temp_sum = sum(fixed_lambda) + dimer_lambdas_sum + background_lambda;
  assert(fabs(temp_sum - 1.0) < 0.001);



  printf("Initial fixed lambdas are %s\n", print_vector(fixed_lambda).c_str());  

  if (use_dimers) {
    for (int r=0; r < my_cob_params.size(); ++r) {
      print_cob(stdout, my_cob_params[r].dimer_lambdas, 
		to_string("Initial dimer lambdas %s:\n", my_cob_params[r].name().c_str()), "%.5f");
      printf("Initial sum of dimer lambdas %s is %.2f\n", my_cob_params[r].name().c_str(), sum(my_cob_params[r].dimer_lambdas));
    }
  }
  
  printf("Initial background lambda is %.4f\n", background_lambda);


  // run the algorithm
  std::vector<dmatrix> res=
    multi_profile_em_algorithm(sequences, 
			       fixed_M,
			       keep_monomer_fixed,
			       background_probability_matrix, background_probabilities, 
			       background_lambda, 
			       fixed_lambda, fixedseedlist, 
			       my_cob_params,
			       epsilon, extension_threshold, my_gapped_kmer_context);
  
  TIME_PRINT("Whole program took %.2f seconds cpu-time\n", t);
  WALL_TIME_PRINT("Whole program took %.2f seconds wall-time\n", t2);

  if (false)
    SSS("dummy");   // Just make sure that code for SSS is included in the executable. For debugging purpose.
  return 0;
}
