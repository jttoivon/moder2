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
#define TIMING
#include "timing.hpp"

#include "multinomial_helper.hpp"
#include "bndm.hpp"
#include "common.hpp"
#include "parameters.hpp"
#include "iupac.hpp"
#include "matrix_tools.hpp"


#include <boost/foreach.hpp>

extern bool use_palindromic_correction;

// Do not reject sequences with multiple occurrences of query strings.
// Compute the counts for the multinomial1 matrix
boost::tuple<dmatrix,int,int>
find_snips_multimer_helper(const std::string& consensus, const std::vector<std::string>& sequences)
{
  
  std::string str1;
  std::string str2;

  //typedef boost::tuple<int,int,int> triple;  // (seqno, position, direction)
  //std::vector<triple> alignment;
  std::vector<std::string> alignment;

  //  int lines = sequences.size();

  str1=join(sequences, '#');
  str2=join_rev(sequences, '#');

  int k = consensus.length();
  char nucs[] = "ACGT";
  bool is_palindrome = is_palindromic(consensus);
  int consensus_count = BNDM_with_joker(str1, consensus);
  if (use_two_strands && (count_palindromes_twice || not is_palindrome))
    consensus_count += BNDM_with_joker(str2, consensus);

  if (print_alignment) {
    for (int t=0; t < consensus_count; ++t)
      alignment.push_back(consensus);
  }

  matrix<double> result(4, k);
  for (int j=0; j < k; ++j) {        // iterate through all string positions
    std::string temp = consensus;
    for (int a=0; a < 4; ++a) {      // iterate through all characters
      if (nucs[a] == consensus[j])   // the consensus count was already computed. NOTE: no iupac_match here, on purpose
	continue;
      temp[j]=nucs[a];
      is_palindrome = is_palindromic(temp);
      
      result(a,j) = BNDM_with_joker(str1,temp);
      if (use_two_strands && (count_palindromes_twice || not is_palindrome))
	result(a,j) += BNDM_with_joker(str2,temp);

      if (print_alignment) {
	for (int t=0; t < result(a,j); ++t)
	  alignment.push_back(temp);
      }
    }
  }


  int total_count = consensus_count;
  for (int j=0; j < k; ++j) {       // iterate through all string positions
    for (int a=0; a < 4; ++a) {     // iterate through all characters
      if (nucs[a] == consensus[j]) {
	result(to_int(consensus[j]), j) = consensus_count;   // add consensus count to the matrix
	if (print_alignment) {
	  for (int t=0; t < consensus_count; ++t)              // add consensus to alignment also
	    alignment.push_back(consensus);                    // for other columns
	}
      }
      else if (not iupac_match(nucs[a], consensus[j])) // the consensus count was already added to the total_count
	total_count += result(a, j);   // number of sequences used for the matrix
    }
  }


  // print the alignment to file descriptor 3, if it is open
  if (print_alignment) {
    FILE* fp = fdopen(3, "a");
    if (fp != NULL) {
      for (int t = 0; t < alignment.size(); ++t)
	fprintf(fp, "%s\n", alignment[t].c_str());
      fclose(fp);
    }
  }

  return boost::make_tuple(result, consensus_count, total_count);
} // find_snips_multimer


dmatrix
find_snips_multimer(const std::string& consensus, const std::vector<std::string>& sequences, int hamming_distance)
{
  TIME_START(t);
  assert(hamming_distance == 1);
  dmatrix result;
  int seed_count;
  int multinomial_count;
  boost::tie(result, seed_count, multinomial_count) = find_snips_multimer_helper(consensus, sequences);
  TIME_PRINT("Multinomial-1 algorithm took %.2f seconds.\n", t);
  printf("Seed count = %i\n", seed_count);
  printf("Total multinomial1 count is %d\n", multinomial_count);
  return result;
}

string_to_tuple_type
get_n_neighbourhood(const std::string&seed, int n)
{
  const int k = seed.length();
  //const int L = sequences[0].length();
  assert(n >= 0);
  assert(n <= k);
  char nucs[] = "ACGT";

  //code_to_tuple_type code_to_tuple;
  string_to_tuple_type string_to_tuple;

  std::vector<std::string> complements(k);  // These are set complements, not nucleotide complements
  std::vector<int> bases(k);
  unsigned long long N_mask=0;                // bitmask for positions that contain 'N'. Those positions cannot contain an error
  for (int i=0; i < k; ++i) {
    complements[i] = complement_set(seed[i]);
    bases[i]=complements[i].length() - 1;
    if (seed[i]=='N')
      N_mask |=  (1 << (k-1-i));
  }
  for (int j=0; j < k; ++j) {        // iterate through all string positions
    std::string temp = seed;         // Variable temp is not really needed. It's just for sanity check.
    //packed_string myid2(seed);
    for (int a=0; a < 4; ++a) {      // this handles hamming distances 0 and 1
      if (n == 0 and not iupac_match(nucs[a], seed[j]))
	continue;
      temp[j]=nucs[a];
      // myid2[j] = a;
      // my_assert(myid2.get_bits(), dna_to_number(temp));
      // code_to_tuple[myid2.get_bits()].push_back(boost::make_tuple(j, a));
      string_to_tuple[temp].push_back(boost::make_tuple(j, a));
    }

    for (int error=1; error < n; ++error) { // errors outside position j, handles hamming distances 2 <= hd <= n
      // bitvector c has 1-bit for each member of the subset, rightmost bit is bit number k-1
      unsigned long long c = (1ull<<error)-1;  // initially rightmost 'error' bits are 1
      int mycount = 0;
      // iterate through all subsets c of {0, ..., k-1} that have size 'error'
      while (c < (1ull<<k)) {   // Superset has only k elements
	assert(__builtin_popcountll(c) == error);
	if (((c & (1ull << (k-1-j)))) == 0 and ((c & N_mask) == 0))  { // j doesn't belong to the subset, and subset positions don't contain 'N'
	  ++mycount;
	  std::vector<int> P;  // positions that contain error
	  int number_of_combinations = 1;
	  std::string temp = seed;
	  for (int pos=0; pos < k; ++pos) {
	    if ((c & (1ull << (k-1-pos))) != 0) {   // pos belongs to set c
	      P.push_back(pos);
	      temp[pos] = complements[pos][0];  // initialize to first string that has mismatches at these positions
	      number_of_combinations *= complements[pos].length();
	    }
	  }
	  //packed_string myid4(seed);
	  
	  std::vector<int> y(error, 0);
	  y[error-1]=-1;  // Initialize
	  for (int r=0; r < number_of_combinations; ++r) {
      
	    int i;
	    for (i=error-1; y[i] == bases[P[i]]; --i) {
	      y[i]=0;
	      // temp[P[i]] = nucs[myskip(y[i], skip[i])];
	      // myid4[P[i]] = myskip(y[i], skip[i]);
	      temp[P[i]] = complements[P[i]][y[i]];
	    }
	    y[i]++;
	    temp[P[i]] = complements[P[i]][y[i]];
	    // temp[P[i]] = nucs[myskip(y[i], skip[i])];
	    // myid4[P[i]] = myskip(y[i], skip[i]);

	    for (int a=0; a < 4; ++a){
	      temp[j] = nucs[a];
	      //myid4[j] = a;
	      //my_assert(myid4.get_bits(), dna_to_number(temp));
	      //code_to_tuple[myid4.get_bits()].push_back(boost::make_tuple(j, a));
	      string_to_tuple[temp].push_back(boost::make_tuple(j, a));
	    }
	  }  // end for r

	}    // end if j not in c

	unsigned long long a = c&-c;
	unsigned long long b = c+a;   // update bitvector c. This is "Gosper's hack"
	c = (c^b)/4/a|b;

      } // end foreach subset c


    }  // end for error
  } // end for j

  return string_to_tuple;
}

// Returns a vector with the seed pattern at the first index
std::vector<std::pair<std::string, std::vector<boost::tuple<int, int> > > >
get_n_neighbourhood_in_vector(const std::string&seed, int n)
{
  std::vector<std::pair<std::string, std::vector<boost::tuple<int, int> > > > result;
  string_to_tuple_type neigh = get_n_neighbourhood(seed, n);
  std::pair<std::string, std::vector<boost::tuple<int, int> > > t;
  string_to_tuple_type::iterator it = neigh.find(seed);   // Make sure that the pair corresponding to the seed is first on the vector
  result.push_back(*it);
  neigh.erase(it);
  BOOST_FOREACH(t, neigh) {
     result.push_back(t);
  }
  return result;
}

class iupac_probability_in_background
{
public:
  iupac_probability_in_background(const std::vector<double>& bg)
    : iupac_probabilities(256)
  {
    BOOST_FOREACH(char iupac_char, iupac_chars) {
      BOOST_FOREACH(char c, iupac_class(iupac_char)) {
	iupac_probabilities[(unsigned char)iupac_char] += bg[to_int(c)];
      }
    }
  }

  double
  operator()(const std::string& s) {
    assert(is_iupac_string(s));
    double result = 1.0;
    for (int i=0; i < s.length(); ++i)
      result *= iupac_probabilities[s[i]];
    return result;
  }
  
private:
  std::vector<double> iupac_probabilities;
};

double
palindromic_correction(const std::string& pattern, const std::string& seed, const std::string& seed_rev)
{
  double t=1.5;
  int hd1=hamming_distance(pattern,seed);
  int hd2=hamming_distance(pattern, seed_rev);
  double correction = pow(t, -hd1) /
    (pow(t, -hd1) + pow(t, -hd2));

  return correction;
}

boost::tuple<dmatrix,int>
find_multinomial_n_background(const std::string& seed, const std::vector<std::string>& sequences, const std::vector<double>& bg,
			      int n, bool use_multimer)
{
  const int k = seed.length();
  const int L = sequences[0].length();
  assert(n >= 0);
  assert(n <= k);
  dmatrix result(4, k);

  string_to_tuple_type string_to_tuple;
  string_to_tuple = get_n_neighbourhood(seed, n);

  //printf("Number of patterns %lu\n", string_to_tuple.size());
  unsigned long seed_count=0;
  unsigned long total_count=0;
  int lines = sequences.size();
  int sites = lines * (L-k+1)*2;
  iupac_probability_in_background iupac_prob(bg);
  
  if (use_multimer) {
    seed_count = sites * iupac_prob(seed);
    total_count += seed_count;

    std::vector<boost::tuple<int, int> > pairs;
    std::string pattern;
    std::string seed_rev = reverse_complement(seed);
    BOOST_FOREACH(boost::tie(pattern, pairs), string_to_tuple) {
      unsigned long count = sites * iupac_prob(pattern);
      if (use_palindromic_correction)
	count *= palindromic_correction(pattern, seed, seed_rev);
      
      if (not iupac_string_match(pattern, seed))
	total_count += count;
      int j, a;
      BOOST_FOREACH(boost::tie(j,a), pairs) {
	result(a, j) += count;
      }
    }
  } else {  // allow only single occurrence per read
    error(true, "Not implemented");
  }
  return boost::make_tuple(result, total_count);
} // find_multinomial_n_background



boost::tuple<dmatrix,int>
find_multinomial_n_suffix_array(const std::string& seed, const std::vector<std::string>& sequences, const suffix_array& sa, int n, bool use_multimer)
{
  const int k = seed.length();
  const int L = sequences[0].length();
  assert(n >= 0);
  assert(n <= k);
  //char nucs[] = "ACGT";
  dmatrix result(4, k);
  //code_to_tuple_type code_to_tuple;

  string_to_tuple_type string_to_tuple;
  string_to_tuple = get_n_neighbourhood(seed, n);

  printf("Number of patterns %lu\n", string_to_tuple.size());
  unsigned long seed_count=0;
  unsigned long total_count=0;
  //unsigned long number_of_sites=0;

  std::string seed_rev = reverse_complement(seed);
  if (use_multimer) {
    seed_count = sa.count_iupac(seed);
    total_count += seed_count;

    std::vector<boost::tuple<int, int> > pairs;
    std::string pattern;
    printf("#String\tColumn\tHamming distance\tPalindrome\tCount\tMatches at col\n");
    BOOST_FOREACH(boost::tie(pattern, pairs), string_to_tuple) {
      unsigned long count = sa.count_iupac(pattern);
      bool is_palindrome = is_palindromic(pattern);
      int hd = hamming_distance(seed, pattern);
      if (is_palindrome and use_two_strands and not count_palindromes_twice)
	count /= 2;
      if (use_palindromic_correction)
	count *= palindromic_correction(pattern, seed, seed_rev);

      if (not iupac_string_match(pattern, seed))
	total_count += count;
      int j, a;
      BOOST_FOREACH(boost::tie(j,a), pairs) {
	printf("#%s\t%i\t%i\t%s\t%zu\t%s\n", pattern.c_str(), j, hd,
	       yesno(is_palindrome), count, yesno(seed[j]==pattern[j]));
	result(a, j) += count;
      }
    }
  } else {  // allow only single occurrence per read
    typedef string_to_tuple_type::iterator iterator;
    int lines = sequences.size();
    std::vector<std::set<int> > hit_positions(lines);
    std::vector<std::vector<iterator> > hit_patterns(lines);
    int divisor = L + 1;    // Includes the separator '#'
    std::vector<boost::tuple<int, int> > pairs;
    std::string pattern;
    for (iterator it=string_to_tuple.begin(); it != string_to_tuple.end(); ++it) {
      std::vector<long int> positions;
      pattern = it->first;
      bool is_palindrome = is_palindromic(pattern);
      sa.locate_iupac(pattern, positions);
      BOOST_FOREACH(int pos, positions) {
	int i = pos / divisor;  // index of the read containing the pos
	int j = pos % divisor;
	if (i >= lines) {       // handle reverse complement
	  i = 2*lines - i - 1;
	  j = L - (j + k - 1) - 1;
	}
	if (hit_positions[i].count(j) < 1 or not is_palindrome or count_palindromes_twice) {
	  hit_positions[i].insert(j);
	  //	  if (hits[i] == 1)
	  hit_patterns[i].push_back(it);
	}
      }
    }
    for (int i=0; i < lines; ++i) {
      if (hit_positions[i].size() != 1)                    // Because hit_positions[i] is a set, this doesn't exclude, for instance, palindromes
	continue;
      for (int t=0; t < hit_patterns[i].size(); ++t) {
	boost::tie(pattern, pairs) = *(hit_patterns[i][t]);
	if (not iupac_string_match(pattern, seed))
	  total_count += 1;
	else
	  ++seed_count;
	int j, a;
	BOOST_FOREACH(boost::tie(j,a), pairs) {
	  result(a, j) += 1;
	}
      }
    }
    // for (int a=0; a < 4; ++a)             // get the seed count                      // THIS PROBABLY ISN'T CORRECT, CONTAINS ALSO EXTRA COUNTS, WORKS ONLY WHEN n=1
    //   if (iupac_match(nucs[a], seed[0]))
    // 	seed_count += result(a,0);
    total_count += seed_count;
  }

  printf("Seed %s count = %lu\n", seed.c_str(), seed_count);
  printf("Total multinomial-n count is %lu\n", total_count);

  return boost::make_tuple(result, total_count);
}


dmatrix
align_all(const std::vector<std::string>& sequences)
{
  int L=sequences[0].length();
  int lines = sequences.size();
  int k = L;
  dmatrix result(4, k);
  for (int i=0; i < lines; ++i) {
    const std::string& line = sequences[i];
    for (int j=0; j < k; ++j) 
      result(to_int(line[j]), j) += 1;
    // if (use_two_strands) {
    //   const std::string& line_rev = reverse_complement(sequences[i]);
    //   for (int j=0; j < k; ++j) 
    // 	result(to_int(line_rev[j]), j) += 1;
    // }
  }

  return result;
}
