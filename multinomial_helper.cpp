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
#include "multinomial_helper.hpp"
#include "bndm.hpp"
#include "common.hpp"
#include "parameters.hpp"
#include "iupac.hpp"
#include "matrix_tools.hpp"

#define TIMING 1
#include "timing.hpp"

#include <boost/foreach.hpp>



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

  int consensus_count = BNDM_with_joker(str1, consensus);
  if (use_two_strands)
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
      
      result(a,j) = BNDM_with_joker(str1,temp);
      if (use_two_strands)
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
  TIME_PRINT("Multinomial-1 algorithm tooks %.2f seconds.\n", t);
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


boost::tuple<dmatrix,int>
find_multinomial_n_suffix_array(const std::string& seed, const std::vector<std::string>& sequences, const suffix_array& sa, int n, bool use_multimer)
{
  const int k = seed.length();
  const int L = sequences[0].length();
  assert(n >= 0);
  assert(n <= k);
  char nucs[] = "ACGT";
  dmatrix result(4, k);
  //code_to_tuple_type code_to_tuple;

  string_to_tuple_type string_to_tuple;
  string_to_tuple = get_n_neighbourhood(seed, n);

  printf("Number of patterns %lu\n", string_to_tuple.size());
  unsigned long seed_count=0;
  unsigned long total_count=0;
  //unsigned long number_of_sites=0;

  if (use_multimer) {
    seed_count = sa.count_iupac(seed);
    total_count += seed_count;

    std::vector<boost::tuple<int, int> > pairs;
    std::string pattern;
    BOOST_FOREACH(boost::tie(pattern, pairs), string_to_tuple) {
      unsigned long count = sa.count_iupac(pattern);
      if (not iupac_string_match(pattern, seed))
	total_count += count;
      //printf("Pattern %s has count %lu\n", pattern.c_str(), count);
      for (int i=0; i < pairs.size(); ++i) {
	int j, a;
	boost::tie(j,a) = pairs[i];
	// if (j==4 and count != 0)
	//   printf("#%s %lu\n", pattern.c_str(), count);
	result(a, j) += count;
      }
    }
  } else {  // allow only single occurrence per read
    typedef string_to_tuple_type::iterator iterator;
    int lines = sequences.size();
    std::vector<std::set<int> > hit_positions(lines);
    std::vector<std::vector<iterator> > hit_patterns(lines);
    int divisor = L + 1;
    std::vector<boost::tuple<int, int> > pairs;
    std::string pattern;
    for (iterator it=string_to_tuple.begin(); it != string_to_tuple.end(); ++it) {
      std::vector<long int> positions;
      pattern = it->first;
      sa.locate_iupac(pattern, positions);
      BOOST_FOREACH(int pos, positions) {
	int i = pos / divisor;  // index of the read containing the pos
	int j = pos % divisor;
	if (i >= lines) {       // handle reverse complement
	  i -= lines;
	  j = L - j - 1;
	}
	hit_positions[i].insert(j);
	//	  if (hits[i] == 1)
	hit_patterns[i].push_back(it);
      }
    }
    for (int i=0; i < lines; ++i) {
      if (hit_positions[i].size() != 1)
	continue;
      for (int t=0; t < hit_patterns[i].size(); ++t) {
	boost::tie(pattern, pairs) = *(hit_patterns[i][t]);
	if (not iupac_string_match(pattern, seed))
	  total_count += 1;
	for (int s=0; s < pairs.size(); ++s) {
	  int j, a;
	  boost::tie(j,a) = pairs[s];
	  result(a, j) += 1;
	}
      }
    }
    for (int a=0; a < 4; ++a)             // get the seed count
      if (iupac_match(nucs[a], seed[0]))
	seed_count += result(a,0);
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
