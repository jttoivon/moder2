#define TIMING
#include "timing.hpp"

#include "multinomial_helper.hpp"
#include "bndm.hpp"
#include "common.hpp"
#include "parameters.hpp"
#include "probabilities.hpp"
#include "iupac.hpp"
#include "huddinge.hpp"
#include "matrix_tools.hpp"
#include "kmer_tools.hpp"


#include <boost/foreach.hpp>

int extra_flank;   // This is only for testing. Remove later


//counting_type data_counting = all_occurrences;
counting_type data_counting = neighbourhood_contains_one;
counting_type background_counting = data_counting;

//extern bool use_palindromic_correction;

int palindromic_index_limit=0;
//bool use_palindromic_correction = false;
int low_count_limit = 20;
int cluster_threshold=4;
bool use_one_per_cluster=true;

int
conflict_free_palindromic_index(int hamming_radius)
{
  int limit;
  switch (hamming_radius) {
  case 0: limit = 0; break;
  case 1: limit = 2; break;
  default: limit = 2*hamming_radius + 1;
  }
  return limit;
}

std::vector<std::string>
remove_masked_areas(const std::vector<std::string>& sequences, int k)
{
  std::vector<std::string> result;
  int lines = sequences.size();
  for (int i=0; i < lines; ++i) {
    const std::string& line = sequences[i];
    std::string current;
    for (int j=0; j < line.length(); ++j) {
      if (line[j] == 'N') {
	if (current.length() >= k)
	  result.push_back(current);
	if (current.length() > 0)
	  current.clear();
      }
      else
	current.push_back(line[j]);
    }
    if (current.length() >= k)
      result.push_back(current);
      
  }
  return result;
}

// The triples are of (code, count, palindromic index) type 
// Primary sort key is the count, the secondary sort key is the palindromic index
bool
triple_comp(boost::tuple<big_int, int, int> a, boost::tuple<big_int, int, int> b)
{
  return a.get<1>() > b.get<1>() or (a.get<1>() == b.get<1>() and a.get<2>() > b.get<2>());
}

boost::tuple<std::string,int>
most_common_pattern_multimer(const std::vector<std::string>& sequences, int k, std::string seed,
			     bool contains_N, int hamming_radius)
{
  assert(k > 1 && k <= max_matrix_len); 
  //int strings = (int)pow(4,k); // if k==12 then this is about 16M

  //  std::vector<unsigned> number_of_occurrences(strings);
  count_container_t number_of_occurrences;
  bool count_palindromes_twice = use_two_strands;  // There is also global version of this
  if (contains_N)
    get_kmer_counts(remove_masked_areas(sequences, k), k, number_of_occurrences, use_two_strands, count_palindromes_twice);
  else
    get_kmer_counts(sequences, k, number_of_occurrences, use_two_strands, count_palindromes_twice);

  
  printf("Size of number_of_occurrences container %zu\n", number_of_occurrences.size());

  //count_t count;
  //code_t code=0;
  /*
  BOOST_FOREACH(boost::tie(code, count), number_of_occurrences)
    printf("*%llu\t%s\t%i\n", (unsigned long long) code, number_to_dna(code,k).c_str(), count);
  
  printf("Size of number_of_occurrences container %zu\n", number_of_occurrences.size());
  */
  
  code_t argmax;
  std::vector<boost::tuple<big_int, int, int> > v;
  code_t code;
  int count;
  BOOST_FOREACH(boost::tie(code, count), number_of_occurrences)
    v.push_back(boost::make_tuple(code, count, palindromic_index(number_to_dna(code, k))));
  std::sort(v.begin(), v.end(), triple_comp);   // Primary sort key is the count, the secondary sort key is the palindromic index

  // These are the count and pi of the most common kmer
  int top_count;
  int top_pi;
  boost::tie(boost::tuples::ignore, top_count, top_pi) = v[0];

  double ratio_cutoff = 0.25;   // Drop of two units in PI can be allowed to drop the corresponding count into one quarter
  if (palindromic_index_limit > 0) {
    // the palindromic index of a seed needs to be at least 'limit' for the n-Hamming-neighbourhood to be conflict-free
    code_t code=0;
    code_t max_code = 0;
    int max_pi = -1;
    std::vector<int> alignments;
    std::string maxseed =  number_to_dna(v[0].get<0>(), k);   // Seed with highest count
    std::string maxseed_revcomp = reverse_complement(maxseed);
    for (int i=0; i < v.size(); ++i) {
      code = v[i].get<0>();
      int count = v[i].get<1>();
      int pi = v[i].get<2>();
      if (count < low_count_limit and max_pi >= 0)   // If we have already one candidate and the current kmer is too small, then quit.
	break;
      std::string temp = number_to_dna(code, k);
      if (huddinge_distance(maxseed, temp) <= huddinge_distance(maxseed_revcomp, temp))
	alignments = huddinge_alignment(maxseed, temp);
      else
	alignments = huddinge_alignment(maxseed_revcomp, temp);
      
      // middle part of the condition checks that candidate is not a shift of maxseed (or its reverse complement)
      if (pi > max_pi and alignments.size() == 1 and alignments[0] == 0 and (float)count/top_count >= pow(ratio_cutoff, (pi - top_pi)/2)) {   
	max_pi = pi;
	max_code = code;
	printf("!New candidate for seed: %s %i %i\n", temp.c_str(), count, pi);
      }
      if (max_pi >= palindromic_index_limit)
	break;                   // found a good seed
    }
    argmax = max_code;
  }
  else {
    argmax = v[0].get<0>();
  }
  
  //  int max_count = number_of_occurrences[argmax];
  std::string result;

  // between string and its reverse complement, choose lexicographically smaller
  if (seed.length() == k) {
    result = seed;
  }
  else {
    std::string result1 = number_to_dna(argmax,k);
    std::string result2 = reverse_complement(result1);
    if (use_two_strands)
      result = (result1 < result2) ? result1 : result2;
    else
      result = result1;
  }

  printf("Seed %s has %i occurences\n",
	 result.c_str(),number_of_occurrences[dna_to_number(result)]);

  return boost::make_tuple(result, number_of_occurrences[dna_to_number(result)]);
  //return result;
}


boost::tuple<std::string,int>
most_common_pattern_monomer(const std::vector<std::string>& sequences, int k, std::string seed,
			    int hamming_radius)
{
  assert(k > 1 && k <= max_matrix_len); 
  //  int strings = (int)pow(4,k); // if k==12 then this is about 16M
  //std::vector<int> number_of_occurences(strings);
  //typedef std::map<big_int, int> my_container;
  typedef boost::unordered_map<big_int, int> my_container; // items are (code,count) pairs
  my_container number_of_occurrences;
  int lines = sequences.size();
  int max_count=-1;
  big_int argmax=-1;
  big_int id;
  big_int id2;
  bool count_palindromes_twice = use_two_strands;  // There is also global version of this
  for  (int i=0; i < lines; ++i) {
    const std::string& line = sequences[i];
    //std::map<big_int, std::vector<int> >occurences;

    boost::unordered_map<big_int, std::vector<int> > occurrences;  // occurrences on this line
    std::set<big_int> ids;

    // find all occurrences in sequence
    for (int j=0; j < line.length()-k+1; ++j) {
      id = dna_to_number(line.substr(j,k));
      occurrences[id].push_back(j);
      ids.insert(id);
      if (use_two_strands) {
	id2 = dna_to_number(reverse_complement(line.substr(j,k)));
	occurrences[id2].push_back(j);
      }
    }
    // accept only subsequences that appear only once per sequence, or is a palindromic occurrence
    for (std::set<big_int>::iterator it=ids.begin(); it != ids.end(); ++it) {
      big_int id = *it;
      std::vector<int>& r = occurrences[id];
      if (r.size() == 1 || (count_palindromes_twice && r.size() == 2 && r[0] == r[1])) {
	big_int id2 = reverse_complement_2bitstring(id, k);
	//	std::string s = number_to_dna(id, k);
	//	big_int id2 = dna_to_number(reverse_complement(s));
	++number_of_occurrences[id];
	if (use_two_strands)
	  ++number_of_occurrences[id2];
	if (number_of_occurrences[id] > max_count) {
	  max_count = number_of_occurrences[id];
	  argmax = id;
	}
	if (use_two_strands && number_of_occurrences[id2] > max_count) {
	  max_count = number_of_occurrences[id2];
	  argmax = id2;
	}	  
      }

    }
  }

  // print top10 of strings
  // triples are (code, count, palindromic index)
  std::vector<boost::tuple<big_int, int, int> > v;
  my_container::iterator it;
  for (it=number_of_occurrences.begin(); it != number_of_occurrences.end(); ++it) {
    v.push_back(boost::make_tuple(it->first, it->second, palindromic_index(number_to_dna(it->first, k))));
  }
  std::sort(v.begin(), v.end(), triple_comp);   // compares according the second member of the pair: the count
  if (palindromic_index_limit > 0) {
    code_t code=0;
    code_t max_code = 0;
    int max_pi = -1;
    std::vector<int> alignments;
    std::string maxseed =  number_to_dna(v[0].get<0>(), k);   // Seed with highest count
    std::string maxseed_revcomp = reverse_complement(maxseed);

    for (int i=0; i < v.size(); ++i) {
      code = v[i].get<0>();
      int pi = v[i].get<2>();
      int count = v[i].get<1>();
      if (count < low_count_limit and max_pi >= 0)   // If we have already one candidate and the current kmer is too small, then quit.
	break;

      std::string temp = number_to_dna(code, k);
      if (huddinge_distance(maxseed, temp) <= huddinge_distance(maxseed_revcomp, temp))
	alignments = huddinge_alignment(maxseed, temp);
      else
	alignments = huddinge_alignment(maxseed_revcomp, temp);
      
      // latter part of the condition checks that candidate is not a shift of maxseed (or its reverse complement)
      if (pi > max_pi and alignments.size() == 1 and alignments[0] == 0) {
	max_pi = pi;
	max_code = code;
      }
      if (max_pi >= palindromic_index_limit)
	break;                   // found a good seed
    }
    argmax = max_code;
  } else
    argmax = v[0].get<0>();

  // between string and its reverse complement, choose lexicographically smaller
  std::string result;
  std::string result1;
  if (seed.length() == k)
    result = seed;
  else {
    result1 = number_to_dna(argmax,k);
    std::string result2 = reverse_complement(result1);
    if (use_two_strands)
      result = (result1 < result2) ? result1 : result2;
    else
      result = result1;
  }
  printf("Seed %s has %i occurences\n",
	 result.c_str(),number_of_occurrences[dna_to_number(result)]);

  return boost::make_tuple(result, number_of_occurrences[dna_to_number(result)]);
}



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


class min_hamming_distance_class
{
public:

  min_hamming_distance_class () {}
  
  min_hamming_distance_class(const std::string& seed)
  {
    init(seed);
  }

  void
  init(const std::string& seed)
  {
    int k=seed.length();
    code_t s = dna_to_number(seed);
    code_t s_rev_comp = reverse_complement_2bitstring(s, k);
    unsigned int number_of_sequences = pow(4, k);
    v.resize(number_of_sequences);
    for (code_t code=0; code < number_of_sequences; ++code) {
      v[code] = std::min(hamming_distance_with_bits(code, s), hamming_distance_with_bits(code, s_rev_comp));
    }
  }
  
  unsigned char
  get(code_t code) const {
    return v[code];
  }
  
private:
  std::vector<unsigned char> v;
};


class cluster_probability_array {
public:

  typedef boost::multi_array<double, 2>::extent_range range;

  cluster_probability_array() {}
  
  // d is the hamming radius, e is the smallest Hamming distance to the seed in the cluster.
  // e must be atteined only in the end of the cluster.
  cluster_probability_array(const std::string& seed, int d_, int e_, const std::vector<double>& q_, int lmax, const min_hamming_distance_class& f_)
  //    : d(d_), e(e_), s(dna_to_number(seed)), q(q_), f(f_)
  {
    init(seed, d_, e_, q_, lmax, f_);
  }
  
  void
  init(const std::string& seed, int d_, int e_, const std::vector<double>& q_, int lmax, const min_hamming_distance_class& f_)
  {
    d = d_;
    e = e_;
    s = dna_to_number(seed);
    q = q_;
    f = f_;
    k = seed.length();
    epsilon = cluster_threshold;
    imax = lmax - (k-epsilon);
    number_of_sequences = pow(4, k);
    a.resize(boost::extents[number_of_sequences][range(k, imax+1)][k-epsilon+1]);
    mask = 1;
    mask = (mask << (2*k)) - 1;

    initialize();
    compute();
  }

  double
  get(code_t code, int i, int hamdist) const
  {
    return a[code][i][hamdist];
  }
  
private:

  void
  initialize()
  {
    int number_of_sequences = pow(4, k);
    for (code_t code=0; code < number_of_sequences; ++code) {
      if (f.get(code) > d)
	a[code][k][1] = compute_bernoulli_probability(code, k, q);
    }
  }

  void
  compute()
  {
    for (int i=k; i < imax; ++i) {
      if (i < (2*k-epsilon)) {
	int hampos = i-k+1;
	for (code_t code=0; code < number_of_sequences; ++code) {
	  if (f.get(code) > d) {
	    double p = a[code][i][hampos];
	    if (p > 0.0) {
	      for (int x=0; x < 4; ++x) {
		code_t newcode = ((code << 2) + x) & mask;
		int j = f.get(newcode) <= d ? 0 : hampos + 1;
		if (j <= k-epsilon)
		  a[newcode][i+1][j] += p * q[x];
	      }
	    }
	  }
	}
      }
      else {
	for (int hampos=0; hampos <= k-epsilon; ++hampos) {
	  for (code_t code=0; code < number_of_sequences; ++code) {
	    double p = a[code][i][hampos];
	    if (p > 0.0 and f.get(code) > e) {
	      for (int x=0; x < 4; ++x) {
		code_t newcode = ((code << 2) + x) & mask;
		int j = f.get(newcode) <= d ? 0 : hampos + 1;
		if (j <= k-epsilon)
		  a[newcode][i+1][j] += p * q[x];
	      }
	    }
	  }
	}
      }
    }
  }



  int d;
  int e;
  code_t s;
  std::vector<double> q;
  min_hamming_distance_class f;
  int k;
  int imax;
  int epsilon;
  int number_of_sequences;
  code_t mask;
  boost::multi_array<double, 3> a;
};

class cluster_probability_type
{
public:

  cluster_probability_type() {}
  
  cluster_probability_type(const std::string& seed, int d, int e, const std::vector<double>& q_, int epsilon_, int max_cluster_len,
			   const min_hamming_distance_class& f, const min_hamming_distance_class& f_rev)
  //    : k(seed.length()), q(q_), epsilon(epsilon_), lmax(max_cluster_len),
  //      prefix(seed, d, e, q, lmax, f), suffix(reverse(seed), d, e, q, lmax, f_rev)
  {
    init(seed, d, e, q_, epsilon_, max_cluster_len, f, f_rev);
  }

  void
  init(const std::string& seed, int d, int e, const std::vector<double>& q_, int epsilon_, int max_cluster_len,
			   const min_hamming_distance_class& f, const min_hamming_distance_class& f_rev)
  {
    k = seed.length();
    q = q_;
    epsilon = epsilon_;
    lmax = max_cluster_len;
    prefix.init(seed, d, e, q, lmax, f);
    suffix.init(reverse(seed), d, e, q, lmax, f_rev);

  }
  
  double
  operator()(const std::string& u, int cluster_len)
  {
    code_t code = dna_to_number(u);
    code_t code_rev = reverse_2bitstring(code, k);
    double p = 0.0;
    double div = compute_bernoulli_probability(code, k, q);
    for (int l1=2*k-epsilon; l1 <= cluster_len - (k-epsilon); ++l1) {
      p += prefix.get(code, l1, 0) * suffix.get(code_rev, cluster_len - l1 + k, 0) / div;
    }
    
    return p;
  }
  
private:
  int k;
  std::vector<double> q;
  int epsilon;
  int lmax;
  cluster_probability_array prefix;
  cluster_probability_array suffix;
};

class neighbourhood_probability_type
{
public:
  
  void
  init(const std::string& seed, int d, const std::vector<double>& q_, int epsilon_)
  {
    k = seed.length();
    q = q_;
    epsilon = epsilon_;
    l = 2*k - epsilon;
    prefix.resize(pow(4, l));
    suffix.resize(pow(4, l));
    result.resize(pow(4, k));
    f.init(seed);
    
    unsigned long long number_of_partial_sequences = pow(4, l);
    code_t maskk = ((code_t)1 << (2*k)) - 1;
    for (code_t code=0; code < number_of_partial_sequences; ++code) {
      double p = compute_bernoulli_probability(code, l, q);
      for (int shift=1; shift <= k - epsilon; ++shift) {
	if (f.get((code >> shift) & maskk) <= d) {
	  p = 0.0;
	  break;
	}
      }
      prefix[code] = p;
    }
    
    for (code_t code=0; code < number_of_partial_sequences; ++code) {
      double p = compute_bernoulli_probability(code, k-epsilon, q);  // Note! The center part is not included in the probability
      for (int shift=0; shift < k - epsilon; ++shift) {
	if (f.get((code >> shift) & maskk) <= d) {
	  p = 0.0;
	  break;
	}
      }
      suffix[code] = p;
    }
    
    code_t maskl = 1;
    maskl <<= (2*l);
    --maskl;
    unsigned long long number_of_full_sequences = pow(4, 3*k-2*epsilon);
    for (code_t code=0; code < number_of_full_sequences; ++code) {
      code_t first = code >> (2*(k-epsilon));
      code_t second = code & maskl;
      result[second >> (2*(k-epsilon))] += prefix[first] * suffix[second];
    }
  }

  double
  operator()(const std::string& u) const
  {
    code_t code = dna_to_number(u);
    return result[code];
  }
  
private:
  int k;
  std::vector<double> q;
  int epsilon;
  int l;
  std::vector<double> prefix;
  std::vector<double> suffix;
  std::vector<double> result;
  min_hamming_distance_class f;
};

boost::tuple<dmatrix,int>
find_multinomial_n_background(const std::string& seed, const std::vector<std::string>& sequences, const std::vector<double>& bg,
			      int n, bool use_multimer)
{
  TIME_START(t);
  const int k = seed.length();
  const int L = sequences[0].length();
  assert(n >= 0);
  assert(n <= k);
  dmatrix result(4, k);

  string_to_tuple_type string_to_tuple;
  string_to_tuple = get_n_neighbourhood(seed, n);

  //printf("Number of patterns %lu\n", string_to_tuple.size());
  //double seed_count=0;
  int epsilon = cluster_threshold;
  double total_count=0;
  int lines = sequences.size();
  int sites = lines * (L-(k+extra_flank*2)+1)*2;
  int neighbour_sites = lines * (L - (3*k-2*epsilon) + 1) * 2;
  iupac_probability_in_background iupac_prob(bg);
  //  int max_cluster_len = k + 4*(k-epsilon);  // This is crude approximation.
  int max_cluster_len = 200;  // This is crude approximation.
  int min_cluster_len = k + 2*(k-epsilon);  // By definition of cluster, the cluster length cannot be shorter than this
  //seed_count = sites * iupac_prob(seed);
  //total_count += seed_count;

  std::vector<boost::tuple<int, int> > pairs;
  std::string pattern;
  std::string seed_rev = reverse_complement(seed);
    
  std::vector<string_to_tuple_type> hamming_neighbours(n+1);   // Bin the patterns according to Hamming distance to the seed
  BOOST_FOREACH(boost::tie(pattern, pairs), string_to_tuple) {
    hamming_neighbours[hamming_distance(seed, pattern)].insert(std::make_pair(pattern, pairs));
  }
  double prob_sum = 0.0;
  min_hamming_distance_class f;
  min_hamming_distance_class f_rev;
  neighbourhood_probability_type neighbourhood_probability;
  if (background_counting == choose_one_per_cluster) {
    f.init(seed);
    f_rev.init(reverse(seed));
  }
  else if (background_counting == neighbourhood_contains_one) {
    neighbourhood_probability.init(seed, n, bg, epsilon);
  }
  for (int e=0; e <= n; ++e) {
    cluster_probability_type cluster_probability;
    if (background_counting == choose_one_per_cluster)
      cluster_probability.init(seed, n, e, bg, epsilon, max_cluster_len, f, f_rev);
    BOOST_FOREACH(boost::tie(pattern, pairs), hamming_neighbours[e]) {

      double p;
      double count=0;
      switch (background_counting) {
      case choose_one_per_cluster:
	for (int cluster_len=min_cluster_len; cluster_len <= max_cluster_len; ++cluster_len) {
	  int cluster_sites = lines * (L-cluster_len+1) * 2;
	  p = cluster_probability(pattern, cluster_len);
	  prob_sum += p;
	  count += cluster_sites * p;
	}
	break;
      case neighbourhood_contains_one:
	p = neighbourhood_probability(pattern);
	prob_sum += p;
	count = neighbour_sites * p;
	break;
      case all_occurrences:
	p = iupac_prob(pattern);
	prob_sum += p;
	count = sites * p;
	break;
      case sequence_contains_one:
	error(true, "Not implemented.");
	break;
      }
      //	if (not iupac_string_match(pattern, seed))
      total_count += count;
	
      //	double count = sites * (p*pow(0.25, extra_flank*2));
      //      if (use_palindromic_correction)
      //	count *= palindromic_correction(pattern, seed, seed_rev);
      
      int j, a;
      BOOST_FOREACH(boost::tie(j,a), pairs) {
	result(a, j) += count;
      }
    }
  }
    
  printf("Total probability of hamming neighbourhood in background is %f\n", prob_sum);
  TIME_PRINT("find_multinomial_n_background took %.2f seconds.\n", t);
  return boost::make_tuple(result, (int)total_count);
} // find_multinomial_n_background


boost::tuple<dmatrix, unsigned long, unsigned long>
count_all_occurrences(const string_to_tuple_type& string_to_tuple, const std::string& seed, const suffix_array& sa)
{
  unsigned long seed_count = 0;
  unsigned long total_count = 0;
  int k = seed.length();
  
  dmatrix result(4, k);
  std::vector<boost::tuple<int, int> > pairs;
  std::string pattern;
  if (print_alignment) 
    printf("#String\tColumn\tHamming distance\tPalindrome\tCount\tMatches at col\n");
  BOOST_FOREACH(boost::tie(pattern, pairs), string_to_tuple) {
    unsigned long count = sa.count_iupac(pattern);
    bool is_palindrome = is_palindromic(pattern);
    int hd = hamming_distance(seed, pattern);
    if (is_palindrome and use_two_strands and not count_palindromes_twice)
      count /= 2;
    //      if (use_palindromic_correction)
    //	count *= palindromic_correction(pattern, seed, seed_rev);

    if (iupac_string_match(pattern, seed))
      seed_count += count;
    total_count += count;
    int j, a;
    BOOST_FOREACH(boost::tie(j,a), pairs) {
      if (print_alignment) 
	printf("#%s\t%i\t%i\t%s\t%zu\t%s\n", pattern.c_str(), j, hd,
	       yesno(is_palindrome), count, yesno(seed[j]==pattern[j]));
      result(a, j) += count;
    }
  }

  return boost::make_tuple(result, seed_count, total_count);
};


// Note! This assumes that all sequences are of equal length, for efficiency
boost::tuple<dmatrix, unsigned long, unsigned long>
count_sequence_contains_one(const string_to_tuple_type& string_to_tuple, const std::string& seed, const suffix_array& sa,
			    const std::vector<std::string>& sequences)
{
  // allow reads that contain exactly one hit.
  int lines = sequences.size();

  int L = sequences[0].length();
  for (int i=0; i < lines; ++i)
    assert(sequences[i].length() == L);
  
  unsigned long seed_count = 0;
  unsigned long total_count = 0;
  int k = seed.length();
  
  dmatrix result(4, k);
  typedef string_to_tuple_type::const_iterator iterator;
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
  
		    
  return boost::make_tuple(result, seed_count, total_count);
};

boost::tuple<dmatrix, unsigned long, unsigned long>
count_choose_one_per_cluster(const string_to_tuple_type& string_to_tuple, const std::string& seed, const suffix_array& sa,
			     const std::vector<std::string>& sequences)
{
  typedef string_to_tuple_type::const_iterator iterator;
  int lines = sequences.size();
  int L = sequences[0].length();
  for (int i=0; i < lines; ++i)
    assert(sequences[i].length() == L);

  unsigned long seed_count = 0;
  unsigned long total_count = 0;
  int k = seed.length();
  
  dmatrix result(4, k);
  std::vector<std::vector<int> > hit_positions(lines);
  std::vector<std::vector<iterator> > hit_patterns(lines);
  int divisor = L + 1;    // Includes the separator '#'
  std::vector<boost::tuple<int, int> > pairs;
  std::string pattern;

  // Bin the occurrences according to the sequence they are in.
  for (iterator it=string_to_tuple.begin(); it != string_to_tuple.end(); ++it) {  // iterator through string in Hamming neighbourhood
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
      if (not is_palindrome or count_palindromes_twice) {
	hit_positions[i].push_back(j);
	hit_patterns[i].push_back(it);
      }
      else {
	if (std::find(hit_positions[i].begin(), hit_positions[i].end(), j) == hit_positions[i].end()) {
	  hit_positions[i].push_back(j);
	  hit_patterns[i].push_back(it);
	}
      }
    }
  }
  
  for (int i=0; i < lines; ++i) {
    int hit_count = hit_positions[i].size();
    std::vector<int> best_from_each_cluster;
    if (hit_count <= 1)
      best_from_each_cluster.push_back(0);
    else {

      ////////////////////////
      // Form the clusters

      int max_cluster_size=0;
      int max_cluster_length=0;
      std::vector<std::vector<int> > clusters;
      std::vector<int> hits;               // This will contain indices to hit_positions/patterns vector
      for (int index=0; index < hit_count; ++index)
	hits.push_back(index);
      std::sort(hits.begin(), hits.end(), [i, &hit_positions](int a, int b) { return hit_positions[i][a] < hit_positions[i][b];} );
      std::vector<int> new_cluster;
      int previous_end_position=sequences[i].length();  // Sure to give an overlap in the following test
      //int distance_threshold = -1000;  // each occurrence will be considered its own cluster
      for (int t=0; t < hit_count; ++t) {
	if (previous_end_position - hit_positions[i][hits[t]] + 1 >= cluster_threshold) {   // do two consequent hits overlap?
	  new_cluster.push_back(hits[t]);           // extend the old cluster
	  previous_end_position = hit_positions[i][hits[t]] + k - 1;
	} else {
	  clusters.push_back(new_cluster);   // store the old cluster
	  new_cluster.clear();               // Start a new cluster
	  new_cluster.push_back(hits[t]);    // new cluster now contains the current hit
	  previous_end_position = hit_positions[i][hits[t]] + k - 1;
	}
      }
      clusters.push_back(new_cluster);  // last cluster

      std::vector<int> cluster;         // cluster contains indices to hit_patterns vector
      BOOST_FOREACH(cluster, clusters) {
	int size = cluster.size();
	if (size > max_cluster_size)
	  max_cluster_size = size;
	int len = hit_positions[i][cluster[size-1]] + k - hit_positions[i][cluster[0]];
	if (len > max_cluster_length)
	  max_cluster_length = len;
      }
	  
      ////////////////////////
      // Iterate the clusters

      std::vector<int> cluster_size_distribution(max_cluster_size+1);
      std::vector<int> cluster_length_distribution(max_cluster_length+1);
      BOOST_FOREACH(cluster, clusters) {
	int size = cluster.size();
	int len = hit_positions[i][cluster[size-1]] + k - hit_positions[i][cluster[0]];
	++cluster_size_distribution[cluster.size()];
	++cluster_length_distribution[len];
	if (cluster.size() == 1)
	  best_from_each_cluster.push_back(cluster[0]);
	else {
	  std::map<int, int> hds;  // map from index to Hamming distance
	  BOOST_FOREACH(int index, cluster)
	    hds[index] = hamming_distance(seed, hit_patterns[i][index]->first);

	  // sort cluster by Hamming distance of occurrences to the seed
	  std::sort(cluster.begin(), cluster.end(), [i, &hit_positions, &hds](int a, int b) { 
	      // a and b are indices to hit_patterns vector

	      // Sort primarily by Hamming index and secondarily by start position
	      return hds[a] < hds[b] or (hds[a] == hds[b] and hit_positions[i][a] < hit_positions[i][b]);
	    } );
	    
	  int best_hd = hds[cluster[0]];
	  int best_hd_count = 0;

	  for (int t = 0; t < cluster.size() and hds[cluster[t]] == best_hd; ++t)
	    ++best_hd_count;

	  if (best_hd_count <= 2)
	    best_from_each_cluster.push_back(cluster[0]); // Add the element (index) of cluster with smallest Hamming distance
	  else {
	    int best_theoretical_starting_point = floor((hit_positions[i][cluster[best_hd_count-1]] - hit_positions[i][cluster[0]])/2); // center point
	    int min_dist=std::numeric_limits<int>::max();
	    int min_arg=0;
	    for (int t=0; t < best_hd_count; ++t) {
	      int dist = abs(best_theoretical_starting_point - hit_positions[i][cluster[t]]);
	      if (dist < min_dist) {
		min_dist = dist;
		min_arg = t;
	      }
	    }
	    best_from_each_cluster.push_back(cluster[min_arg]); // Add the element of cluster with smallest Hamming distance, if not unique
	    // use the centermost of those having the best Hamming distance
	  }
	      
	}
      }  // end BOOST_FOREACH cluster
      printf("#Cluster size\tCount\n");
      for (int i=0; i<= max_cluster_size; ++i) {
	printf("#%i\t%i\n", i, cluster_size_distribution[i]);
      }
	  
      printf("$Cluster length\tCount\n");
      for (int i=0; i<= max_cluster_length; ++i) {
	printf("$%i\t%i\n", i, cluster_length_distribution[i]);
      }
    }
    for (int t=0; t < best_from_each_cluster.size(); ++t) {
      boost::tie(pattern, pairs) = *(hit_patterns[i][best_from_each_cluster[t]]);
      if (not iupac_string_match(pattern, seed))
	total_count += 1;
      else
	++seed_count;
      int j, a;
      BOOST_FOREACH(boost::tie(j,a), pairs) {
	result(a, j) += 1;
      }
    }
  }  // end i
      
  return boost::make_tuple(result, seed_count, total_count);
};


boost::tuple<dmatrix, unsigned long, unsigned long>
count_neighbourhood_contains_one(const string_to_tuple_type& string_to_tuple, const std::string& seed, const suffix_array& sa,
				  const std::vector<std::string>& sequences)
{
  typedef string_to_tuple_type::const_iterator iterator;
  int lines = sequences.size();
  int L = sequences[0].length();
  for (int i=0; i < lines; ++i)
    assert(sequences[i].length() == L);

  unsigned long seed_count = 0;
  unsigned long total_count = 0;
  int k = seed.length();
  
  dmatrix result(4, k);
  std::vector<std::vector<int> > hit_positions(lines);
  std::vector<std::vector<iterator> > hit_patterns(lines);
  int divisor = L + 1;    // Includes the separator '#'
  std::vector<boost::tuple<int, int> > pairs;
  std::string pattern;

  // Bin the occurrences according to the sequence they are in.
  for (iterator it=string_to_tuple.begin(); it != string_to_tuple.end(); ++it) {  // iterator through string in Hamming neighbourhood
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
      if (not is_palindrome or count_palindromes_twice) {
	hit_positions[i].push_back(j);
	hit_patterns[i].push_back(it);
      }
      else {
	if (std::find(hit_positions[i].begin(), hit_positions[i].end(), j) == hit_positions[i].end()) {
	  hit_positions[i].push_back(j);
	  hit_patterns[i].push_back(it);
	}
      }
    }
  }
  
  for (int i=0; i < lines; ++i) {
    int hit_count = hit_positions[i].size();
    std::vector<int> non_intersecting_occurences;
    if (hit_count <= 1)
      non_intersecting_occurences.push_back(0);
    else {

      ////////////////////////
      // Form the clusters, and store clusters with size one to non_intersecting_occurences container

      std::vector<std::vector<int> > clusters;
      std::vector<int> hits;               // This will contain indices to hit_positions/patterns vector
      for (int index=0; index < hit_count; ++index)
	hits.push_back(index);
      std::sort(hits.begin(), hits.end(), [i, &hit_positions](int a, int b) { return hit_positions[i][a] < hit_positions[i][b];} );
      std::vector<int> new_cluster;
      int previous_end_position=sequences[i].length();  // Sure to give an overlap in the following test
      //int distance_threshold = -1000;  // each occurrence will be considered its own cluster
      for (int t=0; t < hit_count; ++t) {
	if (previous_end_position - hit_positions[i][hits[t]] + 1 >= cluster_threshold) {   // do two consequent hits overlap?
	  new_cluster.push_back(hits[t]);           // extend the old cluster
	  previous_end_position = hit_positions[i][hits[t]] + k - 1;
	} else {
	  if (new_cluster.size() == 1)
	    clusters.push_back(new_cluster);   // store the old cluster
	  new_cluster.clear();               // Start a new cluster
	  new_cluster.push_back(hits[t]);    // new cluster now contains the current hit
	  previous_end_position = hit_positions[i][hits[t]] + k - 1;
	}
      }
      if (new_cluster.size() == 1)
	clusters.push_back(new_cluster);  // last cluster



	  
      ////////////////////////
      // Iterate the clusters

      std::vector<int> cluster;         // cluster contains indices to hit_patterns vector
      BOOST_FOREACH(cluster, clusters) {
	non_intersecting_occurences.push_back(cluster[0]);
      }  // end BOOST_FOREACH cluster
 
    }

    for (int t=0; t < non_intersecting_occurences.size(); ++t) {
      boost::tie(pattern, pairs) = *(hit_patterns[i][non_intersecting_occurences[t]]);
      if (iupac_string_match(pattern, seed))
	++seed_count;
      
      ++total_count;
      int j, a;
      BOOST_FOREACH(boost::tie(j,a), pairs) {
	result(a, j) += 1;
      }
    }
  }  // end i
      
  return boost::make_tuple(result, seed_count, total_count);
};

boost::tuple<dmatrix,int>
find_multinomial_n_suffix_array(const std::string& seed, const std::vector<std::string>& sequences, const suffix_array& sa, int n, bool use_multimer)
{
  const int k = seed.length();
  //  const int L = sequences[0].length();
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
  switch (data_counting) {
  case all_occurrences:
    boost::tie(result, seed_count, total_count) =
      count_all_occurrences(string_to_tuple, seed, sa);
    break;
  case choose_one_per_cluster:  
    boost::tie(result, seed_count, total_count) = count_choose_one_per_cluster(string_to_tuple, seed, sa, sequences);
    break;
  case neighbourhood_contains_one:
    boost::tie(result, seed_count, total_count) = count_neighbourhood_contains_one(string_to_tuple, seed, sa, sequences);
    break;
  case sequence_contains_one:
    boost::tie(result, seed_count, total_count) = count_sequence_contains_one(string_to_tuple, seed, sa, sequences);
    break;
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
