#include "dinucleotide.hpp"
#include "multinomial_helper.hpp"
#include "iupac.hpp"
#include "common.hpp"
#include "matrix.hpp"
#include "parameters.hpp"
#include "matrix_tools.hpp"
#include "kmer_tools.hpp"
#ifndef TIMING
#define TIMING 1
#endif
#include "timing.hpp"
#include "data.hpp"
#include "probabilities.hpp"
#include "my_assert.hpp"
//#include "aho_corasick_wrapper.hpp"

#include <cstring>
#include <set>
#include <fstream>
#include <boost/tuple/tuple.hpp>
#include <boost/unordered_map.hpp>
#include <boost/foreach.hpp>

//typedef boost::unordered_map<int, int> id_to_count_type;
//id_to_count_type id_to_count;

//typedef boost::unordered_map<big_int, std::vector<boost::tuple<int, int, int, int> > > id_to_tuple_type;
//id_to_tuple_type id_to_tuple;

//typedef boost::unordered_map<std::string, std::vector<boost::tuple<int, int> > > string_to_tuple_type;

//extern bool use_multimer;
bool use_grouped_dinucleotide_method = false;

extern bool require_directional_seed;

//static char nucs[] = "ACGT";

std::vector<std::string> dinuc_headers =
  {"AA","AC","AG","AT",
   "CA","CC","CG","CT",
   "GA","GC","GG","GT",
   "TA","TC","TG","TT"};

const char* mono_headers[] = {"ADM_MONO_A", "ADM_MONO_C", "ADM_MONO_G", "ADM_MONO_T"};


/*
void
write_dinuc_model(FILE* fp, const dinuc_model<double>& adm, const std::string& tag, 
		  const std::string& format)
{
  fprintf(fp, "%s", tag.c_str());

  adm.write(fp, format);
}
*/

bool
almost_equal(const dmatrix& dm1, const dmatrix& dm2, double threshold)
{
  return distance(dm1, dm2) < threshold;
}

bool
almost_equal(double d1, double d2, double threshold)
{
  return fabs(d1-d2) < threshold;
}


bool
operator==(const dinuc_model<double>& dm1, const dinuc_model<double>& dm2)
{
  return dm1.k == dm2.k and almost_equal(dm1.ip, dm2.ip) and almost_equal(dm1.dm, dm2.dm);
}






dinuc_model<double>
pwm_to_dinucleotide(const dmatrix& pwm)
{
  int rows, cols;
  boost::tie(rows, cols) = pwm.dim();
  assert(rows == 4);

  dmatrix result(16, cols);

  for (int j=0; j < cols-1; ++j) {
    for (int c = 0; c < 16; ++c) {
      result(c, j+1) = pwm(c>>2,j) * pwm(c&3,j+1);
    }
  }
  for (int a=0; a < 4; ++a)
    result(a, 0) = pwm(a, 0);
  
  dinuc_model<double> dinuc(result);
  
  return dinuc;
}





 



// this computes the dinucleotide-n matrix counts by scanning all possible windows
std::vector<dmatrix>
dinucleotide_counts_scan(const std::string& seed, const std::vector<std::string>& sequences, int n)
{
  TIME_START(t);
  const int k = seed.length();
  const int L = sequences[0].length();
  assert(n >= 0);
  assert(n <= k);
  //  char nucs[] = "ACGT";
  std::vector<dmatrix> result;
  result.push_back(dmatrix(16, k-1));
  for (int i=0; i < sequences.size(); ++i) {
    int max_dir = use_two_strands ? 2 : 1;
    for (int dir=0; dir < max_dir; ++dir) {
      const std::string& line = dir == 0 ? sequences[i] : reverse_complement(sequences[i]);
      for (int j=0; j < L-k+1; ++j) {
	int hd = iupac_hamming_dist(line.substr(j, k), seed, n);
	if (hd > n)
	  continue;
	for (int pos=0; pos < k-1; ++pos) {
	  char c1 = line[j+pos];
	  char c2 = line[j+pos+1];
	  int cc = iupac_match(c1, seed[pos]) ? 0 : 1;
	  cc += iupac_match(c2, seed[pos+1]) ? 0 : 1;
	  if (hd-cc <= n-2)
	    ++result[0]((to_int(c1)<<2) + to_int(c2), pos);
	}
      }
    }
  }

  TIME_PRINT("Dinucleotide-n scanning algorithm took %.2f seconds\n", t);
  return result;
}




void
print_ics(const dinuc_model<double>& adm, std::vector<double> bg = std::vector<double>(4, 0.25))
{
  int k = adm.get_length();

  std::vector<double> ic(k);
  for (int i=0; i < k; ++i) 
    ic[i] = information_content_of_adm_column(adm, i, bg);
  printf("Information content by columns\n");
  printf("%s\n", print_vector(ic, "\t", 2).c_str());
  printf("Average information content is %.2f\n", sum(ic) / k);
}







/*
double
dinuc_model::bias_string_probability(int prev_nucleotide, int start_pos,
				     const std::string& bias_string, int bias_string_start_pos)
{
  if (bias_string_start_pos == bias_string.length())
    return 1.0;
  double result = 0.0;
  std::string c = iupac_class(bias_string[bias_string_start_pos]);
  for (int j=0; j < c.length(); ++j) {
    int a = to_int(c[j]);
    result += cp(prev_nucleotide*4+a, start_pos-1)*bias_string_probability(a, start_pos+1, bias_string, bias_string_start_pos+1);
  }
  return result;
}
*/





/*
template<typename T>
void
dinuc_model<T>::add_pseudo_count(double p)
{
  assert(p < 1.0);
  dm = p + dm;  // add pseudo counts
  dm = conditional_probabilities(dm);
  ip = dmatrix(4,k);
  // transforms dinucleotide counts in the dmatrix(16,k) to 
  // initial_probabilities dmatrix(4,k)

  generate_initial_probabilities(dm, ip);
}
*/

/*
double
dinuc_model::bias_helper(int prev_nucleotide_orig, int start_pos, const std::string& seed, int d)
{
  static character_to_values<bool> isnuc("ACGT", true);

  int len=k-start_pos;
  std::string bias_string = seed.substr(start_pos, len);
  std::vector<std::string> neighbourhood = get_n_neighbourhood(bias_string, d);
  double sum=0.0;
  BOOST_FOREACH(std::string s, neighbourhood) {
    //    sum += bias_string_probability(prev_nucleotide, start_pos, s, 0);
    double temp = 1.0;
    int prev_nucleotide = prev_nucleotide_orig;
    for (int i=0; i < s.length(); ++i) {
      temp *= cp(prev_nucleotide*4+to_int(s[i]), i+start_pos);
      prev_nucleotide = to_int(s[i]);
    }
    sum += temp;
  }

  return sum;
}
*/

/*

void
dinuc_model::correct_for_seed_bias(const std::string& seed, int hd)
{
  assert(seed.length() == k);

  int count = 0;
  dvector temp(4, 0.0);
  for (int i=k-2; i >= 0; --i) {
    for (int a=0; a<4; ++a) {            // character in pos i-1
      if (i == 0 and a != 0)             // Initial symbol in pos -1 is 'A'
	break;
      bool a_matches_seed = i > 0 ? iupac_match(nucs[a], seed[i-1]) : true;  // handles the implisit 'A' in pos -1
      int d = std::min(hd-1-(a_matches_seed ? 0 : 1), k - (i+1));
      for (int b=0; b<4; ++b) {          // character in pos i
	double divisor = bias_helper(b, i+1, seed, d);
	if (divisor > 0.0)
	  temp[b] = cp(a*4+b, i) / divisor;  // transition from i-1 to i
	else {
	  temp[b] = 0.0;
	  ++count;
	}
      }
      if (sum(temp) > 0.0)
	normalize_vector(temp);
      else
	++count;
      for (int b=0; b<4; ++b)
	cp(4*a+b, i) = temp[b];
    }
  }

  generate_initial_probabilities(cp, ip);
  if (count != 0)
    printf("In correct_for_seed_bias: division by zero %i times\n", count);
}
*/


/*
template<typename T>
dinuc_model<T>::dinuc_model(const dmatrix& dm_, bool normalize)
{
  init(dm_, normalize);
}
*/



/*
template<typename T>
void
dinuc_model<T>::print(const std::string& format) const
{
  this->write(stdout, format);

}
*/

/*
template<typename T>
void
dinuc_model<T>::write(const std::string& filename, const std::string& format) const
{
  FILE* fp;
  fp = fopen(filename.c_str(), "w");
  error(fp == NULL, to_string("Could not open file %s\n", filename.c_str()));

  dinuc_model::write("", format, fp);
  fclose(fp);
}
*/


std::ostream&
operator<<(std::ostream& str, const dinuc_model<double>& adm)
{

  int k = adm.get_length();
  bool col_headers = false;
  str.precision(6);
  str.setf(std::ios_base::fixed, std::ios_base::floatfield);
  //  str.width(8);
  str << std::endl;
  if (col_headers) {
    for (int i=0; i < k-1; ++i) {
      str << "\t" << i;
      //fprintf(fp, "\t%i", i);
    }
    //fprintf(fp, "\n");
    str << std::endl;
  }
   
  for (int c=0; c < 16; ++c) {
    for (int i=0; i < k-1; ++i) {
      str << adm.dm(c, i+1);
      //      fprintf(fp, format.c_str(), dm(c, i+1));
      //      fprintf(fp, "\t");
      str << "\t";
    }
    str << "ADM_DI\t" << dinuc_headers[c] << std::endl;
    //    fprintf(fp, "%s\t%s\n", "ADM_DI", dinuc_headers[c]);
  }
  for (int c=0; c < 4; ++c) {
    for (int i=0; i < k; ++i) {
      //      fprintf(fp, format.c_str(), ip(c, i));
      //      fprintf(fp, "\t");
      str << adm.ip(c, i);
      str << "\t";
    }
    //    fprintf(fp, "%s\n", mono_headers[c]);
    str << mono_headers[c] << std::endl;
  }

  return str;
}	   





/*

string_to_tuple_type
get_n_neighbourhood_dinucleotide_contributions(const std::string&seed, int n)
{
  const int k = seed.length();
  assert(n >= 0);
  assert(n <= k);

  string_to_tuple_type string_to_contributions;

  std::vector<std::string> complements(k);  // These are set complements, not nucleotide complements
  std::vector<int> bases(k);
  unsigned long long N_mask=0;                // bitmask for positions that contain 'N'. Those positions cannot contain an error
  for (int i=0; i < k; ++i) {
    complements[i] = complement_set(seed[i]);
    bases[i]=complements[i].length() - 1;
    if (seed[i]=='N')
      N_mask |=  (1 << (k-1-i));
  }
  for (int j=0; j < k-1; ++j) {        // iterate through all possible dinucleotide positions
    std::string temp = seed;
    for (int a=0; a < 16; ++a) {      // number of errors outside positions j and j+1 is zero, handles Hamming distances 0,1,2
      int number_of_mismatches = not iupac_match(nucs[a/4], seed[j]) + not iupac_match(nucs[a%4], seed[j+1]);
      assert(number_of_mismatches <= 2);  // Just to be sure that bools are converted to either 0 or 1
      if (number_of_mismatches > n)
	continue;
      temp[j]   = nucs[a/4];
      temp[j+1] = nucs[a%4];
      string_to_contributions[temp].push_back(boost::make_tuple(j, a));  // ei mene oikein, tai menee sittenkin
    }

    for (int error=1; error <= n-2; ++error) { // errors outside positions j and j+1, handles hamming distances 1 <= hd <= n
      // bitvector c has 1-bit for each member of the subset, rightmost bit is bit number k-1
      unsigned long long c = (1ull<<error)-1;  // initially rightmost 'error' bits are 1
      //int mycount = 0;
      // iterate through all subsets c of {0, ..., k-1} that have size 'error'
      while (c < (1ull<<k)) {   // Superset has only k elements
	assert(__builtin_popcountll(c) == error);
	if (((c & (1ull << (k-1-j)))) == 0 and     // pos j not included in c
	    ((c & (1ull << (k-1-(j+1))))) == 0 and // pos j+1 not included in c
	    ((c & N_mask) == 0))  { // j and j+1 don't belong to the subset, and subset positions don't contain 'N'
	  //++mycount;
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
	  
	  std::vector<int> y(error, 0);
	  y[error-1]=-1;  // Initialize
	  for (int r=0; r < number_of_combinations; ++r) {
      
	    int i;
	    for (i=error-1; y[i] == bases[P[i]]; --i) {
	      y[i]=0;
	      temp[P[i]] = complements[P[i]][y[i]];
	    }
	    y[i]++;
	    temp[P[i]] = complements[P[i]][y[i]];

	    for (int a=0; a < 16; ++a){
	      if (n == 1 and not iupac_match(nucs[a/4], seed[j]) and not iupac_match(nucs[a%4], seed[j+1]) )
		continue;
	      temp[j]   = nucs[a/4];
	      temp[j+1] = nucs[a%4];
	      string_to_contributions[temp].push_back(boost::make_tuple(j, a));
	    }
	  }  // end for r

	}    // end if j not in c

	unsigned long long a = c&-c;
	unsigned long long b = c+a;   // update bitvector c. This is "Gosper's hack"
	c = (c^b)/4/a|b;

      } // end foreach subset c


    }  // end for error
  } // end for j

  
  // Just a sanity check for the case when the seed does not contain iupac characters
  // if (is_nucleotide_string(seed)) {
  //   std::string pattern;
  //   std::vector<boost::tuple<int, int> > pairs;
    
  //   BOOST_FOREACH(boost::tie(pattern, pairs), string_to_contributions) {
  //     assert(pairs.size() == k or pairs.size() == n);
  //   }
  // }
  
  return string_to_contributions;
}

*/

  
//typedef string_to_tuple_type::const_iterator iterator;
typedef boost::tuple<int, int, std::string> mytuple;

boost::tuple<std::vector<dmatrix>, unsigned long, unsigned long>
helper(std::vector<std::vector<mytuple > >& hit_info, int k, const std::string& seed, int hamming_radius,
       bool use_grouped_dinucleotide_method_local = use_grouped_dinucleotide_method)
{ 
  int lines = hit_info.size();
  int rows = 16;
  std::vector<dmatrix> result;
  if (use_grouped_dinucleotide_method_local)
    for (int i=0; i < hamming_radius; ++i)
      result.push_back(dmatrix(rows, k));
  else
    result.push_back(dmatrix(rows, k));
    
  unsigned long seed_count = 0;
  unsigned long total_count = 0;
  std::vector<boost::tuple<int, int> > pairs;
  std::string neighbour;
  const int mask = 15;

    
  if (print_alignment) 
    printf("#String\tColumn\tHamming distance\tPalindrome\tCount\tMatches at col\n");
  for (int i=0; i < lines; ++i) {
 
    std::sort(hit_info[i].begin(), hit_info[i].end(),
	      [](mytuple& a, mytuple& b) { return a.get<0>() < b.get<0>();} );
    int start_pos;
    int dir;

    BOOST_FOREACH(boost::tie(start_pos, dir, neighbour), hit_info[i]) {
      //printf("%s\n", neighbour.c_str());
      if (print_alignment) {
	std::string pal = is_palindromic(neighbour) ? "Palindrome" : "-";
	int hd = hamming_distance(seed, neighbour);
	printf("*%i\t%i\t%i\t%i\t%s\t%s\t%i\n", i, start_pos, dir, k, neighbour.c_str(), pal.c_str(), hd);
      }

      if (iupac_string_match(neighbour, seed))
	seed_count += 1;
      total_count += 1;
	
      std::vector<int> positions = iupac_hamming_mismatches(neighbour, seed);
      int hd = positions.size();

      // j refers to the last character of a dinucleotide !!!!!!!!!!!!!!!!!
	
      code_t code = dna_to_number<code_t>(neighbour);
      int r = 0; // Number of mismatches before position j
      if (hd <= hamming_radius - 1) {
	for (int j=0; j < k; ++j) {      // initial probibilities will be in result(.,0), trans. prop are for j>0
	  int a = (code >> ((k-j-1)*2)) & mask;              // get dinucleotides
	  result[r](a, j) += 1;
	  if (use_grouped_dinucleotide_method_local and not iupac_match(neighbour[j], seed[j]))
	    ++r;
	}
      }
      else {  // hd == hamming_radius
	//	  BOOST_FOREACH(int j, positions) {                  // only use mismatch positions
	for (int j=0; j < k; ++j) {
	  if (not iupac_match(neighbour[j], seed[j])) {
	    int a = (code >> ((k-j-1)*2)) & mask;              // get dinucleotides
	    result[r](a, j) += 1;
	    if (use_grouped_dinucleotide_method_local)
	      ++r;
	  }
	}
      }
    
	
    } // end for index
  } // end i

  return boost::make_tuple(result, seed_count, total_count);
}
    
boost::tuple<std::vector<dmatrix>, unsigned long, unsigned long>
count_all_occurrences(const std::vector<std::string>& neighbourhood, const std::string& seed, const suffix_array& sa,
		      const std::vector<std::string>& sequences, int hamming_radius)
{
  int k = seed.length();
  int lines = sequences.size();
  int L = sequences[0].length();


  int rows = 16;
  //  const int mask = 15;
  std::vector<dmatrix> result;
  if (use_grouped_dinucleotide_method)
    for (int i=0; i < hamming_radius; ++i)
      result.push_back(dmatrix(rows, k));
  else
    result.push_back(dmatrix(rows, k));
    
  //unsigned long seed_count = 0;
  //unsigned long total_count = 0;
  std::string neighbour;


  //int sites = lines * (L-k+1);
  //double lambda = 0.0;
  //std::vector<double> bg(4, 0.25);
  //iupac_probability_in_background iupac_prob(bg);

    
  std::vector<std::vector<mytuple > > hit_info(lines);
  
  int divisor = L + 1;    // Includes the separator '#'
  
  std::vector<boost::tuple<int, int> > pairs;
  std::string pattern;

  // All the following hassle is just to categorize the hits by the sequence they appear in,
  // and to make sure that palindromes are counted correctly.
  BOOST_FOREACH(std::string neighbour, neighbourhood) {

    
    std::vector<long int> positions;
    sa.locate_iupac(neighbour, positions);
    BOOST_FOREACH(int pos, positions) {
      int i = pos / divisor;  // index of the read containing the pos
      int j = pos % divisor;
      int dir = 1;
      if (i >= lines) {       // handle reverse complement
	i = 2*lines - i - 1;
	j = L - (j + k - 1) - 1;
	dir = -1;
      }
      bool is_palindrome = is_palindromic(sequences[i].substr(j, k));
      if (not is_palindrome or count_palindromes_twice) { // This branch allows the same site to be counted twice, for both orientations
	hit_info[i].push_back(boost::make_tuple(j, dir, neighbour));
      }
      else {  // count palindrome only once
	if (std::find_if(hit_info[i].begin(), hit_info[i].end(),
			 [j](mytuple& x) {return x.get<0>()==j;}) == hit_info[i].end()) {
	  hit_info[i].push_back(boost::make_tuple(j, dir, neighbour));
	}
      }
    } // end foreach pos

 
	
  } // end for neighbour

   

  
  return helper(hit_info, k, seed, hamming_radius);
};




std::vector<dmatrix>
dinucleotide_counts_suffix_array(const std::string& seed, const std::vector<std::string>& sequences,
				 const suffix_array& sa, int hamming_radius)
{
  TIME_START(t);
  unsigned long seed_count = 0;
  unsigned long total_count = 0;
  const int k = seed.length();
  assert(hamming_radius >= 0);
  assert(hamming_radius <= k);
  //  std::string str1;
  //  std::string str2;

    
  std::vector<std::string> neighbourhood = get_n_neighbourhood(seed, hamming_radius);
  //  BOOST_FOREACH(std::string s, neighbourhood)
  //    printf("%s\n", s.c_str());
  std::vector<dmatrix> result;
  //  suffix_array sa(str1);
  boost::tie(result, seed_count, total_count) =
    count_all_occurrences(neighbourhood, seed, sa,
			  sequences, hamming_radius);
  TIME_PRINT("Dinucleotide-n algorithm took %.2f seconds\n", t);
  printf("Seed %s count = %lu\n", seed.c_str(), seed_count);
  printf("Total dinucleotide-n count is %lu\n", total_count);

  return result;
} // dinucleotide_counts_suffix_array


boost::tuple<std::vector<dmatrix>, double>
find_dinucleotide_n_background(const std::string& seed, const std::vector<std::string>& sequences, const std::vector<double>& bg,
			       int hamming_radius,
			       bool local_use_grouped_dinucleotide_method)
{

  std::vector<std::string> neighbourhood = get_n_neighbourhood(seed, hamming_radius);
  int k = seed.length();
  int rows = 16;
  std::vector<dmatrix> result;
  if (local_use_grouped_dinucleotide_method)
    for (int i=0; i < hamming_radius; ++i)
      result.push_back(dmatrix(rows, k));
  else
    result.push_back(dmatrix(rows, k));
  
  unsigned long seed_count = 0;
  unsigned long total_count = 0;
  std::vector<boost::tuple<int, int> > pairs;
  std::string pattern;
  double prob_sum = 0.0;
  iupac_probability_in_background iupac_prob(bg);   // Normal probability computation should be enough since
                                                    // since the iupac string were expanded earlier
  
  const int mask = 15;

  BOOST_FOREACH(std::string neighbour, neighbourhood) {
    double p = 0.0;
    switch (background_counting) {
    case all_occurrences:
      p = iupac_prob(neighbour);
      if (use_two_strands)
	p += iupac_prob(reverse_complement(neighbour));
      prob_sum += p;
    case sequence_contains_one:
    case neighbourhood_contains_one:
    case choose_one_per_cluster:
      error(true, "Not implemented");
    }
      
    if (iupac_string_match(neighbour, seed))
      seed_count += p;
    total_count += p;

    std::vector<int> positions = iupac_hamming_mismatches(neighbour, seed);
    int hd = positions.size();
    
    code_t code = dna_to_number<code_t>(neighbour);
    int r = 0;
    if (hd <= hamming_radius - 1) {
      int shift = (k-1)*2;
      for (int j=0; j < k; ++j, shift -= 2) {      // initial probibilities will be in result(.,0), trans. prop are for j>0
	int a = (code >> shift) & mask;              // get dinucleotides
	result[r](a, j) += p;
	if (local_use_grouped_dinucleotide_method and not iupac_match(neighbour[j], seed[j]))
	  ++r;
      }
    }
    else {  // hd == hamming_radius
      //BOOST_FOREACH(int j, positions) {                  // only use mismatch positions
      for (int j=0; j < k; ++j) {
	if (not iupac_match(neighbour[j], seed[j])) {
	  int a = (code >> ((k-j-1)*2)) & mask;              // get dinucleotides
	  result[r](a, j) += p;
	  if (local_use_grouped_dinucleotide_method)
	    ++r;
	}
      }
    }
    
  } // end for neighbour

  return boost::make_tuple(result, prob_sum);
}

/*  
dmatrix
align_all_dinucleotide(const std::vector<std::string>& sequences)
{
  int L=sequences[0].length();
  int lines = sequences.size();
  int k = L;
  dmatrix result(16, k);
  for (int i=0; i < lines; ++i) {
    const std::string& line = sequences[i];
    code_t code = dna_to_number<code_t>(line);
    for (int j=k-1; j >= 0; --j, code >>= 2) 
      result(code & 15, j) += 1;
  }

  return result;
}
*/

/*
boost::tuple<std::vector<dmatrix>, unsigned long, unsigned long>
adm_count_neighbourhood_contains_one(const std::string& seed, int hamming_radius, const suffix_array& sa,
				     const std::vector<std::string>& sequences,
				     bool use_grouped_dinucleotide_method_local)
{
  std::vector<std::string> neighbourhood = get_n_neighbourhood(seed, hamming_radius);

  //  typedef string_to_tuple_type::const_iterator iterator;
  int lines = sequences.size();
  int L = sequences[0].length();
  for (int i=0; i < lines; ++i)
    assert(sequences[i].length() == L);

  //unsigned long seed_count = 0;
  //unsigned long total_count = 0;
  int k = seed.length();
  std::string seed_rev = reverse_complement(seed);

  std::vector<std::vector<int> > hit_positions(lines);
  std::vector<std::vector<int> > hit_directions(lines);
  std::vector<std::vector<std::string> > hit_contributions(lines);
  std::vector<std::vector<mytuple > > hit_info(lines);
  int divisor = L + 1;    // Includes the separator '#'
  //std::vector<boost::tuple<int, int> > pairs;
  //std::string pattern;

  // Bin the occurrences according to the sequence they are in.
  //  for (iterator it=string_to_contributions.begin(); it != string_to_contributions.end(); ++it) {  // iterator through strings in Hamming neighbourhood
  BOOST_FOREACH(std::string pattern, neighbourhood) {
    std::vector<long int> positions;
    bool is_palindrome = is_palindromic(pattern);
    sa.locate_iupac(pattern, positions);
    BOOST_FOREACH(int pos, positions) {
      int i = pos / divisor;  // index of the read containing the pos
      int j = pos % divisor;
      int dir = 1;
      if (i >= lines) {       // handle reverse complement
	i = 2*lines - i - 1;
	j = L - (j + k - 1) - 1;
	dir = -1;
      }
      if (not is_palindrome or count_palindromes_twice) {  // This branch allows a site to be counted twice, for both orientations
	hit_positions[i].push_back(j);
	hit_directions[i].push_back(dir);
	hit_contributions[i].push_back(pattern);
      }
      else {  // count palindrome only once
	if (std::find(hit_positions[i].begin(), hit_positions[i].end(), j) == hit_positions[i].end()) {
	  hit_positions[i].push_back(j);
	  hit_directions[i].push_back(dir);
	  hit_contributions[i].push_back(pattern);
	}
      }
    }
  }

  
  for (int i=0; i < lines; ++i) {
    int hit_count = hit_positions[i].size();
    std::vector<int> non_intersecting_occurrences;
    if (hit_count == 1)
      non_intersecting_occurrences.push_back(0);
    else {

      std::vector<int> hits(hit_count);               // This will contain indices to hit_positions/patterns vector
      std::vector<int> hamming_distances(hit_count);
      for (int index=0; index < hit_count; ++index) {
	hits[index] = index;
	if (use_two_strands)
	  hamming_distances[index] = std::min(hamming_distance(hit_contributions[i][index], seed),
					      hamming_distance(hit_contributions[i][index], seed_rev));
	else
	  hamming_distances[index] = hamming_distance(hit_contributions[i][index], seed);
      }
      std::sort(hits.begin(), hits.end(), [i, &hit_positions](int a, int b) { return hit_positions[i][a] < hit_positions[i][b];} );

      for (int current=0; current < hit_count; ++current) {
	int current_hd = hamming_distances[hits[current]];
	int current_pos = hit_positions[i][hits[current]];
	bool stop = false;
	int prev_index = current - 1;
	// Check the previous occurrences that are within cluster_threshold distance
	while (prev_index >= 0 and current_pos - hit_positions[i][hits[prev_index]] <= cluster_threshold) {
	  if (hamming_distances[hits[prev_index]] <= current_hd and hit_positions[i][hits[prev_index]] != current_pos) {
	    stop = true;
	    break;
	  }
	  --prev_index;
	}
	if (stop)
	  continue;
	int next_index = current + 1;
	// Check the next occurrences that are within cluster_threshold distance
	while (next_index < hit_count and hit_positions[i][hits[next_index]] - current_pos <= cluster_threshold) {
	  if (hamming_distances[hits[next_index]] <= current_hd and hit_positions[i][hits[next_index]] != current_pos) {
	    stop = true;
	    break;
	  }
	  ++next_index;
	}
	if (stop)
	  continue;
	int index = hits[current];
	hit_info[i].push_back(boost::make_tuple(hit_positions[i][index], hit_directions[i][index], hit_contributions[i][index]));
      }

    }

  }  // end i

  return helper(hit_info, k, seed, hamming_radius, use_grouped_dinucleotide_method_local);

}; // end adm_count_neighbourhood_contains_one

*/


dinuc_model<double>
right_extend_adm(const dinuc_model<double>& orig_adm, const std::vector<double>& bg, int extension)
{
  assert(extension >= 0);
  int orig_k = orig_adm.get_length();
  int k = orig_k + extension;
  dmatrix result(16, k);
  dmatrix matrix_bg(4,1);
  for (int b=0; b < 4; ++b)
    matrix_bg(b, 0) = bg[b];
  
  result.inject(orig_adm.dm, 0, 0);
  for (int j=orig_k; j < k; ++j) {
    for (int a=0; a < 4; ++a) {
      result.inject(matrix_bg, 4*a, j);
    }
  }

  dinuc_model<double> adm(result);
  return adm;
}


dinuc_model<double>
left_extend_adm(const dinuc_model<double>& orig_adm, const std::vector<double>& bg, int extension)
{
  assert(extension >= 0);
  int orig_k = orig_adm.get_length();
  int k = orig_k + extension;
  dmatrix result(16, k);
  dmatrix matrix_bg(4,1);
  for (int b=0; b < 4; ++b)
    matrix_bg(b, 0) = bg[b];
  
  result.inject(orig_adm.dm, 0, extension);
  for (int b=0; b < 4; ++b) {
    result(4+b, extension) = result(b, extension);
    result(8+b, extension) = result(b, extension);
    result(12+b, extension) = result(b, extension);
  }
    
  for (int j=0; j < extension; ++j) {
    int amax = j==0 ? 1 : 4;
    for (int a=0; a < amax; ++a) {
      result.inject(matrix_bg, 4*a, j);
    }
  }

  dinuc_model<double> adm(result);
  return adm;
}


dinuc_model<double>
force_adms_equal(const dinuc_model<double>& adm1, const dinuc_model<double>& adm2)
{
  int k = adm1.get_length();
  assert(k == adm2.get_length());

  dmatrix product(16, k);
  for (int row=0; row < 16; ++row) {
    for (int j=0; j < k; ++j) {
      product(row, j) = adm1.dm(row, j) * adm2.dm(row, j);
    }
  }

  dmatrix r(4, k+1);
  for (int a=0; a < 4; ++a) {
    r(a, k) = 1.0;
  }
    

  for (int j=k-1; j >= 0; --j) {
    int amax = j == 0 ? 1 : 4;
    for (int a=0; a < amax; ++a) {
      for (int b=0; b < 4; ++b) {
	r(a, j) += product(a*4+b, j) * r(b, j+1);
      }
    }
  }

  dmatrix result(16, k);
  for (int j=0; j < k; ++j) {
    int amax = j == 0 ? 1 : 4;
    for (int a=0; a < amax; ++a) {
      for (int b=0; b < 4; ++b) {
	if (r(a, j) > 0.0)
	  result(4*a+b, j) = product(4*a+b, j) * r(b, j+1) / r(a, j);
      }
    }
  }
  

  dinuc_model<double> adm(result);

  return adm;
}

dinuc_model<double>
dinuc_model_product(const dinuc_model<double>& adm1, const dinuc_model<double>& adm2, int d)
{
  int k1 = adm1.get_length();
  int k2 = adm2.get_length();
  int dimer_len = k1 + k2 + d;
  std::vector<double> bg(4, 0.25);
  dinuc_model<double> a1 = right_extend_adm(adm1, bg, dimer_len - k1);
  dinuc_model<double> a2 = left_extend_adm(adm2, bg, dimer_len - k2);

  return force_adms_equal(a1, a2);
}

template<typename T>
std::string
dinuc_model<T>::string_giving_max_probability(bool use_rna, bool use_iupac) const
{
  const char* nucs = use_rna ? "ACGU" : "ACGT";
  //  int k = length();
  std::string result(k, '-');
  //  if (require_directional_seed) {
  if (false) {                   // This did not seem to help with convergence of ID4
    // Tries to choose sequence with palindromic index greater than one as seed
    typedef std::pair<double, std::string> value_t;
    std::vector<value_t > probabilities;
    code_t size = pow(4,k);
    probabilities.reserve(size);
    for (code_t code=0; code < size; ++code) {
      const std::string& s = number_to_dna(code, k);
      probabilities.push_back(std::make_pair(probability(s), s));
    }
    std::sort(probabilities.begin(), probabilities.end(),
	      [](value_t a, value_t b) { return a.first >= b.first;}
	      );
    double probability;
    int counter=0;
    BOOST_REVERSE_FOREACH(boost::tie(probability, result), probabilities) {
      if (palindromic_index(result) > 1)
	break;
      ++counter;
    }
    printf("Counter is %i\n", counter);
  }
  else {
    if (use_iupac) {
      result = iupac_string_giving_max_probability(ip, use_rna);
      if (std::count(result.begin(), result.end(), 'N') / (double) k <= 1.0/3.0) // not too many Ns
	return result;
    }
    matrix<int> prev(4,k);
    for (int i=0; i < 4; ++i)
      prev(i, 0) = -1;
  
    std::vector<T> score_prev = ip.column(0);
    std::vector<T> score(4);
    for (int i=1; i<k; ++i) {
      for (int current=0; current < 4; ++current) {
	std::vector<T> temp(4);
	for (int p=0; p < 4; ++p)  // previous nucleotide
	  temp[p] = score_prev[p] * cond(i, p, current);
	int best = arg_max(temp);  // best previous nucleotide for this currect nucleotide
	score[current] = max_element(temp);
	prev(current, i) = best;
      }
      score_prev = score;
    }
    int current = arg_max(score);
    for (int i=k-1; i >= 0; --i) {
      result[i] = nucs[current];
      current = prev(current, i);
    }
  }
  
  return result;
}


/*
dinuc_model
adm_from_background(const std::vector<double>& bg, int k)
{
  assert(k >= 1);
  dmatrix result(16, k);
  dmatrix matrix_bg(4,1);
  for (int b=0; b < 4; ++b)
    matrix_bg(b, 0) = bg[b];
      
  for (int j=0; j < k; ++j) {
    int amax = j==0 ? 1 : 4;
    for (int a=0; a < amax; ++a) {
      result.inject(matrix_bg, 4*a, j);
    }
  }

  dinuc_model adm;
  adm.init(result);
  return adm;
}
*/

dinuc_model<double>
join_adms(const dinuc_model<double>& left, const dinuc_model<double>& right)
{

  int k1 = left.get_length();
  int k2 = right.get_length();
  int k = k1 + k2 - 1;
  dmatrix result(16, k);

  
  result.inject(right.representation(), 0, k1-1);
  result.inject(left.representation(), 0, 0);

  dinuc_model<double> adm(result);
  return adm;
}
