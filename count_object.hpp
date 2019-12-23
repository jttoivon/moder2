#include "matrix.hpp"
#include "probabilities.hpp"
#include "iupac.hpp"
#include "matrix_tools.hpp"
#include "kmer_tools.hpp"
#include "pwm_model.hpp"
#include "dinucleotide.hpp"

#include <vector>
#include <cassert>

bool esko_tau_method=true;

extern bool use_multinomial;
extern int hamming_radius;
extern bool extra_debug;
// use_multinomial, hamming_radius, adm, seed_bias


const char* model_type_strings[] = {"ppm", "adm-unfixed", "adm-fixed"};

dmatrix
correct_seed_bias(const std::vector<dmatrix>& grouped_dinucleotide_counts, const std::string& seed, int hamming_radius);

std::map<std::string, std::vector<dmatrix> > pseudo_count_table_cache;

class count_object
{
public:

  count_object() {type = ppm;}
  
  count_object(model_type type_, int width_)
  {
    int rows = type_ == ppm ? 4 : 16;
    int r2 = type_ == adm_fixed ? hamming_radius : 1;
    counts.assign(r2, dmatrix(rows, width_));
    type = type_;
    length = width_;
  }

  count_object(const std::vector<dmatrix>& counts_)
    : counts(counts_)
  {
    assert(counts.size() > 0);
    int rows, cols;
    boost::tie(rows, cols) = counts[0].dim();
    if (counts.size() > 1)
      type = adm_fixed;
    else if (rows==4)
      type = ppm;
    else
      type = adm_unfixed;
    length = cols;

      
  }

  /*
  boost::shared_ptr<binding_model<> >
  normalized_pwm()
  {
    assert(type == ppm);
    return boost::shared_ptr<binding_model<> >(new pwm_model<>(normalize_matrix_columns_copy(counts[0])));
  }
  */
  
  boost::shared_ptr<binding_model<> >
  normalized(const std::string& seed)
  {
    assert(type == adm_unfixed or type == adm_fixed or type == ppm);
    
    if (type==ppm)
      return boost::make_shared<pwm_model<> >(normalize_matrix_columns_copy(counts[0]));
    else if (type==adm_unfixed or seed.length()==0)  // for gapped motifs no seed is used.
      return boost::make_shared<dinuc_model<> >(counts[0]);
    else {
      /*
      if (false) {
	int count=0;
	// Zero-out too small counts
	double absolute_min = 5.0;
	for (int r2=0; r2 < counts.size(); ++r2) {
	  for (int j=0; j < length; ++j) {
	    double s = ::sum(counts[r2].cut(0, j, 16, 1));
	    for (int ab=0; ab < 16; ++ab) {
	      double& c = counts[r2](ab, j);
	      //	      if ((c < 0.01*s or c < absolute_min) and c > 0.0) {
	      if ((c < 0.01*s) and c > 0.0) {
		++count;
		printf("r=%i, ab=%i, j=%i, count=%e\n", r2, ab, j, c);
		c = 0.0;
	      }
	    }
	  }
	}
	printf("Cut %i numbers.\n", count);
      }
      */
      return boost::make_shared<dinuc_model<> >(correct_seed_bias(counts, seed, hamming_radius));
    }
  }

  /*
  dinuc_model
  normalized_dinuc_model() {
    if (type == count_adm_fixed)
      ;
    else if (type == count_adm)
      return dinuc_model(counts[0]);
    else
      error(true, "Unknown parameter");
  }
  */

  void
  fill_with(double value)
  {
    for (int r2=0; r2 < counts.size(); ++r2) {
      counts[r2].fill_with(value);
    }
  }
  
  // this creates pseudo count tables for adm-fixed model
  template <typename BitString>
  std::vector<dmatrix>
  //void
  create_pseudo_count_tables(const std::string& seed, const std::vector<double>& background_probabilities)
  {
    std::vector<std::string> dinucleotides = {
      "AA", "AC", "AG", "AT",
      "CA", "CC", "CG", "CT",
      "GA", "GC", "GG", "GT",
      "TA", "TC", "TG", "TT",
    };
    std::vector<double> dinucleotide_probability(16);
    for (int a=0; a < 4; ++a) {
      for (int b=0; b < 4; ++b) {
	dinucleotide_probability[4*a+b] = background_probabilities[a] * background_probabilities[b];
      }
    }
    
    std::string nucs="ACGT";
    bool force_multinomial = true;
    const int w = seed.length();
    //double z = pow(0.25, w); // weight of a sequence
    std::vector<dmatrix> pscounts(hamming_radius, dmatrix(16, length));
    //pcounts.assign(hamming_radius, dmatrix(16, length));
    std::vector<std::string> neighbourhood = get_n_neighbourhood(seed, hamming_radius, false);
    BOOST_FOREACH(std::string substr, neighbourhood) {
      std::vector<double> v = compute_bernoulli_iupac_probability<double>(substr, background_probabilities);
      double z = std::accumulate(v.begin(), v.end(), 1.0, std::multiplies<double>());
      BitString mismatches = iupac_mismatch_positions<BitString>(substr, seed);
      int hd = mypopcount(mismatches);
      BitString positions = 0;
      if (not force_multinomial or hd < hamming_radius)
	positions = ~static_cast<BitString>(0);  // update all
      else if (hd == hamming_radius)  // update only mismatch positions
	positions = mismatches;
      else
	continue; //positions = 0;   // update nothing
      //  printf("HD is %i %i %s %s\n", hd, mismatches, print_bitvector(mismatches).c_str(), print_bitvector(positions).c_str());
      BitString mask = static_cast<BitString>(1)<<(w-1);  // contains 1 in the wth position from right
      //int shift = 2*(w-1);
      //BitString code = dna_to_number<BitString>(substr);
      int r=0;
      if (positions & mask) {
	for (int b=0; b < 4; ++b) {
	  if (iupac_match(nucs[b], substr[0])) {
	    pscounts[r](b, 0) += 0.01*z/v[0]*background_probabilities[b];
	  }
	}
	if (mismatches&mask)
	  ++r;
      }
      mask>>=1;
      for (int pos=1; pos < w; ++pos, mask>>=1) {
	if (positions & mask) {
	  double z2 = z / (v[pos-1] * v[pos]);
	  std::string di = substr.substr(pos-1, 2);
	  for (int ab=0; ab < 16; ++ab) {
	    if (iupac_string_match(dinucleotides[ab], di))
	      pscounts[r](ab, pos) += 0.01*z2*dinucleotide_probability[ab]; // update columns of pwm marked by bit vector positions
	  }
	}
	if (mismatches&mask)
	  ++r;
      }
    }
    /*
    // find the average column sum
    std::vector<double> average_column_sums;
    int column_counter = 0; // number of nonzero columns
    for (int r=0; r < hamming_radius; ++r) {
      for (int pos=0; pos < w; ++pos) {
	double s = ::sum(pscounts[r].column(pos));
	if (s > 0.0) {
	  average_column_sums.push_back(s);
	  ++column_counter;
	}
      }
    }
    double average = ::sum(average_column_sums) / (float) column_counter;
    double scaling_factor = 0.01 / average;
    assert(scaling_factor > 0.0);
    for (int r=0; r < hamming_radius; ++r) {
      pscounts[r] = scaling_factor * pscounts[r];
    }
    */
    return pscounts;
  }
  
  // this creates pseudo count tables for adm-fixed model
  template <typename BitString>
  std::vector<dmatrix>
  //void
  create_pseudo_count_tables_slow(const std::string& seed, const std::vector<double>& background_probabilities)
  {
    bool force_multinomial = true;
    const int w = seed.length();
    //double z = pow(0.25, w); // weight of a sequence
    std::vector<dmatrix> pscounts(hamming_radius, dmatrix(16, length));
    //pcounts.assign(hamming_radius, dmatrix(16, length));
    std::vector<std::string> neighbourhood = get_n_neighbourhood(seed, hamming_radius, true);
    BOOST_FOREACH(std::string substr, neighbourhood) {
      double z = 0.01 * compute_bernoulli_probability<double>(substr, background_probabilities);
      BitString mismatches = iupac_mismatch_positions<BitString>(substr, seed);
      int hd = mypopcount(mismatches);
      BitString positions = 0;
      if (not force_multinomial or hd < hamming_radius)
	positions = ~static_cast<BitString>(0);  // update all
      else if (hd == hamming_radius)  // update only mismatch positions
	positions = mismatches;
      else
	continue; //positions = 0;   // update nothing
      //  printf("HD is %i %i %s %s\n", hd, mismatches, print_bitvector(mismatches).c_str(), print_bitvector(positions).c_str());
      BitString mask = static_cast<BitString>(1)<<(w-1);  // contains 1 in the wth position from right
      int shift = 2*(w-1);
      BitString code = dna_to_number<BitString>(substr);
      int r=0;
      for (int pos=0; pos < w; ++pos, mask>>=1, shift -= 2) {
	if (positions & mask) {
	  pscounts[r]((code>>shift)&15, pos) += z; // update columns of pwm marked by bit vector positions
	}
	if (mismatches&mask)
	  ++r;
      }
    }
    return pscounts;
  }
  
  void
  add_pseudo_counts(const prior<double>& pseudo_counts,
		    const dinucleotide_prior<double>& dinucleotide_pseudo_counts,
		    const std::string& seed)
  {
    bool use_new_pseudo_counts = true;
    if (type == ppm)
      pseudo_counts.add(counts[0]);
    else {
      if (type == adm_fixed && use_new_pseudo_counts && seed != "") {
	auto iter = pseudo_count_table_cache.find(seed);
	if (iter != pseudo_count_table_cache.end()) {
	  printf("Found %s in cache\n", seed.c_str());
	  pcounts = iter->second;
	}
	else {
	  printf("Did not found %s in cache\n", seed.c_str());
	  std::vector<dmatrix> temp = create_pseudo_count_tables<myuint128>(seed, dinucleotide_pseudo_counts.background_probabilities);
	  pseudo_count_table_cache[seed] = temp;
	  pcounts = temp;
	}
	for (int r2=0; r2 < counts.size(); ++r2) {
	  counts[r2] += pcounts[r2];
	}
      }
      else {
	for (int r2=0; r2 < counts.size(); ++r2) {
	  dinucleotide_pseudo_counts.add(counts[r2]);
	}
      }
    }
  }

  
  template <typename BitString>
  void
  add_sequence(const std::string& substr, const std::string& seed, bool force_multinomial, double z)
  {
    int w = substr.length();
    BitString mismatches = iupac_mismatch_positions<BitString>(substr, seed);
    int hd = mypopcount(mismatches);
    BitString positions = 0;
    if (not force_multinomial or hd < hamming_radius)
      positions = ~static_cast<BitString>(0);  // update all
    else if (hd == hamming_radius)  // update only mismatch positions
      positions = mismatches;
    else {
      return;
      //positions = 0;   // update nothing, SHOULDN'T THIS JUST RETURN
    }
    //  printf("HD is %i %i %s %s\n", hd, mismatches, print_bitvector(mismatches).c_str(), print_bitvector(positions).c_str());
    BitString mask = static_cast<BitString>(1)<<(w-1);  // contains 1 in the wth position from right
    if (type == ppm) {
      for (int pos=0; pos < w; ++pos, mask>>=1) {
	if (positions & mask)
	  counts[0](to_int(substr[pos]), pos) += z; // update columns of pwm marked by bit vector positions
      }
    }
    else {
      int shift = 2*(w-1);
      BitString code = dna_to_number<BitString>(substr);
      if (type == adm_unfixed) {
	// If pos == 0, then the dinucleotide is (A, substr[0])
	for (int pos=0; pos < w; ++pos, mask>>=1, shift -= 2) {
	  if (positions & mask) {
	    counts[0]((code>>shift)&static_cast<BitString>(15), pos) += z; // update columns of pwm marked by bit vector positions
	  }
	}
      }
      else {  // type == adm_fixed
	int r=0;
	for (int pos=0; pos < w; ++pos, mask>>=1, shift -= 2) {
	  if (positions & mask) {
	    counts[r]((code>>shift)&15, pos) += z; // update columns of pwm marked by bit vector positions
	  }
	  if (mismatches&mask)
	    ++r;
	}
      }
    }
  }

  template <typename BitString>
  void
  add_gap_sequence(const std::string& substr, const std::string& seed, int d, int w1, int w2,
		   bool force_multinomial, double z)
  {
    assert(substr.length() == seed.length());
  
    BitString mismatches = iupac_mismatch_positions<BitString>(substr, seed);
    int hd = mypopcount(mismatches);
    BitString positions = 0;
    if (not force_multinomial or hd <= hamming_radius)
      positions = ~static_cast<BitString>(0);  // update all
    else {
      return;
      //positions = 0;   // update nothing, SHOULDN'T THIS JUST RETURN
    }
    int first = w1-1;  // last position of the first half-site
    int last = w1+d;   // first position of the second half-site
    BitString mask = static_cast<BitString>(1)<<(d+w2);  // 1-bit is in the 'first' position
    if (type == ppm) {
      for (int pos=first; pos <= last; ++pos, mask>>=1) {
	if (positions & mask)
	  counts[0](to_int(substr[pos]), pos) += z; // update all columns
      }
    }
    else {
      BitString code = dna_to_number<BitString>(substr);
      int w = d+w2+1;
      //      BitString mask2 = (static_cast<BitString>(1) << (w*2)) - 1;
      //      code &= mask2;  // zero-out the prefix bits
      int shift = 2*(w-1);   // shifting this much to right puts the bits of char in the 'first' position to bits 1 and 0
      for (int pos=first; pos <= last; ++pos, mask>>=1, shift -= 2) {
	if (positions & mask)
	  counts[0]((code>>shift)&15, pos) += z; // update all columns
      }
    }
  }

  int
  get_length() const
  { return length; }
  
  double
  sum() const
  {
    double s = 0.0;
    for (int r2=0; r2 < counts.size(); ++r2) {
      s += ::sum(counts[r2]);
    }
    return s;
  }

  void
  write_counts(FILE* f, const std::string& str, const std::string& format)
  {
    for (int r2=0; r2 < counts.size(); ++r2) {
      if (type == adm_fixed)
	fprintf(f, "r=%i\n", r2);
      write_matrix(f, counts[r2], str, format);
    }
  }
  
  count_object&
  operator+=(const count_object& rhs);
  
  //private:
  std::vector<dmatrix> counts;
  std::vector<dmatrix> pcounts;  // pseudo counts
  model_type type;
  int length;
};


count_object&
count_object::operator+=(const count_object& rhs)
{
  assert(counts.size() == rhs.counts.size());
  for (int r2=0; r2 < counts.size(); ++r2) {
    counts[r2] += rhs.counts[r2];
  }
  return *this;
}

double
count_helper(double count)
{
  return std::min(count/10.0, 1.0);
}

double
esko_count_helper(double count)
{
  if (count <= 1.0)
    return 0.0;
  double relative_standard_deviation = 1.0/sqrt(count);
  double result = 1.0 - relative_standard_deviation;
  assert(result <= 1.0);
  assert(result >= 0.0);
  return result;
}

double
lowest_dinucleotide_count(const std::string& w, const dmatrix& all_counts)
{
  assert(is_nucleotide_string(w));
  int k = all_counts.get_columns();
  double result = std::numeric_limits<double>::max();
  int start_pos = k-w.length()-1;
  for (int i=1; i < w.length(); ++i) {
    result = std::min(all_counts(to_int(w[i-1])*4+to_int(w[i]), start_pos+i+1), result);
  }
  return result;
}


// computes the maximum of minimum over all sequences defined by the iupac
std::vector<double>
lowest_dinucleotide_count_iupac(const std::string& w, const dmatrix& all_counts)
{
  assert(is_iupac_string(w));
  int k = all_counts.get_columns();
  std::vector<double> result(4, std::numeric_limits<double>::lowest());
  int start_pos = k-w.length()-1;
  //    std::vector<double> result(w.length(), std::numeric_limits<double>::max());
  dmatrix A(4, w.length());
  A.fill_with(std::numeric_limits<double>::max());
  for (int i=w.length()-1; i >= 1; --i) {
    for (char a : std::string(iupac_class(w[i-1]))) {
      double temp = std::numeric_limits<double>::lowest();
      for (char b : std::string(iupac_class(w[i]))) {
	double dcount = all_counts(to_int(a)*4+to_int(b), start_pos+i+1);
	temp = std::max(std::min(dcount, A(to_int(b), i)),
			temp);
      }
      A(to_int(a), i-1) = temp;
    }
  }
  //int last = w.length()-1;
  for (char b : std::string(iupac_class(w[0])))
    result[to_int(b)] = A(to_int(b), 0);
  return result;
}

std::tuple<std::vector<double>, std::vector<double> >
compute_bias_and_low_counts(const std::string& s, int t, int j, const dmatrix& all_counts, const dmatrix& corrected)
{
  dmatrix array(4, s.length());
  std::vector<double> tau(4, 0.0);
  std::vector<double> low_counts(4, std::numeric_limits<double>::lowest());
  std::vector<std::string> neighbourhood = get_n_neighbourhood(s, t, true);
  BOOST_FOREACH(std::string w, neighbourhood) {
    double low_count = lowest_dinucleotide_count(w, all_counts);
    double temp = 1.0;
    for (int i=1; i < w.length(); ++i) {
      temp *= corrected(to_int(w[i-1])*4+to_int(w[i]), i+j+1);
    }
    for (int b=0; b < 4; ++b) {
      int a = (b << 2) + to_int(w[0]);  // bw_0
      tau[b] += corrected(a, j+1) * temp;
      double low_count2 = std::min(low_count, all_counts(a, j+1));
      low_counts[b] = std::max(low_count2, low_counts[b]);
    }
	  
  }  // end BOOST_FOREACH
  return std::make_tuple(tau, low_counts);
}

std::tuple<std::vector<double>, std::vector<double> >
compute_bias_and_low_counts_iupac(const std::string& s, int t, int j, const dmatrix& all_counts, const dmatrix& corrected)
{
  dmatrix array(4, s.length());
  std::vector<double> tau(4, 0.0);
  std::vector<double> low_counts(4, std::numeric_limits<double>::lowest());
  std::vector<std::string> neighbourhood = get_n_neighbourhood(s, t, false);
  BOOST_FOREACH(std::string w, neighbourhood) {
    std::vector<double> low_count = lowest_dinucleotide_count_iupac(w, all_counts);
    //double temp = 1.0;
    array.fill_with(0.0);
    for (int a=0; a < 4; ++a)
      array(a, w.length()-1) = 1.0;
    for (int i=w.length()-2; i >= 0; --i) {
      for (char a : std::string(iupac_class(w[i]))) {
	for (char b : std::string(iupac_class(w[i+1]))) {
	  int ab = to_int(a)*4 + to_int(b);
	  array(to_int(a), i) += corrected(ab, j+i+2) * array(to_int(b), i+1);
	}
      }
    }

    for (int b=0; b < 4; ++b) {
      for (char c : std::string(iupac_class(w[0]))) {
	int bc = (b << 2) + to_int(c);  // bw_0
	tau[b] += corrected(bc, j+1) * array(to_int(c), 0);
	double low_count2 = std::min(all_counts(bc, j+1), low_count[to_int(c)]);
	low_counts[b] = std::max(low_count2, low_counts[b]);
      }
    }
	  
  }  // end BOOST_FOREACH
  return std::make_tuple(tau, low_counts);
}

dmatrix
correct_seed_bias(const std::vector<dmatrix>& grouped_dinucleotide_counts, const std::string& seed, int hamming_radius)
{
  int k = seed.length();
  dmatrix all_counts(16, k);
  for (int r=0; r < grouped_dinucleotide_counts.size(); ++r)
    all_counts += grouped_dinucleotide_counts[r];
  
  //  write_matrix(stdout, motif, to_string("%s dinucleotide motif matrix counts:\n", name.c_str()), "%.0f");
  dmatrix corrected(16, k);
  for (int j=k-1; j >= 0; --j) {
    if (extra_debug) {
      printf("\n");
      printf("j=%i\n", j);
      printf("=============================\n");
    }
    
    int tmax = std::min(hamming_radius-1, k-j-1);          // Compute the correction factors tau[b][t]
    boost::multi_array<double, 2> tau(boost::extents[4][tmax+1]);
    boost::multi_array<double, 2> tau2(boost::extents[4][tmax+1]);   // testing, REMOVE
    boost::multi_array<double, 2> low_counts(boost::extents[4][tmax+1]);
    boost::multi_array<double, 2> low_counts2(boost::extents[4][tmax+1]);   // testing, REMOVE
    if (j == k-1) {
      for (int t=0; t <= tmax; ++t)
	for (int b=0; b < 4; ++b) 
	  tau[b][t] = 1.0;
    }
    else {
      for (int t=0; t <= tmax; ++t) {
	for (int b=0; b < 4; ++b) {
	  tau[b][t] = 0.0;
	  low_counts[b][t] = std::numeric_limits<double>::lowest();
	  low_counts2[b][t] = std::numeric_limits<double>::lowest();   // testing, REMOVE
	}
	int suffix_len = k-j-1;
	std::vector<double> my_tau;
	std::vector<double> my_low_counts;
	std::tie(my_tau, my_low_counts) = compute_bias_and_low_counts_iupac(seed.substr(j+1, suffix_len), t, j, all_counts, corrected);
	for (int b=0; b < 4; ++b) {
	  tau[b][t] = my_tau[b];
	  low_counts[b][t] = my_low_counts[b];
	}
      }
    }

    if (extra_debug) {
      //print_array_with_default_headers(stdout, tau);
      int row_begin, row_end;
      int col_begin, col_end;
      boost::tie(row_begin, row_end, col_begin, col_end) = get_ranges(tau);
      //    std::vector<std::string> row_headers = integer_range(row_begin, row_end);
      std::vector<std::string> row_headers = { "A", "C", "G", "T"};
      //std::vector<std::string> col_headers = integer_range(col_begin, col_end);
      std::vector<std::string> col_headers;
      for (int t=0; t < col_end; ++t)
	col_headers.push_back(to_string("t=%i", t));
      printf("tau:\n");
      print_array(stdout, tau, 
		  row_headers,
		  col_headers, "%e");
    }
    
    std::vector<double> N(16);                            // Correct for seed bias
    int rmax = std::min(hamming_radius-1, j);
    dmatrix quotient(16, rmax+1);
    for (int a=0; a < 16; ++a) {
      for (int r=0; r <= rmax; ++r) {
	int tr = std::min(hamming_radius-r-1, k-j-1);
	double divisor = tau[a & 3][tr];
	//	double alpha = j == k-1 ? 1.0 : count_helper(all_counts((a&3)*4 + to_int(seed[j+1]), j+1));
	double alpha;
	if (esko_tau_method)
	  alpha = j == k-1 ? 1.0 : esko_count_helper(low_counts[a&3][tr]);
	else
	  alpha = j == k-1 ? 1.0 : count_helper(low_counts[a&3][tr]);
	//printf("Beta=%.4f, j=%i, row=%i, r=%i\n", alpha, j, a, r);
	N[a] += quotient(a, r) = grouped_dinucleotide_counts[r](a, j) / (divisor*alpha + (1-alpha));
	/*
	if (divisor != 0.0)
	  N[a] += quotient(a, r) = grouped_dinucleotide_counts[r](a, j) / divisor;
	else
	  N[a] += quotient(a, r) = grouped_dinucleotide_counts[r](a, j);
	*/
      }
    }
    if (extra_debug) {
      std::vector<std::string> col_headers2;
      for (int r=0; r <= rmax; ++r)
	col_headers2.push_back(to_string("r=%i", r));
      printf("Quotient\n");
      print_matrix(stdout, quotient, dinuc_headers, col_headers2, "%f");
      //write_matrix(stdout, quotient, to_string("Quotient ", j));
    }

    int amax = j == 0 ? 1 : 4;
    for (int a=0; a < amax; ++a) {                           // Normalize
      double sum = 0.0;
      double m = 0.0;
      for (int b=0; b < 4; ++b) {
	int ab = (a << 2) + b;
	m = std::max(m, N[ab]);
      }
      for (int b=0; b < 4; ++b) {
	int ab = (a << 2) + b;
	//N[ab] += pseudo_count * m * background_probabilities[b];
	sum += N[ab];
      }
      for (int b=0; b < 4; ++b) {
	int ab = (a << 2) + b;
	if (sum == 0.0)
	  corrected(ab, j) = 0;
	else 
	  corrected(ab, j) = N[ab] / sum;
      }
    }
    if (extra_debug)
      write_matrix(stdout, corrected, "Corrected");
  } // end for j
  return corrected;
  
} // correct_seed_bias
