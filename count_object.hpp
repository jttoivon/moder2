#include "matrix.hpp"
#include "probabilities.hpp"
#include "iupac.hpp"
#include "matrix_tools.hpp"
#include "kmer_tools.hpp"
#include "pwm_model.hpp"
#include "dinucleotide.hpp"

#include <vector>
#include <cassert>

extern bool use_multinomial;
extern int hamming_radius;
extern bool extra_debug;
// use_multinomial, hamming_radius, adm, seed_bias
enum model_type { ppm=0, adm=1, adm_fixed=2};

const char* model_type_strings[] = {"ppm", "adm", "adm-fixed"};

dmatrix
correct_seed_bias(const std::vector<dmatrix>& grouped_dinucleotide_counts, const std::string& seed, int hamming_radius);

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
      type = adm;
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
    assert(type == adm or type == adm_fixed or type == ppm);
    
    if (type==ppm)
      return boost::make_shared<pwm_model<> >(normalize_matrix_columns_copy(counts[0]));
    else if (type==adm or seed.length()==0)  // for gapped motifs no seed is used.
      return boost::make_shared<dinuc_model<> >(counts[0]);
    else {
      if (false) {
	// Zero-out too small counts
	double absolute_min = 5.0;
	for (int r2=0; r2 < counts.size(); ++r2) {
	  for (int j=0; j < length; ++j) {
	    double s = ::sum(counts[r2].cut(0, j, 16, 1));
	    for (int ab=0; ab < 16; ++ab) {
	      double& c = counts[r2](ab, j);
	      if (c < 0.01*s or c < absolute_min)
		c = 0.0;
	    }
	  }
	}
      }
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

  void
  add_pseudo_counts(const prior<double>& pseudo_counts, const dinucleotide_prior<double>& dinucleotide_pseudo_counts)
  {
    if (type == ppm)
      pseudo_counts.add(counts[0]);
    else 
      for (int r2=0; r2 < counts.size(); ++r2) {
	dinucleotide_pseudo_counts.add(counts[r2]);
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
    else
      positions = 0;   // update nothing
    //  printf("HD is %i %i %s %s\n", hd, mismatches, print_bitvector(mismatches).c_str(), print_bitvector(positions).c_str());
    BitString mask = static_cast<BitString>(1)<<(w-1);
    if (type == ppm) {
      for (int pos=0; pos < w; ++pos, mask>>=1) {
	if (positions & mask)
	  counts[0](to_int(substr[pos]), pos) += z; // update columns of pwm marked by bit vector positions
      }
    }
    else {
      int shift = 2*w-2;
      BitString code = dna_to_number<BitString>(substr);
      if (type == adm) {
	for (int pos=0; pos < w; ++pos, mask>>=1, shift -= 2) {
	  if (positions & mask) {
	    counts[0]((code>>shift)&static_cast<BitString>(15), pos) += z; // update columns of pwm marked by bit vector positions
	  }
	}
      }
      else {
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
    else
      positions = 0;   // update nothing
    int first = w1-1;  // last position of the first half-site
    int last = w1+d;   // first position of the second half-site
    BitString mask = static_cast<BitString>(1)<<(d+w2);
    if (type == ppm) {
      for (int pos=first; pos <= last; ++pos, mask>>=1) {
	if (positions & mask)
	  counts[0](to_int(substr[pos]), pos) += z; // update all columns
      }
    }
    else {
      BitString code = dna_to_number<BitString>(substr);
      int w = d+w2;
      BitString mask2 = (static_cast<BitString>(1) << (w*2)) - 1;
      int shift = 2*w-2;
      code &= mask2;  // zero-out the prefix bits
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

dmatrix
correct_seed_bias(const std::vector<dmatrix>& grouped_dinucleotide_counts, const std::string& seed, int hamming_radius)
{
  int k = seed.length();
  
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
    if (j == k-1) {
      for (int t=0; t <= tmax; ++t)
	for (int b=0; b < 4; ++b) 
	  tau[b][t] = 1.0;
    }
    else {
      for (int t=0; t <= tmax; ++t) {
	for (int b=0; b < 4; ++b)
	  tau[b][t] = 0.0;
	int suffix_len = k-j-1;
	std::vector<std::string> neighbourhood = get_n_neighbourhood(seed.substr(j+1, suffix_len), t);
	BOOST_FOREACH(std::string w, neighbourhood) {
	  double temp = 1.0;
	  for (int i=1; i < w.length(); ++i) {
	    temp *= corrected(to_int(w[i-1])*4+to_int(w[i]), i+j+1);
	  }
	  for (int b=0; b < 4; ++b) {
	    int a = (b << 2) + to_int(w[0]);
	    tau[b][t] += corrected(a, j+1) * temp;
	  }
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
	double divisor = tau[a & 3][std::min(hamming_radius-r-1, k-j-1)];
	if (divisor != 0.0)
	  N[a] += quotient(a, r) = grouped_dinucleotide_counts[r](a, j) / divisor;
	else
	  N[a] += quotient(a, r) = grouped_dinucleotide_counts[r](a, j);
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
  }
  return corrected;
  
} // correct_seed_bias
