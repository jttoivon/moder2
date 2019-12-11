#ifndef DINUCLEOTIDE_HPP
#define DINUCLEOTIDE_HPP

/*

  Both in count matrix and in internal representation of the dinucleotide/adm model
  in the first column there are nonzero values only on the first four rows.

 */


#include "type.hpp"
#include "matrix.hpp"
#include "common.hpp"
#include "probabilities.hpp"
#include "matrix_tools.hpp"
#include "suffix_array_wrapper.hpp"

#include <string>
#include <vector>
#include <fstream>

#include <boost/make_shared.hpp>

extern std::vector<std::string> dinuc_headers;
extern const char* mono_headers[];

// Position dependent first order model for a site of length k
template <typename T=double>
class dinuc_model : public binding_model<T>
{
public:


  // n tells the size of the hamming neighbourhood
  
  //dinuc_model(const std::vector<std::string>& sequences, const std::string& seed, int n=2); 
  dinuc_model(const std::string& filename);

  dinuc_model(int k_);
  
  dinuc_model() : binding_model<T>() { }

  dinuc_model(const matrix<T>& dm_, bool normalize=true) { init(dm_, normalize); }

  // This is mainly for the log2 method. No normalization!
  dinuc_model(const matrix<T>& dm_, const matrix<T>& ip_)
    : k(dm_.get_columns()), dm(dm_), ip(ip_) { }

  boost::shared_ptr<binding_model<T> >
  clone() const;

  // initialize the position dependent first order model using dinucleotide count array
  void
  init(const matrix<T>& dm_, bool normalize=true);     

  matrix<T>
  representation() const;

  double
  cond(int i, int c1, int c2) const;  // If character in position i is c1, the what is the probability of getting c2 in the next pos 
  double
  cond(int i, char a, char b) const;

  double
  initial(int c, int i) const;  // What is the probability of having character c in pos i


  double
  initial(char a, int i) const;  // What is the probability of having character a in pos i
  
  int
  get_length() const;   // returns k, the width of the binding site

  bool
  is_probability_model() const;

  std::vector<double>
  information_content(std::vector<T> bg = std::vector<T>(4, 0.25)) const;

  std::pair<int,int>
  dim() const;
  
  void
  print(const std::string& header, const std::string& format, FILE* f) const;  

  void
  read(const std::string& filename);
  
  T
  probability(const std::string& s, int start_pos = 0) const; // start_pos is the starting position in the model

  T
  probability(std::string::const_iterator begin,
	      std::string::const_iterator end,
	      int start_pos = 0) const; // start_pos is the starting position in the model

  T
  log_probability(const std::string& s, int start_pos = 0) const; // start_pos is the starting position in the model

  T
  log_probability(std::string::const_iterator begin,
		  std::string::const_iterator end,
		  int start_pos = 0) const; // start_pos is the starting position in the model

  
  T
  score(const std::string& s, int start_pos = 0) const; // start_pos is the starting position in the model

  double
  distance(const binding_model<T>&) const;

  boost::shared_ptr<binding_model<T> >
  cut(int start_pos, int width) const;

  boost::shared_ptr<binding_model<T> >
  reverse_complement() const;

  boost::shared_ptr<binding_model<FloatType> >
  log2() const;

  
  std::string
  string_giving_max_probability(bool use_rna, bool use_iupac) const;

  int k;
  matrix<T> dm;  // 16 x k   transition probabilities
  matrix<T> ip;  // 4 x k    initial probabilities to each state
};



bool
almost_equal(const dmatrix& dm1, const dmatrix& dm2, double threshold = 0.0000001);

bool
almost_equal(double d1, double d2, double threshold = 0.0000001);

bool
operator==(const dinuc_model<double>& dm1, const dinuc_model<double>& dm2);

std::ostream&
operator<<(std::ostream& str, const dinuc_model<double>& adm);

std::vector<dmatrix>
dinucleotide_counts_scan_better(const std::string& seed, const std::vector<std::string>& sequences, int n,
				model_type model_type);

std::vector<dmatrix>
dinucleotide_counts_suffix_array(const std::string& seed, const std::vector<std::string>& sequences,
				 const suffix_array& sa, int n);


template <typename T>
void
generate_initial_probabilities(const matrix<T>& dm, matrix<T>& ip)
{
  assert(dm.get_columns() == ip.get_columns());
  
  int k = dm.get_columns();
  for (int c=0; c < 4; ++c) {
    ip(c, 0) = dm(c, 0);
  }
  
  // generate rest of the initial probabilities
  for (int i=1; i < k; ++i) {
    std::vector<T> temp(4, 0.0);
    for (int a=0; a<4; ++a) {
      for (int b=0; b<4; ++b) {
	temp[b] += ip(a, i-1) * dm(a*4+b, i);
      }
    }
    if (sum(temp) > 0.0)
      ip.set_column(i, normalize_vector_copy(temp));
    else
      ip.set_column(i, std::vector<T>(4, 0.0));
  }
}

// transforms dinucleotide counts in the dmatrix(16,k) to 
// conditional_probabilities dmatrix(16,k)
template <typename T>
matrix<T>
conditional_probabilities(const matrix<T>& dm)
{
  assert(dm.get_rows() == 16);
  int k = dm.get_columns();
  matrix<T> result(16, k);

  int count = 0;
  // do normalization
  for (int i=0; i < k; ++i) {     
    for (int first=0; first < 4; ++first) {
      if (i==0 and first != 0)    // Initial probabilities in the first column. Start symbols is also 'A'
	break;
      double sum = 0;
      for (int second=0; second < 4; ++second)
	sum += dm((first<<2) + second, i);
      if (sum > 0) {
	for (int second=0; second < 4; ++second)
	  result((first<<2) + second, i) = dm((first<<2) + second, i) / sum;
      }
      else {
	for (int second=0; second < 4; ++second)
	  result((first<<2) + second, i) = 0.0;//0.25;
	++count;

      }
    }
  }

  //  if (count != 0)
  //    printf("In conditional_probabilities: division by zero %i times\n", count);
  return result;
}


template <typename T>
boost::shared_ptr<binding_model<FloatType> >
dinuc_model<T>::log2() const
{
  // The first column in dm contains zeros!
  return boost::make_shared<dinuc_model<FloatType> >(::log2_special<FloatType>(dm), ::log2<FloatType>(ip));
}


template <typename T>
boost::shared_ptr<binding_model<T> >
dinuc_model<T>::cut(int position, int len) const
{
  return boost::make_shared<dinuc_model<T> >(sub_adm(*this, position, len));
}

template <typename T>
T
weighted_distance(const dinuc_model<T>& model1, const dinuc_model<T>& model2)
{
  T distance = 0.0;
  int k = model1.get_length();
  for (int c=1; c < k; ++c) {
    for (int r=0; r < 16; ++r) {
      distance = std::max(distance, std::abs(model1.ip(r/4, c) * model1.dm(r, c) - model2.ip(r/4, c) * model2.dm(r, c)));
    }
  }
  for (int r=0; r < 4; ++r)    // initial probabilities
    distance = std::max(distance, std::abs(model1.dm(r, 0) - model2.dm(r, 0)));
  return distance;
}

template <typename T>
double
dinuc_model<T>::distance(const binding_model<T>& other) const
{
  //return ::distance(dm, other.representation());
  return weighted_distance(dynamic_cast<const dinuc_model&>(*this), dynamic_cast<const dinuc_model&>(other));
}



template <typename T>
matrix<T>
dinuc_model<T>::representation() const
{
  return dm;
}


template <typename T>
boost::shared_ptr<binding_model<T> >
dinuc_model<T>::clone() const
{
  return boost::make_shared<dinuc_model<T> >(*this);
}

template <typename T>
std::pair<int,int>
dinuc_model<T>::dim() const
{
  return std::make_pair(16, k);
}

template <typename T>
dinuc_model<T>
reverse_complement(const dinuc_model<T>& dm)
{
  int k = dm.k;
  
  // reverse complement of the dinucleotides
  int transform[] = {15, 11, 7, 3, 14, 10, 6, 2,
		     13, 9, 5, 1, 12, 8, 4, 0};

  dinuc_model<T> result(k);
  result.ip = reverse_complement(dm.ip);
  for (int j=1; j < k; ++j) {
    for (int ab=0; ab < 16; ++ab) {
      int a = ab / 4;
      int b = ab % 4;
      T divisor = dm.ip(b, j);
      if (divisor != 0.0)
	result.dm(transform[ab], k-j) = dm.dm(ab, j) * dm.ip(a, j-1) / divisor;
      else
	result.dm(transform[ab], k-j) = 0.0;
    }
  }

  for (int a=0; a < 4; ++a)   // Put the initial probabilities in the first column of the transition matrix
    result.dm(a, 0) = result.ip(a, 0);

  return result;
}


template <typename T>
boost::shared_ptr<binding_model<T> >
dinuc_model<T>::reverse_complement() const
{
  return boost::make_shared<dinuc_model<T> >(::reverse_complement(*this));
}

template <typename T>
double
information_content_of_adm_column(const dinuc_model<T>& adm, int col, const std::vector<T>& bg)
{
  int k = adm.get_length();
  assert(col >= 0);
  assert(col < k);
  double ic=0.0;
  std::vector<T> v(4);
  if (col == 0) {
    for (int b=0; b < 4; ++b)
      v[b] = adm.dm(b, 0);
    ic = information_content(v, bg);
  }
  else {
    for (int a=0; a < 4; ++a) {
      for (int b=0; b < 4; ++b)
	v[b] = adm.dm(4*a+b, col);
      ic += adm.ip(a, col-1) * information_content(v, bg);
    }
  }
  return ic;
}

template<typename T>
std::vector<double>
dinuc_model<T>::information_content(std::vector<T> bg) const
{
  std::vector<double> ic(k);
  for (int i=0; i < k; ++i) 
    ic[i] = information_content_of_adm_column(*this, i, bg);
  return ic;
}

template<typename T>
bool
dinuc_model<T>::is_probability_model() const
{
  int k = this->k;
  if (ip.dim() != std::make_pair(4, k))
    return false;
  if (dm.dim() != std::make_pair(16, k))
    return false;
  
  // Check initial probabilities
  for (int i=0; i < k; ++i) {
    if (not almost_equal(sum(ip.column(i)), 1.0)) {
      printf("Error in initial probabilities\n");
      return false;
    }
  }

  // Check conditional probabilities
  for (int i=0; i < k; ++i) {
    for (int a=0; a < 4; ++a) {
      if (i==0 and a != 0)   // 'A' is also the begin symbol
	break;
      dvector temp(4, 0.0);
      for (int b=0; b < 4; ++b)
	temp[b]=dm(4*a+b, i);
      if (not almost_equal(sum(temp), 1.0) and sum(temp) != 0) {
	printf("Error in conditional probabilities\n");
	return false;
      }
    }
  }  
  return true;
}

template<typename T>
void
dinuc_model<T>::init(const matrix<T>& dm_, bool normalize)
{
  assert(dm_.get_rows() == 16);
  k = dm_.get_columns();
  //dm=dm_;
  dm = conditional_probabilities(dm_);
  ip = matrix<T>(4,k);
  // transforms dinucleotide counts in the dmatrix(16,k) to 
  // initial_probabilities dmatrix(4,k)

  generate_initial_probabilities(dm, ip);
  //normalize_matrix_columns(result);
}


template<typename T>
dinuc_model<T>::dinuc_model(int k_)
  : k(k_), dm(16, k_), ip(4,k_)
{}


template<typename T>
void
dinuc_model<T>::read(const std::string& filename) 
{
  std::vector<std::string> lines;
  std::ifstream stream(filename);
  if (not stream.is_open()) {
    fprintf(stderr, "Could not open file %s\n", filename.c_str());
    exit(1);
  }
    
  std::string line;
  while (getline(stream, line)) {
    lines.push_back(line);
  }
  assert(lines.size() == 16+4);
  int k = split(lines[0], '\t').size() - 1;
  dmatrix ip2(4, k);
  dmatrix dm2(16, k);
  //  double pseudo_count = 0.000001;
  double pseudo_count = 0.0;
  for (int i=0; i < 16; ++i) {
    std::vector<std::string> fields = split(lines[i], '\t');
    assert(fields.size() == k+1);
    assert(fields[k-1] == "ADM_DI");
    assert(fields[k] == dinuc_headers[i]);
    for (int j=0; j < k-1; ++j)
      dm2(i, j+1) = atof(fields[j]) + pseudo_count;
  }
  for (int i=0; i < 4; ++i) {
    std::vector<std::string> fields = split(lines[16+i], '\t');
    assert(fields.size() == k+1);
    assert(fields[k] == mono_headers[i]);
    for (int j=0; j < k; ++j)
      ip2(i, j) = atof(fields[j]);
  }
  //ip = ip2;  // don't use the original initial probabilities (except in the first column), but regenerate them.
  for (int a = 0; a < 4; ++a)
    dm2(a, 0) = ip2(a, 0);
  
  init(dm2);
  
  //  init(read_matrix_file(filename));
}

template<typename T>
dinuc_model<T>::dinuc_model(const std::string& filename)
{
  read(filename);
}

/*
template<typename T>
dinuc_model<T>::dinuc_model(const std::vector<std::string>& sequences, const std::string& seed, int n) 
{
  typedef std::vector<dmatrix> (*func_ptr_t)(const std::string&, const std::vector<std::string>&,
					     const suffix_array&, int);
  int k = seed.length();
  func_ptr_t func_ptr;
  if (k >= 20 or n >= 4) 
    ;//func_ptr = dinucleotide_counts_scan;
  else
    func_ptr = dinucleotide_counts_suffix_array;


  init(func_ptr(seed, sequences, n)[0]);
}
*/

template<typename T>
// i refers to the index of the second character 'c2'
double
dinuc_model<T>::cond(int i, int c1, int c2) const
{
  int dinuc = (c1 << 2) + c2;
    
  return dm(dinuc, i);
}

template<typename T>
// i refers to the index of the second character 'b'
double
dinuc_model<T>::cond(int i, char a, char b) const
{
  return dinuc_model::cond(i, to_int(a), to_int(b));
}

template<typename T>
double
dinuc_model<T>::initial(int c, int i) const  // What is the probability of having character c in pos i
{
  return ip(c, i);
}

template<typename T>
double
dinuc_model<T>::initial(char a, int i) const  // What is the probability of having character a in pos i
{
  return dinuc_model::initial(to_int(a), i);
}

template<typename T>
int
dinuc_model<T>::get_length() const
{ return k; }


template<typename T>
void
dinuc_model<T>::print(const std::string& header, const std::string& format, FILE* fp) const
{
  
  fprintf(fp, "%s", header.c_str());
  
  int k = get_length();
  bool col_headers = false;
  if (col_headers) {
    for (int i=0; i < k-1; ++i) {
      fprintf(fp, "\t%i", i);
    }
    fprintf(fp, "\n");
  }
   
  for (int c=0; c < 16; ++c) {
    for (int i=0; i < k-1; ++i) {
      fprintf(fp, format.c_str(), dm(c, i+1));
      fprintf(fp, "\t");
    }
    fprintf(fp, "%s\t%s\n", "ADM_DI", dinuc_headers[c].c_str());
  }
  for (int c=0; c < 4; ++c) {
    for (int i=0; i < k; ++i) {
      fprintf(fp, format.c_str(), ip(c, i));
      fprintf(fp, "\t");
    }
    fprintf(fp, "%s\n", mono_headers[c]);
  }
}

template<typename T>
// start_pos tells the position in the dinucleotide model where scoring the kmer begins
T
dinuc_model<T>::score(const std::string& kmer, int start_pos) const
{
  assert(start_pos < get_length());
  std::string temp;
  double score = 0.0;
  if (start_pos < 0) {
    start_pos = -start_pos;
    temp = kmer.substr(start_pos);
    start_pos = 0;
  }
  else
    temp = kmer;
  int len = temp.length();
  score += std::log2(ip(to_int(temp[0]), start_pos) / 0.25);
  for (int i=0; i < std::min(len, get_length() - start_pos)-1; ++i)
    score += std::log2(cond(start_pos+i+1, to_int(temp[i]), to_int(temp[i+1])) / 0.25);

  return score;
} 

template <typename T>
T
dinuc_model<T>::probability(std::string::const_iterator begin,
			    std::string::const_iterator end,
			    int start_pos) const // start_pos is the starting position in the model
{
  assert(start_pos < get_length());
  std::string temp;
  T probability = 1.0;
  if (start_pos < 0) {
    begin += -start_pos;
    start_pos = 0;
  }
  typedef std::string::const_iterator iterator;
  std::string::const_iterator it=begin;
  
  probability *= ip(to_int(*it++), start_pos);
  iterator end2 = std::min(end, begin + get_length() - start_pos);
  for (; it < end2; ++it)
    probability *= cond(start_pos + it - begin, to_int(*(it-1)), to_int(*it));

  assert(probability > 0.0);
 
  return probability;

}

template<typename T>
T
dinuc_model<T>::probability(const std::string& kmer, int start_pos) const
{
  assert(start_pos < get_length());
  std::string temp;
  double prob = 1.0;
  if (start_pos < 0) {
    start_pos = -start_pos;
    temp = kmer.substr(start_pos);
    for (int i=0; i < start_pos; ++i)
      prob *= 0.25;
    start_pos = 0;
  }
  else
    temp = kmer;
  int len = temp.length();
  prob *= ip(to_int(temp[0]), start_pos);
  for (int i=0; i < std::min(len, get_length() - start_pos)-1; ++i)
    prob *= cond(start_pos+i+1, to_int(temp[i]), to_int(temp[i+1]));
  for (int i=std::min(len, get_length() - start_pos); i < len; ++i)
    prob *= 0.25;
  return prob;
} 

template <typename T>
T
dinuc_model<T>::log_probability(std::string::const_iterator begin,
				std::string::const_iterator end,
				int start_pos) const // start_pos is the starting position in the model
{
  assert(start_pos < get_length());
  std::string temp;
  T probability = 0.0;
  if (start_pos < 0) {
    begin += -start_pos;
    start_pos = 0;
  }
  typedef std::string::const_iterator iterator;
  std::string::const_iterator it=begin;
  
  probability += ip(to_int(*it++), start_pos);
  iterator end2 = std::min(end, begin + get_length() - start_pos);
  for (; it < end2; ++it)
    probability += cond((it - begin)+start_pos, to_int(*(it-1)), to_int(*it));

  return probability;

}

template<typename T>
T
dinuc_model<T>::log_probability(const std::string& kmer, int start_pos) const
{
  return log_probability(kmer.begin(), kmer.end(), start_pos);
}

dinuc_model<double>
right_extend_adm(const dinuc_model<double>& orig_adm, const std::vector<double>& bg, int extension);


dinuc_model<double>
left_extend_adm(const dinuc_model<double>& orig_adm, const std::vector<double>& bg, int extension);


template<typename T>
dinuc_model<T>
sub_adm(const dinuc_model<T>& orig_adm, int position, int length)
{
  int orig_k = orig_adm.get_length();
  assert(length > 0);
  assert(position >= 0);
  assert(position + length <= orig_k);

  matrix<T> result = orig_adm.dm.cut(0, position, 16, length);
  for (int a=0; a < 4; ++a)
    result(a, 0) = orig_adm.ip(a, position);
  for (int a=4; a < 16; ++a)
    result(a, 0) = 0.0;

  dinuc_model<T> adm(result);
  return adm;
  
}

dinuc_model<double>
force_adms_equal(const dinuc_model<double>& adm1, const dinuc_model<double>& adm2);

dinuc_model<double>
dinuc_model_product(const dinuc_model<double>& adm1, const dinuc_model<double>& adm2, int d);

dinuc_model<double>
join_adms(const dinuc_model<double>& left, const dinuc_model<double>& right);

/*
dinuc_model
reverse_complement(const dinuc_model& dm);

dmatrix
dinucleotide_reverse_complement(const dmatrix& dm);


dinuc_model
pwm_to_dinucleotide(const dmatrix& pwm);


dmatrix
dinucleotide_counts_scan(const std::string& seed, const std::vector<std::string>& sequences, int n);
*/


#endif // DINUCLEOTIDE_HPP
