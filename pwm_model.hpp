#ifndef PWM_MODEL_HPP
#define PWM_MODEL_HPP

#include "probabilities.hpp"
#include "common.hpp"
#include "matrix.hpp"
#include "matrix_tools.hpp"

#include <boost/shared_ptr.hpp>
#include <boost/make_shared.hpp>

#include <vector>
#include <string>


template <typename T=double>
class pwm_model : public binding_model<T>
{
public:


  pwm_model(const std::vector<std::string>& sequences, const std::string& seed, int n=2); // n tells the size of the hamming neighbourhood
  pwm_model(const std::string& filename);

  pwm_model() { }

  pwm_model(int k) : dm(4, k) { }

  pwm_model(const matrix<T>& dm_, bool normalize=true)
  { init(dm_, normalize); }

  boost::shared_ptr<binding_model<T> >
  clone() const;


  void
  init(const matrix<T>& dm_, bool normalize=true);     // initialize the position dependent first order model using mononucleotide count array

  matrix<T>
  representation() const;

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

  T
  probability(std::string::const_iterator begin,
	      std::string::const_iterator end,
	      int start_pos = 0) const; // start_pos is the starting position in the model

  T
  probability(const std::string& s, int start_pos = 0) const; // start_pos is the starting position in the model

  T
  log_probability(const std::string& s, int start_pos = 0) const; // start_pos is the starting position in the model

  T
  log_probability(std::string::const_iterator begin,
		  std::string::const_iterator end,
		  int start_pos = 0) const; // start_pos is the starting position in the model

  
  T
  score(const std::string& s, int start_pos = 0) const; // start_pos is the starting position in the model


  double
  distance(const binding_model<T>& other) const;


  boost::shared_ptr<binding_model<T> >
  cut(int start_pos, int width) const;
  

  boost::shared_ptr<binding_model<T> >
  reverse_complement() const;

  boost::shared_ptr<binding_model<FloatType> >
  log2() const;

  
  std::string
  string_giving_max_probability(bool use_rna) const;

  
  matrix<T> dm;  // 4 x k   mononucleotide array
};



double
compute_probability(const std::string& line_orig, int pos, int direction,
		    const matrix<double>& m, 
		    const std::vector<double>& q,
		    const matrix<double>& q2);

template <typename T>
void
pwm_model<T>::init(const matrix<T>& dm_, bool normalize)
{
  dm = dm_;
  if (normalize)
    normalize_matrix_columns(dm);
}



template <typename T>
matrix<T>
pwm_model<T>::representation() const
{
  return dm;
}



template <typename T>
int
pwm_model<T>::get_length() const
{
  return dm.get_columns();
}

template <typename T>
std::vector<double>
pwm_model<T>::information_content(std::vector<T> bg_model) const
{
  int k=get_length();
  std::vector<double> ic(k);
  for (int i=0; i < k; ++i)
    ic[i] = ::information_content(dm.column(i), bg_model);

  return ic;
}


template <typename T>
std::pair<int,int>
pwm_model<T>::dim() const
{
  return dm.dim();
}

template <typename T>
bool
pwm_model<T>::is_probability_model() const
{
  return is_column_stochastic_matrix(dm);
}



template <typename T>
void
pwm_model<T>::print(const std::string& header, const std::string& format, FILE* f) const    // Prints the counts in a 4 x k) matrix
{
  write_matrix(f, dm, header, format, false);
}



template <typename T>
T
pwm_model<T>::probability(std::string::const_iterator begin,
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
  std::string::const_iterator it;
  for (it=begin; it < std::min(end, begin + get_length() - start_pos); ++it)
    probability *= dm(to_int(*it), it - begin);

  assert(probability > 0.0);
 
  return probability;

}



template <typename T>
T
pwm_model<T>::log_probability(std::string::const_iterator begin,
			      std::string::const_iterator end,
			      int start_pos) const // start_pos is the starting position in the model
{
  assert(start_pos < get_length());
  std::string temp;
  T log_probability = 0.0;
  if (start_pos < 0) {
    begin += -start_pos;
    start_pos = 0;
  }
  std::string::const_iterator it;
  for (it=begin; it < std::min(end, begin + get_length() - start_pos); ++it)
    log_probability += dm(to_int(*it), it - begin);

  //  assert(probability > 0.0);
 
  return log_probability;

}

template <typename T>
T
pwm_model<T>::probability(const std::string& s, int start_pos) const
{
  return probability(s.begin(), s.end(), start_pos);
}

template <typename T>
T
pwm_model<T>::log_probability(const std::string& s, int start_pos) const
{
  return log_probability(s.begin(), s.end(), start_pos);
}

template <typename T>
T
pwm_model<T>::score(const std::string& s, int start_pos) const // start_pos is the starting position in the model
{
  return -100;
}

template <typename T>
double
pwm_model<T>::distance(const binding_model<T>& other) const
{
  return ::distance(dm, other.representation());
}

template <typename T>
boost::shared_ptr<binding_model<T> >
pwm_model<T>::cut(int start_pos, int width) const
{
  return boost::make_shared<pwm_model<T> >(dm.cut(0, start_pos, 4, width));
}

template <typename T>
boost::shared_ptr<binding_model<T> >
pwm_model<T>::reverse_complement() const
{
  return boost::make_shared<pwm_model<T> >(::reverse_complement(dm));
}

template <typename T>
boost::shared_ptr<binding_model<FloatType> >
pwm_model<T>::log2() const
{
  return boost::make_shared<pwm_model<FloatType> >(::log2<FloatType>(dm), false);
}

template <typename T>
boost::shared_ptr<binding_model<T> >
pwm_model<T>::clone() const
{
  return boost::make_shared<pwm_model<T> >(*this);
}


template <typename T>
std::string
pwm_model<T>::string_giving_max_probability(bool use_rna) const
{
  const char* nucs = use_rna ? "ACGU" : "ACGT";
  int k = dm.get_columns();
  std::string result(k, '-');
  for (int i=0; i<k; ++i)
    result[i] = nucs[arg_max(dm.column(i))];

  return result;
}

/*
double
compute_probability(const std::string& line_orig, int pos, int direction,
		    const matrix<double>& m, 
		    const std::vector<double>& q,
		    const matrix<double>& q2)
{
  return compute_probability(line_orig, reverse_complement(line_orig), pos, direction, pwm_model(m), q, q2);
}
*/

#endif // PWM_MODEL_HPP
