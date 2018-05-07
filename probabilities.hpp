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
#ifndef PROBABILITIES_HPP
#define PROBABILITIES_HPP

#include "matrix.hpp"
#include "data.hpp"
#include "parameters.hpp"
#include "common.hpp"

#include <string>
#include <vector>
#include <boost/tuple/tuple.hpp>

enum {atleast, atmost, exactly};

typedef double (*score_function_t)(const std::string&, const dmatrix&);

std::vector<double>
binomial_distribution(double p, int n);


std::vector<double>
binomial_distribution2(double p, int n);

std::vector<double>
poisson_binomial_distribution(const std::vector<double>& p);

double
poisson_binomial_expectation(const std::vector<double>& p);

struct condition
{
  condition(int n_, int type_) : n(n_), type(type_) {}

  int n;
  int type;
};



// compute the tail of a random variable X \in 0,1,2,...,n with distribution p
double
tail(const std::vector<double>& p, condition cond);



// compute the expectation of a random variable X \in 0,1,2,...,n with distribution p
double
expectation(const std::vector<double>& p);



double
entropy(const std::vector<double>& v);

double
KL_distance(const std::vector<double>& p, const std::vector<double>& q);

double
symmetric_KL_distance(const std::vector<double>& p, const std::vector<double>& q);

double
information_content(const std::vector<double>& p, 
		    std::vector<double> q = std::vector<double>());

double
average_information_content(const matrix<double>& m, std::vector<double> q = std::vector<double>());

double
matrix_information_content(const matrix<double>& m);

double
matrix_entropy(const matrix<double>& m);

double
matrix_KL_distance(const matrix<double>& m, const std::vector<double>& q);

double
matrix_KL_distance(const std::vector<double>& q, const matrix<double>& m);

double
aic_score(double maximum_log_likelihood, int lines, int k);

matrix<double>
matrix_to_logodds(const matrix<double>& m, const std::vector<double>& bg);

matrix<double>
matrix_to_weighted_logodds(const matrix<double>& m, const std::vector<double>& bg, double m_weight, double bg_weight);

dmatrix
matrix_to_affinity(const dmatrix& m);

double
compute_logodds_probability(const std::string& s, const dmatrix& m);

double
compute_normal_probability(const std::string& s, const dmatrix& m);


double
max_matrix_score(const dmatrix& m);

double
max_matrix_probability(const dmatrix& m);


std::string
string_giving_max_score(const dmatrix& m, bool use_rna);

matrix<double>
count_positional_background(const std::vector<std::string>& sequences);



template <typename T>
T
compute_bernoulli_probability(const std::string& line, const std::vector<double>& q)
{
  T prob = 1.0;
  for (int i=0; i < line.length(); ++i) {
    prob *= q[to_int(line[i])];
    assert(prob > 0);
  }

  return prob;
}

template <typename T, typename S>
T
compute_log_bernoulli_probability(const std::string& line, const std::vector<S>& q)
{
  T log_prob = 0.0;
  for (int i=0; i < line.length(); ++i) {
    log_prob += q[to_int(line[i])];
  }

  return log_prob;
}



template <typename T>
T
compute_bernoulli_probability(big_int code, int k, const std::vector<double>& q)
{
  T prob = 1.0;
  for (int i=0; i < k; ++i, code >>= 2)
    prob *= q[code & 3];
  
  assert(prob > 0);

  return prob;
}


// double
// compute_background_probability(const std::string& line, const std::vector<double>& q,
// 			       const matrix<double>& q2 = dmatrix());


// computes the probability P(line|\gamma)
template <typename T>
T
compute_background_probability(const std::string& line, 
			       const std::vector<double>& q,
			       const dmatrix& q2)
{
  assert(line.length() > 0);
  assert(q.size() == 4);

  T prob = 1.0;

  if (use_markov_background) {
    error(true, "Not implemented");
    // if (use_two_strands)
    //   return (compute_markov_probability(line, q, q2) + compute_markov_probability(reverse_complement(line), q, q2)) / 2.0;
    // else
    //   return compute_markov_probability(line, q, q2);
  } 
  else if (use_positional_background) {
    for (int i=0; i < line.length(); ++i)
      prob *= positional_background(to_int(line[i]),i);
    assert(prob > 0);
  }
  else { // 0th order model
    return compute_bernoulli_probability<T>(line, q); 
  }
    
  return prob;
}


// computes the logarithm of probability P(line|\gamma)
template <typename T, typename S>
T
compute_log_background_probability(const std::string& line, 
				   const std::vector<S>& q,
				   const dmatrix& q2)
{
  assert(line.length() > 0);
  assert(q.size() == 4);

  T log_prob = 0.0;

  if (use_markov_background) {
    error(false, "Not implemented!");
    // if (use_two_strands)
    //   return (compute_markov_probability(line, q, q2) + compute_markov_probability(reverse_complement(line), q, q2)) / 2.0;
    // else
    //   return compute_markov_probability(line, q, q2);
  } 
  else if (use_positional_background) {
    for (int i=0; i < line.length(); ++i)
      log_prob += positional_background(to_int(line[i]),i);
  }
  else { // 0th order model
    return compute_log_bernoulli_probability<T>(line, q); 
  }
    
  return log_prob;
}


// wrapper for the one below

// double
// compute_probability(const std::string& line, int pos, int direction,
// 		    const matrix<double>& m, 
// 		    const std::vector<double>& q,
// 		    const matrix<double>& q2 = dmatrix());

// Parameter direction tells the direction of the factor. It does NOT
// tell whether we chose the upper or lower strand.
// computes the probability P(line|\mathcal M(pos))
template <typename T>
T
compute_probability(const std::string& line_orig, const std::string& line_orig_rev,
		    int pos, int direction,
		    const matrix<double>& m, 
		    const std::vector<double>& q,
		    const matrix<double>& q2)
{
  assert(q[0]>0 && q[0]<=1);

  assert(q[1]>0 && q[1]<=1);
  assert(q[2]>0 && q[2]<=1);
  assert(q[3]>0 && q[3]<=1);
  assert(q.size() == 4);
  int k = m.get_columns();

  const std::string& line = direction == 1 ? line_orig : line_orig_rev;
  
  int L = line.length();
  if (direction == -1)
    pos = L - (pos + k);

  T prob = 1.0;

  // compute probability of motif part
  std::string seq = line.substr(pos,k);
  for (int i=0; i < k; ++i) {
    assert(pos+i < L);
    prob *= m(to_int(line[pos+i]),i);
  }

  if (use_markov_background) {
    // compute probability on the left side of motif
    prob *= q[to_int(line[0])];
    for (int i=0; i < pos - 1; ++i) {
      prob *= q2(to_int(line[i]),to_int(line[i+1]));
    }
    if (direction == -1) // This correction is because ADM(rev(X)) != rev(ADM)(X)
      prob *= (q[to_int(line[pos-1])] / q[to_int(line[0])]);
    // compute probability on the right side of motif
    prob *= q[to_int(line[pos+k])];
    for (int i=pos+k; i < L-1; ++i) {
      prob *= q2(to_int(line[i]), to_int(line[i+1]));
    }  
    if (direction == -1) // This correction is because ADM(rev(X)) != rev(ADM)(X)
      prob *= (q[to_int(line[L-1])] / q[to_int(line[pos+k])]);
  }
  else if (use_positional_background) {
    for (int i=0; i < line.length(); ++i) {
      if (i >= pos && i < pos + k)
	continue;
      prob *= positional_background(to_int(line[i]),i);
    }
  }
  else { // use 0th order model
    // compute probability on the left side of motif
    for (int i=0; i < pos; ++i) {
      prob *= q[to_int(line[i])];
    }
    // compute probability on the right side of motif
    for (int i=pos+k; i < line.length(); ++i) {
      prob *= q[to_int(line[i])];
    }  
  }

  assert(prob > 0.0);
  return prob;

}


// Parameter direction tells the direction of the factor. It does NOT
// tell whether we chose the upper or lower strand.
// computes the logarithm of probability P(line|\mathcal M(pos))
template <typename T, typename S>
T
compute_log_probability(const std::string& line_orig, const std::string& line_orig_rev,
			int pos, int direction,
			const matrix<T>& m, 
			const std::vector<S>& q,
			const matrix<double>& q2)
{
  /*
  assert(q[0]>0 && q[0]<=1);

  assert(q[1]>0 && q[1]<=1);
  assert(q[2]>0 && q[2]<=1);
  assert(q[3]>0 && q[3]<=1);
  */
  assert(q.size() == 4);
  int k = m.get_columns();

  const std::string& line = direction == 1 ? line_orig : line_orig_rev;
  
  int L = line.length();
  if (direction == -1)
    pos = L - (pos + k);

  T log_prob = 0.0;

  // compute probability of motif part
  std::string seq = line.substr(pos,k);
  for (int i=0; i < k; ++i) {
    assert(pos+i < L);
    log_prob += m(to_int(line[pos+i]),i);
  }

  if (use_markov_background) {
    // compute probability on the left side of motif
    log_prob += q[to_int(line[0])];
    for (int i=0; i < pos - 1; ++i) {
      log_prob += q2(to_int(line[i]),to_int(line[i+1]));
    }
    if (direction == -1) // This correction is because ADM(rev(X)) != rev(ADM)(X)
      log_prob += (q[to_int(line[pos-1])] / q[to_int(line[0])]);
    // compute probability on the right side of motif
    log_prob += q[to_int(line[pos+k])];
    for (int i=pos+k; i < L-1; ++i) {
      log_prob += q2(to_int(line[i]), to_int(line[i+1]));
    }  
    if (direction == -1) // This correction is because ADM(rev(X)) != rev(ADM)(X)
      log_prob += (q[to_int(line[L-1])] / q[to_int(line[pos+k])]);
  }
  else if (use_positional_background) {
    for (int i=0; i < line.length(); ++i) {
      if (i >= pos && i < pos + k)
	continue;
      log_prob += positional_background(to_int(line[i]),i);
    }
  }
  else { // use 0th order model
    // compute probability on the left side of motif
    for (int i=0; i < pos; ++i) {
      log_prob += q[to_int(line[i])];
    }
    // compute probability on the right side of motif
    for (int i=pos+k; i < line.length(); ++i) {
      log_prob += q[to_int(line[i])];
    }  
  }

  return log_prob;

}

template <typename T>
T
compute_probability(const std::string& line_orig, int pos, int direction,
		    const matrix<double>& m, 
		    const std::vector<double>& q,
		    const matrix<double>& q2)
{
  return compute_probability<T>(line_orig, reverse_complement(line_orig), pos, direction, m, q, q2);
}



template <typename T>
T
compute_dimer_probability(const std::string& line_orig, const std::string& line_orig_rev,
			  int pos, int direction,
			  const matrix<double>& m1, 
			  const matrix<double>& m2,
			  int d,
			  const std::vector<double>& q,
			  const matrix<double>& q2)
{
  assert(q.size() == 4);
  assert(q[0]>0 && q[0]<=1);
  assert(q[1]>0 && q[1]<=1);
  assert(q[2]>0 && q[2]<=1);
  assert(q[3]>0 && q[3]<=1);
  assert(d >= 0);

  int L = line_orig.length();
  
  T prob = 1.0;
  int k1 = m1.get_columns();
  int k2 = m2.get_columns();
  int dimer_len = k1 + k2 + d;

  const std::string& line = direction == 1 ? line_orig : line_orig_rev;
  
  if (direction == -1)
    pos = L - (pos + dimer_len);

  int pos2 = pos + dimer_len - k2;  // Starting position of the second motif
  
  // compute probability of the first motif part
  for (int i=0; i < k1; ++i) {
    assert(pos+i < line.length());
    prob *= m1(to_int(line[pos+i]), i);
  }

  // compute probability of the second motif part
  for (int i=0; i < k2; ++i) {
    assert(pos2+i < L);
    prob *= m2(to_int(line[pos2+i]), i);
  }

  if (use_markov_background) {
    // compute probability on the left side of motif
    prob *= q[to_int(line[0])];
    for (int i=0; i < pos-1; ++i) {
      prob *= q2(to_int(line[i]),to_int(line[i+1]));
    }  
    if (direction == -1) // This correction is because ADM(rev(X)) != rev(ADM)(X)
      prob *= (q[to_int(line[pos-1])] / q[to_int(line[0])]);

    // compute probability between the motifs
    if (d > 0) {
      prob *= q[to_int(line[pos+k1])];
      for (int i=pos+k1; i < pos2-1; ++i) {
	prob *= q2(to_int(line[i]), to_int(line[i+1]));
      }  
      if (direction == -1) // This correction is because ADM(rev(X)) != rev(ADM)(X)
	prob *= (q[to_int(line[pos2-1])] / q[to_int(line[pos+k1])]);
    }

    // compute probability on the right side of motif
    prob *= q[to_int(line[pos2+k2])];
    for (int i=pos2+k2; i < L-1; ++i) {
      prob *= q2(to_int(line[i]),to_int(line[i+1]));
    }  
    if (direction == -1) // This correction is because ADM(rev(X)) != rev(ADM)(X)
      prob *= (q[to_int(line[L-1])] / q[to_int(line[pos2+k2])]);

  }
  else if (use_positional_background) {
    // compute probability on the left side of the first motif
    for (int i=0; i < pos; ++i) 
      prob *= positional_background(to_int(line[i]), i);

    // compute probability between the motifs
    for (int i=pos+k1; i < pos2; ++i) 
      prob *= positional_background(to_int(line[i]), i);

    // compute probability on the right side of the second motif
    for (int i=pos2+k2; i < L; ++i)
      prob *= positional_background(to_int(line[i]), i);
  }
  else { // use 0th order model
    // compute probability on the left side of the first motif
    for (int i=0; i < pos; ++i) 
      prob *= q[to_int(line[i])];

    // compute probability between the motifs
    for (int i=pos+k1; i < pos2; ++i) 
      prob *= q[to_int(line[i])];

    // compute probability on the right side of the second motif
    for (int i=pos2+k2; i < L; ++i)
      prob *= q[to_int(line[i])];
  }

  assert(prob > 0.0);

  return prob;

} // compute_dimer_probability



// Parameter direction tells the direction of the dimer. It does NOT
// tell whether we chose the upper or lower strand.
// computes the spaced dimer probability P(line|\mathcal dimer(pos))
// Assumes that the PWMs are already correctly oriented
template <typename T>
T
compute_log_dimer_probability(const std::string& line_orig, const std::string& line_orig_rev,
			      int pos, int direction,
			      const matrix<T>& m1, 
			      const matrix<T>& m2,
			      int d,
			      const std::vector<T>& q,
			      const matrix<double>& q2)
{
  assert(q.size() == 4);
  /*
  assert(q[0]>0 && q[0]<=1);
  assert(q[1]>0 && q[1]<=1);
  assert(q[2]>0 && q[2]<=1);
  assert(q[3]>0 && q[3]<=1);
  */
  assert(d >= 0);

  int L = line_orig.length();
  
  T log_prob = 0.0;
  int k1 = m1.get_columns();
  int k2 = m2.get_columns();
  int dimer_len = k1 + k2 + d;

  const std::string& line = direction == 1 ? line_orig : line_orig_rev;
  
  if (direction == -1)
    pos = L - (pos + dimer_len);

  int pos2 = pos + dimer_len - k2;  // Starting position of the second motif
  
  // compute probability of the first motif part
  for (int i=0; i < k1; ++i) {
    assert(pos+i < line.length());
    log_prob += m1(to_int(line[pos+i]), i);
  }

  // compute probability of the second motif part
  for (int i=0; i < k2; ++i) {
    assert(pos2+i < L);
    log_prob += m2(to_int(line[pos2+i]), i);
  }

  if (use_markov_background) {
    // compute probability on the left side of motif
    log_prob += q[to_int(line[0])];
    for (int i=0; i < pos-1; ++i) {
      log_prob += q2(to_int(line[i]),to_int(line[i+1]));
    }  
    if (direction == -1) // This correction is because ADM(rev(X)) != rev(ADM)(X)
      log_prob += (q[to_int(line[pos-1])] / q[to_int(line[0])]);

    // compute probability between the motifs
    if (d > 0) {
      log_prob += q[to_int(line[pos+k1])];
      for (int i=pos+k1; i < pos2-1; ++i) {
	log_prob += q2(to_int(line[i]), to_int(line[i+1]));
      }  
      if (direction == -1) // This correction is because ADM(rev(X)) != rev(ADM)(X)
	log_prob += (q[to_int(line[pos2-1])] / q[to_int(line[pos+k1])]);
    }

    // compute probability on the right side of motif
    log_prob += q[to_int(line[pos2+k2])];
    for (int i=pos2+k2; i < L-1; ++i) {
      log_prob += q2(to_int(line[i]),to_int(line[i+1]));
    }  
    if (direction == -1) // This correction is because ADM(rev(X)) != rev(ADM)(X)
      log_prob += (q[to_int(line[L-1])] / q[to_int(line[pos2+k2])]);

  }
  else if (use_positional_background) {
    // compute log_probability on the left side of the first motif
    for (int i=0; i < pos; ++i) 
      log_prob += positional_background(to_int(line[i]), i);

    // compute log_probability between the motifs
    for (int i=pos+k1; i < pos2; ++i) 
      log_prob += positional_background(to_int(line[i]), i);

    // compute log_probability on the right side of the second motif
    for (int i=pos2+k2; i < L; ++i)
      log_prob += positional_background(to_int(line[i]), i);
  }
  else { // use 0th order model
    // compute log_probability on the left side of the first motif
    for (int i=0; i < pos; ++i) 
      log_prob += q[to_int(line[i])];

    // compute log_probability between the motifs
    for (int i=pos+k1; i < pos2; ++i) 
      log_prob += q[to_int(line[i])];

    // compute log_probability on the right side of the second motif
    for (int i=pos2+k2; i < L; ++i)
      log_prob += q[to_int(line[i])];
  }


  return log_prob;

} // compute_dimer_log_probability




template <typename T>
T
compute_dimer_probability(const std::string& line_orig, int pos, int direction,
			  const matrix<double>& m1, 
			  const matrix<double>& m2,
			  int d,
			  const std::vector<double>& q,
			  const matrix<double>& q2)
{
  return compute_dimer_probability<T>(line_orig, reverse_complement(line_orig), pos, direction, m1, m2, d, q, q2);
}


// double
// compute_probability(const std::string& line_orig, const std::string& line_orig_rev,
// 		    int pos, int direction,
// 		    const matrix<double>& m, 
// 		    const std::vector<double>& q,
// 		    const matrix<double>& q2 = dmatrix());





/*

// wrapper for the one below
double
compute_dimer_probability(const std::string& line_orig,
			  int pos, int direction,
			  const matrix<double>& m1, 
			  const matrix<double>& m2,
			  int d,
			  const std::vector<double>& q,
			  const matrix<double>& q2 = dmatrix());

double
compute_dimer_probability(const std::string& line_orig, const std::string& line_orig_rev,
			  int pos, int direction,
			  const matrix<double>& m1, 
			  const matrix<double>& m2,
			  int d,
			  const std::vector<double>& q,
			  const matrix<double>& q2 = dmatrix());

*/

// double
// compute_bernoulli_probability(const std::string& line, const std::vector<double>& q);

// double
// compute_bernoulli_probability(big_int code, int k, const std::vector<double>& q);




boost::tuple<std::vector<double>, dmatrix, std::vector<int> >
count_background(const std::vector<std::string>& sequences, bool use_rna=false);

#endif // PROBABILITIES_HPP
