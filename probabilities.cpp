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
#include "probabilities.hpp"
#include "parameters.hpp"
#include "common.hpp"
#include "iupac.hpp"
#include "combinatorics.hpp"
#include "my_assert.hpp"

#include <cassert>
#include <cfloat>

std::vector<double>
binomial_distribution(double p, int n)
{
  std::vector<double> result(n+1);
  for (int k=0; k < n+1; ++k) {
    result[k] = choose(n,k)*pow(p,k)*pow(1-p,n-k);
  }
  return result;
}

std::vector<double>
binomial_distribution2(double p, int n)
{
  std::vector<double> V(n+1,0);
  std::vector<double> Vnew(n+1);
  V[0] = 1;
  for (int l=0; l < n; ++l) {// go through all sequences
    my_assert2(p, 0, >=);
    my_assert2(p, 1, <=);
    Vnew[0] = V[0]*(1-p);  // initialisation
    for (int s=1; s <= l+1; ++s)// go through all numbers of double occurrences
      Vnew[s] = V[s-1]*p + V[s]*(1-p);
    V = Vnew;
  }

  return Vnew;   // distribution of the number of double occurences in the dataset
}

double
poisson_binomial_expectation(const std::vector<double>& p)
{
  double expected = 0.0;
  for (int i=0; i < p.size(); ++i)
    expected += p[i];

  return expected;
}

std::vector<double>
poisson_binomial_distribution(const std::vector<double>& p)
{
  int n = p.size();
  std::vector<double> V(n+1,0);
  std::vector<double> Vnew(n+1);
  V[0] = 1;
  for (int l=0; l < n; ++l) {// go through all sequences
    my_assert2(p[l], 0, >=);
    my_assert2(p[l], 1, <=);
    Vnew[0] = V[0]*(1-p[l]);  // initialisation
    for (int s=1; s <= l+1; ++s)// go through all numbers of double occurrences
      Vnew[s] = V[s-1]*p[l] + V[s]*(1-p[l]);
    std::swap(V, Vnew);
  }

  return V;   // distribution of the number of double occurences in the dataset
}

// compute the tail of a random variable X \in 0,1,2,...,n with distribution p
double
tail(const std::vector<double>& p, condition cond)
{
  assert(cond.n < p.size());
  my_assert2(sum(p), 1.0 + 100000*DBL_EPSILON, <=);

  double sum = 0.0;
  if (cond.type == atmost) {
    for (int i=0; i <= cond.n; ++i)
      sum += p[i];
  }
  else {
    for (int i=cond.n; i < p.size(); ++i)
      sum += p[i];
  }

  my_assert2(sum, 1.0, <=);
  my_assert2(sum, 0.0, >=);
  return sum;
}


// compute the expectation of a random variable X \in 0,1,2,...,n with distribution p
double
expectation(const std::vector<double>& p)
{
  double e = 0.0;
  for (int i=0; i < p.size(); ++i)
    e += i*p[i];
  return e;
}

double
entropy(const std::vector<double>& v)
{
  double sum=0;
  for (int i=0; i < v.size(); ++i)
    if (v[i] != 0)
      sum += log2(v[i]) * v[i];
  return -sum;
}

double
KL_distance(const std::vector<double>& p, const std::vector<double>& q)
{
  assert(p.size() == q.size());
  double sum=0;
  for (int i=0; i < p.size(); ++i)
    if (p[i] !=0)
      sum += p[i]*log2(p[i]/q[i]);
  return sum;
}

double
symmetric_KL_distance(const std::vector<double>& p, const std::vector<double>& q)
{
  return std::min(KL_distance(p,q), KL_distance(q,p));
}

double
information_content(const std::vector<double>& p, 
		    std::vector<double> q)
{
  if (q.size() == 0)
    q = std::vector<double>(p.size(), 0.25);
  assert(p.size() == q.size());
  return KL_distance(p, q);
}

double
average_information_content(const matrix<double>& m, std::vector<double> q)
{
  if (q.size() == 0)
    q = std::vector<double>(m.get_rows(), 0.25);
  assert(m.get_rows() == q.size());
  double sum=0;
  for (int j=0; j < m.get_columns(); ++j)
    sum += information_content(m.column(j), q);
  return sum/m.get_columns();
}

double
matrix_information_content(const matrix<double>& m)
{
  double ic=0;
  for (int j=0; j < m.get_columns(); ++j)
    ic += information_content(m.column(j));
  return ic;
}

double
matrix_entropy(const matrix<double>& m)
{
  double ic=0;
  for (int j=0; j < m.get_columns(); ++j)
    ic += entropy(m.column(j));
  return ic;
}

double
matrix_KL_distance(const matrix<double>& m, const std::vector<double>& q)
{
  double KL=0;
  for (int j=0; j < m.get_columns(); ++j)
    KL += KL_distance(m.column(j), q);
  return KL;
}

double
matrix_KL_distance(const std::vector<double>& q, const matrix<double>& m) 
{
  double KL=0;
  for (int j=0; j < m.get_columns(); ++j)
    KL += KL_distance(q, m.column(j));
  return KL;
}


matrix<double>
matrix_to_logodds(const matrix<double>& m, const std::vector<double>& bg)
{
  int k=m.get_columns();
  matrix<double> logodds_motif(4,k);
  for (int i=0; i < 4; ++i)
    for (int j=0; j < k; ++j)
      logodds_motif(i,j) = log2(m(i,j)/bg[i]);
  return logodds_motif;
}

matrix<double>
matrix_to_weighted_logodds(const matrix<double>& m, const std::vector<double>& bg, double m_weight, double bg_weight)
{
  int k=m.get_columns();
  matrix<double> logodds_motif(4,k);
  for (int i=0; i < 4; ++i)
    for (int j=0; j < k; ++j)
      logodds_motif(i,j) = log2(m(i,j) / bg[i]) + log2(m_weight / bg_weight);
  return logodds_motif;
}

dmatrix
matrix_to_affinity(const dmatrix& m)
{
  int k=m.get_columns();
  dmatrix affinity_motif(4,k);
  for (int j=0; j < k; ++j) {
    double max = max_element(m.column(j));
    for (int i=0; i < 4; ++i) {
      affinity_motif(i,j) = m(i,j) / max;
    }
  }
  return affinity_motif;
}

double
compute_logodds_probability(const std::string& s, const dmatrix& m)
{
  int k = m.get_columns();
  assert(s.length() == k);
  double result=0.0;
  for (int i=0; i < k; ++i) {
    if (to_int(s[i]) == -1)
      result += -DBL_MAX;
    else
      result += m(to_int(s[i]), i);
  }
  return result;
}

double
compute_normal_probability(const std::string& s, const dmatrix& m)
{
  int k = m.get_columns();
  double result=1.0;
  assert(s.length() == k);
  for (int i=0; i < k; ++i)
    result *= m(to_int(s[i]), i);
  return result;
}

double
max_matrix_score(const dmatrix& m)
{
  int k = m.get_columns();
  double max_sum=0.0;
  for (int c=0; c<k; ++c) {
    max_sum += max_element(m.column(c));
  }
  return max_sum;
}

double
max_matrix_probability(const dmatrix& m)
{

  int k = m.get_columns();
  double max_prod=1.0;
  for (int c=0; c<k; ++c) {
    max_prod *= max_element(m.column(c));
  }
  return max_prod;
}


std::string
string_giving_max_score(const dmatrix& m)
{
  const char* nucs = "ACGT";
  int k = m.get_columns();
  std::string result(k, '-');
  for (int i=0; i<k; ++i)
    result[i] = nucs[arg_max(m.column(i))];

  return result;
}


double
aic_score(double maximum_log_likelihood, int lines, int k)
{
  return 2*(2*lines + 4 + 4*k) - 2*maximum_log_likelihood;
}

matrix<double>
count_positional_background(const std::vector<std::string>& sequences)
{
  int L=sequences[0].length();
  assert(L>0);
  matrix<double> freq(4,L);

  for (int i=0; i < sequences.size(); ++i) {
    const std::string& line = sequences[i];
    for (int j=0; j < L; ++j) {
      ++freq(to_int(line[j]),j);
      //      ++freq(to_int(complement(line[j])), L-j-1); // reverse complement
    }
  }

  return freq;

}


double
compute_bernoulli_probability(const std::string& line, const std::vector<double>& q)
{
  double prob = 1.0;
  for (int i=0; i < line.length(); ++i) {
    prob *= q[to_int(line[i])];
    assert(prob > 0);
  }

  return prob;
}

double
compute_log_bernoulli_probability(const std::string& line, const std::vector<double>& q)
{
  double log_prob = 0.0;
  for (int i=0; i < line.length(); ++i) {
    log_prob += q[to_int(line[i])];
    assert(log_prob > 0);
  }

  return log_prob;
}

double
compute_bernoulli_probability(big_int code, int k, const std::vector<double>& q)
{
  double prob = 1.0;
  for (int i=0; i < k; ++i, code >>= 2)
    prob *= q[code & 3];
  
  assert(prob > 0);

  return prob;
}


double
compute_markov_probability(const std::string& line, 
			       const std::vector<double>& q,
			       const matrix<double>& q2)
{
  double prob = q[to_int(line[0])];
  for (int i=0; i < line.length()-1; ++i)
    prob *= q2(to_int(line[i]),to_int(line[i+1]));
  
  assert(prob > 0.0);
  return prob;
}

// computes the probability P(line|\gamma)
double
compute_background_probability(const std::string& line, 
			       const std::vector<double>& q,
			       const dmatrix& q2)
{
  assert(line.length() > 0);
  assert(q.size() == 4);

  double prob = 1.0;

  if (use_markov_background) {
    if (use_two_strands)
      return (compute_markov_probability(line, q, q2) + compute_markov_probability(reverse_complement(line), q, q2)) / 2.0;
    else
      return compute_markov_probability(line, q, q2);
  } 
  else if (use_positional_background) {
    for (int i=0; i < line.length(); ++i)
      prob *= positional_background(to_int(line[i]),i);
    assert(prob > 0);
  }
  else { // 0th order model
    return compute_bernoulli_probability(line, q); 
  }

    
  return prob;
}


// computes the logarithm of probability P(line|\gamma)
double
compute_background_log_probability(const std::string& line, 
				   const std::vector<double>& q,
				   const dmatrix& q2)
{
  assert(line.length() > 0);
  assert(q.size() == 4);

  double log_prob = 0.0;

  if (use_markov_background) {
    error(false, "Not implemented!");
    if (use_two_strands)
      return (compute_markov_probability(line, q, q2) + compute_markov_probability(reverse_complement(line), q, q2)) / 2.0;
    else
      return compute_markov_probability(line, q, q2);
  } 
  else if (use_positional_background) {
    for (int i=0; i < line.length(); ++i)
      log_prob += positional_background(to_int(line[i]),i);
    assert(log_prob > 0);
  }
  else { // 0th order model
    return compute_log_bernoulli_probability(line, q); 
  }

    
  return log_prob;
}

// Parameter direction tells the direction of the factor. It does NOT
// tell whether we chose the upper or lower strand.
// computes the probability P(line|\mathcal M(pos))
double
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

  double prob = 1.0;

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
double
compute_log_probability(const std::string& line_orig, const std::string& line_orig_rev,
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

  double log_prob = 0.0;

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

  assert(log_prob > 0.0);
  return log_prob;

}

double
compute_probability(const std::string& line_orig, int pos, int direction,
		    const matrix<double>& m, 
		    const std::vector<double>& q,
		    const matrix<double>& q2)
{
  return compute_probability(line_orig, reverse_complement(line_orig), pos, direction, m, q, q2);
}


// Parameter direction tells the direction of the dimer. It does NOT
// tell whether we chose the upper or lower strand.
// computes the spaced dimer probability P(line|\mathcal dimer(pos))
// Assumes that the PWMs are already correctly oriented
double
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
  
  double prob = 1.0;
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
double
compute_dimer_log_probability(const std::string& line_orig, const std::string& line_orig_rev,
			      int pos, int direction,
			      const matrix<double>& m1, 
			      const matrix<double>& m2,
			      int d,
			      const std::vector<double>& q,
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
  
  double log_prob = 0.0;
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

  assert(log_prob > 0.0);

  return log_prob;

} // compute_dimer_probability

double
compute_dimer_probability(const std::string& line_orig, int pos, int direction,
			  const matrix<double>& m1, 
			  const matrix<double>& m2,
			  int d,
			  const std::vector<double>& q,
			  const matrix<double>& q2)
{
  return compute_dimer_probability(line_orig, reverse_complement(line_orig), pos, direction, m1, m2, d, q, q2);
}



boost::tuple<std::vector<double>, dmatrix, std::vector<int> >
count_background(const std::vector<std::string>& sequences)
{
  std::vector<int> character_frequencies(256);
  std::vector<double> background_frequencies(4,0);

  dmatrix background_frequency_matrix(4, 4);
  background_frequency_matrix.fill_with(0);

  int character_count = 0;
  int digram_count = 0;
  static character_to_values<bool> isnuc("ACGT", true);
  
  for (int i=0; i < sequences.size(); ++i) {
    const std::string& line = sequences[i];
    int n = line.length();
    character_count += n;
    digram_count += n - 1;
    for (int j=0; j < n-1; ++j) {
      ++character_frequencies[(unsigned char)line[j]];
      if (isnuc(line[j])) {
	++background_frequencies[to_int(line[j])];
	if (isnuc(line[j+1]))
	  ++background_frequency_matrix(to_int(line[j]),
					to_int(line[j+1]));
      }
    }
    ++character_frequencies[(unsigned char)line[n-1]];
    if (isnuc(line[n-1]))
      ++background_frequencies[to_int(line[n-1])];  
  }

  /*
  // balance the freqs: p(A)==p(T), p(C)==p(G)
  int AT=background_frequencies[0]+background_frequencies[3];
  int CG=background_frequencies[1]+background_frequencies[2];
  if (use_two_strands) {
    character_count *= 2;
    digram_count *= 2;
    background_frequencies[0]=AT;
    background_frequencies[1]=CG;
    background_frequencies[2]=CG;
    background_frequencies[3]=AT;

    // Symmetric counts: AG has same count as CT, etc
    dmatrix temp(background_frequency_matrix.dim());
    for (int i=0; i < 4; ++i)
      for (int j=0; j < 4; ++j)
	temp(i, j) = background_frequency_matrix(3-j, 3-i);
    background_frequency_matrix += temp;
  }
  */
  return boost::make_tuple(background_frequencies, background_frequency_matrix, character_frequencies);
}

