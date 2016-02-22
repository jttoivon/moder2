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
#ifndef PROBABILITIES_HPP
#define PROBABILITIES_HPP

#include "matrix.hpp"
#include "data.hpp"

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
max_matrix_logodds_score(const dmatrix& m);

std::string
string_giving_max_score(const dmatrix& m);

double
max_matrix_normal_score(const dmatrix& m);

matrix<double>
count_positional_background(const std::vector<std::string>& sequences);


double
compute_background_probability(const std::string& line, const std::vector<double>& q,
			       const matrix<double>& q2 = dmatrix());

// wrapper for the one below

double
compute_probability(const std::string& line, int pos, int direction,
		    const matrix<double>& m, 
		    const std::vector<double>& q,
		    const matrix<double>& q2 = dmatrix());

double
compute_probability(const std::string& line_orig, const std::string& line_orig_rev,
		    int pos, int direction,
		    const matrix<double>& m, 
		    const std::vector<double>& q,
		    const matrix<double>& q2 = dmatrix());



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


double
compute_bernoulli_probability(const std::string& line, const std::vector<double>& q);

double
compute_bernoulli_probability(big_int code, int k, const std::vector<double>& q);



boost::tuple<std::vector<double>, dmatrix>
count_background(const std::vector<std::string>& sequences);

#endif // PROBABILITIES_HPP
