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
#include "combinatorics.hpp"
#include "common.hpp"
#include <cassert>
#include <cmath>


#include <boost/math/special_functions/factorials.hpp>

#include <boost/random/linear_congruential.hpp>
#include <boost/random/uniform_int.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random/variate_generator.hpp>


#include <boost/random/uniform_int.hpp>
#include <boost/random/mersenne_twister.hpp>

boost::mt19937 random_number_generator;         // produces randomness out of thin air

// returns a random number from range 0...n-1
int
randrange(int n)
{
  boost::uniform_int<> d(0,n-1); // distribution that maps to 0..n-1
  int x = d(random_number_generator);   
  return x;
}

// returns a vector of k integers from range 0...n-1
std::vector<int>
random_sample(int n, int k)
{
  std::vector<int> result(n);
  for (int i=0; i < n; ++i)
    result[i]=i;

  for(int i = 0; i < k-1; i++) {
    int c = randrange(n-i);
    std::swap(result[i], result[i+c]);
  }

  result.resize(k);
  assert(result.size() == k);

  return result;
}

unsigned long
factorial(int n)
{
  assert(n >= 0);
  assert(n <= 20);  // otherwise won't fit in 64 bits
  unsigned long f = 1;
  int i;
  for (i=2; i <= n; ++i)
    f *= i;
  return f;
}

unsigned long
falling_factorial(int n, int k)
{
  assert(n >= 0);
  assert(k <= n);

  unsigned long f = 1;
  for (int i=n, l=0; l < k; --i, ++l)
    f *= i;
  return f;
}

unsigned long
choose(int n, int k)
{
  //  return factorial(n) / (factorial(n-k)*factorial(k));

  return falling_factorial(n, k) / factorial(k);
}

// n elements in k bins
unsigned long
number_of_distributions(int n, int k)
{
  return choose(n+k-1, k-1);
}

// number of sequences of length n with alphabet {a, c, g, t} 
// when we can use distribution v. Sum(v) == n.
// In other words: number of sequences with character distribution v
int
number_of_sequences(const std::vector<int>& v)
{
  int n = sum(v);
  return choose(n, v[0]) * choose(n-v[0], v[1]) * choose(n-v[0]-v[1], v[2]);
}

double
multinomial_coeff(const std::vector<int>& g)
{
  assert(g.size() == 4);

  int s = sum(g);
  double x = 1.0;
  for (int i=0; i < g.size(); ++i)
    x *= boost::math::factorial<double>(g[i]);

  return boost::math::factorial<double>(s) / x;
}

double
boost_factorial(int s)
{
  return boost::math::factorial<double>(s);
}

// Uses the Stirling approximation
long double
ln_factorial(int k)
{
  if (k == 0)
    return 0;

  return (k+0.5)*logl(k) - k + logl(sqrt(2*M_PI)) + logl(1+1.0/(12*k));
}

long double
ln_factorial2(int k)
{
  if (k == 0)
    return 0;
  return lgammal(k+1);
}

// Returns log(factorial(n)/(factorial(k)*factorial(n-k)) * p^k * (1-p)^{n-k})
long double
ln_binomial(int k, int n, long double p)
{
  return ln_factorial2(n) - ln_factorial2(k) - ln_factorial2(n-k) + k*logl(p) + (n-k)*logl(1-p);
}
