/*

    MODER is a program to learn DNA binding motifs from SELEX datasets.
    Copyright (C) 2016, 2017  Jarkko Toivonen

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
#include <vector>

int
randrange(int n);

std::vector<int>
random_sample(int n, int k);

unsigned long
factorial(int n);


unsigned long
falling_factorial(int n, int k);


unsigned long
choose(int n, int k);


// n elements in k bins
unsigned long
number_of_distributions(int n, int k);


// number of sequences of length n with alphabet {a, c, g, t} 
// when we can use distribution v. Sum(v) == n.
// In other words: number of sequences with character distribution v
int
number_of_sequences(const std::vector<int>& v);


double
multinomial_coeff(const std::vector<int>& g);

double
boost_factorial(int s);

long double
ln_binomial(int k, int n, long double p);
