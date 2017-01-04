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
