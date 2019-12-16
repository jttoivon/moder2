#include <boost/random/uniform_01.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random/uniform_int.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/variate_generator.hpp>

class random_sequence
{
public:
  typedef boost::variate_generator<boost::mt19937, boost::uniform_real<> > real_gen;
  typedef boost::variate_generator<boost::mt19937, boost::uniform_int<> > int_gen;

  random_sequence(const std::string& alphabet_, int random_seed=0)
    : rng(random_seed), die(rng, boost::uniform_real<>(0,1)), alphabet(alphabet_)
  {
    bg.assign(alphabet.length(), 1.0 / alphabet.length());
  }

  std::string
  get_sequence(int L)
  {
    std::string s(L, '-');

    for (int j=0; j < L; ++j)  {
      s[j]=character(bg, die());
    }
    return s;
  }

private:
  char
  character(const std::vector<double>& v, double p)
  {
    assert(p >= 0 && p <= 1);
    double cumulative=0;
    for (int i=0; i<alphabet.size(); ++i) {
      cumulative+=v[i];
      if (p < cumulative)
	return alphabet[i];
    }
    if (p <= 1.0)
      return alphabet.back();
  
    assert(0);
    return 'X';
  }

  boost::mt19937 rng;  
  real_gen die;
  std::vector<double> bg;
  std::string alphabet;					    
};

class random_matrix
{
public:
  typedef boost::variate_generator<boost::mt19937, boost::uniform_int<> > int_gen;

  random_matrix(int min, int max, int random_seed=0)
    : rng(random_seed), die(rng, boost::uniform_int<>(min,max))
  {
  }

  dmatrix
  get_matrix(int rows, int cols)
  {
    dmatrix m(rows, cols);
    for (int r=0; r < rows; ++r) {
      for (int c=0; c < cols; ++c) {
	m(r, c) = die();
      }
    }
    return m;
  }

private:

  boost::mt19937 rng;  
  int_gen die;
};
