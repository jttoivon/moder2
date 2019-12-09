#define BOOST_TEST_DYN_LINK
#ifdef STAND_ALONE
#define BOOST_TEST_MODULE Main
#endif

#include <boost/test/unit_test.hpp>
//#include <boost/test/included/unit_test.hpp>
#include <boost/test/test_case_template.hpp>
#include <boost/mpl/list.hpp>

#include <random>
#include <cassert>
#include <vector>
#include <ostream>
#include <fstream>
#include <libgen.h>
#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif
#include <fenv.h>

#include <boost/random/uniform_01.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random/uniform_int.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/variate_generator.hpp>

#include "probabilities.hpp"

typedef boost::variate_generator<boost::mt19937, boost::uniform_real<> > real_gen;
typedef boost::variate_generator<boost::mt19937, boost::uniform_int<> > int_gen;

char
character(const std::vector<double>& v, double p)
{
  char nucleotide_array[] = "ACGT";
  assert(p >= 0 && p <= 1);
  double cumulative=0;
  for (int i=0; i<4; ++i) {
    cumulative+=v[i];
    if (p < cumulative)
      return nucleotide_array[i];
  }
  if (p <= 1.0)
    return nucleotide_array[3];
  
  assert(0);
  return 'X';
}

std::string
random_sequence(const std::vector<double>& v, int L, real_gen& die)
{
  std::string s(L, '-');

  for (int j=0; j < L; ++j)  {    // left from motif
    s[j]=character(v, die());
  }
  return s;
}


// see https://stackoverflow.com/questions/10976130/boost-check-equal-with-pairint-int-and-custom-operator
namespace boost
{
    namespace test_tools
    {
#if BOOST_VERSION >= 105900
      namespace tt_detail
      {
#endif
	template<typename T,typename U>
	struct print_log_value<std::pair<T, U> >
	{
	  void operator()(std::ostream& os, std::pair<T, U> const& pr)
	  {
	    os << "<" << std::get<0>(pr) << "," << std::get<1>(pr) << ">";
	  }
	};
#if BOOST_VERSION >= 105900
      }
#endif 
    }
}

// struct F {
//   F() : hoxb13(get_model_filename())       { BOOST_TEST_MESSAGE( "setup fixture" ); }
//   ~F()         { BOOST_TEST_MESSAGE( "teardown fixture" ); }

//   dinuc_model<double> hoxb13;
  
// private:

//   std::string
//   get_model_filename() const
//   {
//     std::string temp(boost::unit_test::framework::master_test_suite().argv[0]);
//     std::string dir = dirname((char*)temp.c_str());
//     //printf("argv[0]: %s\n", temp.c_str());
//     //printf("dirname(argv[0]): %s\n", dir.c_str());
//     return dir + "/HOXB13-CYMRTAAAA.adm";
//   }
  
// };

//____________________________________________________________________________//


BOOST_AUTO_TEST_SUITE( probabilities)




BOOST_AUTO_TEST_CASE(test_compute_bernoulli_probability_code)
{
  int k = 5;
  int number_of_codes = pow(4, k);
  std::vector<double> bg(4, 0.25);
  double total = 0.0;
  for (int code=0; code < number_of_codes; ++code) {
  
    double prob = compute_bernoulli_probability<double>(code, k, bg);
    BOOST_CHECK_GE(prob, 0.0);
    BOOST_CHECK_LE(prob, 1.0);
    total += prob;
  }
  BOOST_CHECK_EQUAL(total, 1.0);
}

BOOST_AUTO_TEST_CASE(test_compute_bernoulli_probability_overflow)
{
  std::vector<double> bg(4, 0.25);
  std::string s(4000, 'A');
  feclearexcept(FE_ALL_EXCEPT);
  double p = compute_bernoulli_probability<double>(s, bg);
  BOOST_CHECK_EQUAL(p, 0.0);
  int retval = fetestexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW | FE_UNDERFLOW);
  BOOST_CHECK(retval & FE_UNDERFLOW);

  feclearexcept(FE_ALL_EXCEPT);
  std::string s2(40, 'A');
  double p2 = compute_bernoulli_probability<double>(s2, bg);
  BOOST_CHECK_NE(p2, 0.0);
  retval = fetestexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW | FE_UNDERFLOW);
  BOOST_CHECK(not (retval & FE_UNDERFLOW));
}

BOOST_AUTO_TEST_CASE(test_linear_and_log_compute_bernoulli_probability)
{
  int random_seed=0;
  boost::mt19937 rng(random_seed);  
  boost::uniform_real<> dist(0,1); 
  real_gen die(rng, dist);   // to choose the strand
  std::vector<std::vector<double> > bgs = { {0.25, 0.25, 0.25, 0.25},
					    {0.29, 0.37, 0.17, 0.17} };
  for (auto bg : bgs) {
    std::vector<double> log_bg;
    for (double x : bg)
      log_bg.push_back(log2l(x));
    int L=40;
    for (int i=0; i < 100; ++i) {
      std::string s = random_sequence(bg, L, die);
      BOOST_CHECK_CLOSE(compute_bernoulli_probability<double>(s, bg),
			exp2l(compute_log_bernoulli_probability<double>(s, log_bg)),
			pow(10, -12)*100); // this is a percentage
    }
  }
}

BOOST_AUTO_TEST_CASE(test_log_sum)
{
  std::priority_queue<FloatType, std::vector<FloatType>, std::greater<FloatType> > queue;
  std::vector<double> probabilities {0.5, 0.2, 0.05, 0.005, 0.001};
  double sum = 0.0;
  for (double p : probabilities) {
    sum += p;
    queue.push(log2(p));
    log_sum(queue);
    BOOST_CHECK_EQUAL(exp2(log_sum(queue)), sum);
  }
  
  while(not queue.empty())
        queue.pop();

  
  std::vector<double> probabilities2;
  for (int i=0; i < 100; ++i)
    probabilities2.push_back(1.0/200);
  sum = 0.0;
  for (double p : probabilities2) {
    sum += p;
    queue.push(log2(p));
    log_sum(queue);
    BOOST_CHECK_CLOSE(exp2(log_sum(queue)), sum, 1.0e-10);
  }
}


//____________________________________________________________________________//

BOOST_AUTO_TEST_SUITE_END()
