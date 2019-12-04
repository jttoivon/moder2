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

#include "probabilities.hpp"








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
  double p = compute_bernoulli_probability<double>(s, bg);
  BOOST_CHECK_EQUAL(p, 0.0);

  std::string s2(40, 'A');
  double p2 = compute_bernoulli_probability<double>(s2, bg);
  BOOST_CHECK_NE(p2, 0.0);
  
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
