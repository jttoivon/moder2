#define BOOST_TEST_DYN_LINK
#ifdef STAND_ALONE
#define BOOST_TEST_MODULE Main
#endif

#include <boost/test/unit_test.hpp>
//#include <boost/test/included/unit_test.hpp>
#include <boost/test/test_case_template.hpp>
#include <boost/mpl/list.hpp>

#include <boost/random/uniform_01.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random/uniform_int.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/variate_generator.hpp>

#include "common.hpp"
#include "iupac.hpp"
#include "combinatorics.hpp"
#include "test/helper.hpp"

BOOST_AUTO_TEST_SUITE(common)


unsigned int
size_of_hamming_neighbourhood(int len, int hd)
{
  int result = 0;
  for (int r=0; r <= hd; ++r) {
    result += choose(len, r) * pow(3, r);
  }
  return result;
}

    
BOOST_AUTO_TEST_CASE(test_get_n_neighbourhood)
{


  std::string t = "ACGGGAT";
  for (int l=0; l <= t.length(); ++l) {
    std::string s = t.substr(0, l);
    for (int hd=0; hd <= l; ++hd) {
      std::vector<std::string> neighbourhood = get_n_neighbourhood(s, hd, false);
      BOOST_CHECK_EQUAL(neighbourhood.size(), size_of_hamming_neighbourhood(s.length(), hd));
    }
  }

}

BOOST_AUTO_TEST_CASE(test_get_n_neighbourhood_iupac)
{
  random_sequence rand(iupac_chars, 0);
  for (int i=0; i < 10; ++i) {
    std::string s = rand.get_sequence(5);
    //printf("%s\n", s.c_str());
    int hd = 3;
    std::vector<std::string> expanded_neighbourhood = get_n_neighbourhood(s, hd, true);
    std::vector<std::string> unexpanded_neighbourhood = get_n_neighbourhood(s, hd, false);

    int size = 0;
    for (std::string t : unexpanded_neighbourhood) {
      size += number_of_iupac_combinations(t);
    }
    BOOST_CHECK_EQUAL(expanded_neighbourhood.size(), size);

  }
}
//____________________________________________________________________________//

BOOST_AUTO_TEST_SUITE_END()
