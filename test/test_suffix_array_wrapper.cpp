#define BOOST_TEST_DYN_LINK
#ifdef STAND_ALONE
#define BOOST_TEST_MODULE Main
#endif

//#include <boost/test/included/unit_test.hpp>
#include <boost/test/unit_test.hpp>

#include "suffix_array_wrapper.hpp"
#include "iupac.hpp"
#include "common.hpp"

BOOST_AUTO_TEST_SUITE(test_suite_suffix_array_wrapper)

int
search(const std::string& pattern, const std::string& text)
{
  suffix_array sa(text, false);
  if (is_nucleotide_string(pattern))
    return sa.count(pattern);
  else
    return sa.count_iupac(pattern);
}

BOOST_AUTO_TEST_CASE(test_count)
{
  
  BOOST_CHECK_EQUAL(search("A", "CCCGGGTT"), 0);
  BOOST_CHECK_EQUAL(search("A", "CACGGATT"), 2);
  BOOST_CHECK_EQUAL(search("A", "AACCGGTTAA"), 4);
  BOOST_CHECK_EQUAL(search("A", "ACCGGTTTT"), 1);
  BOOST_CHECK_EQUAL(search("A", "CCGGTTTTA"), 1);
  
  BOOST_CHECK_EQUAL(search("T", "CCGGACCGA"), 0);
  
  BOOST_CHECK_EQUAL(search("T", "A"), 0);
  
  BOOST_CHECK_EQUAL(search("T", ""), 0);
  
  BOOST_CHECK_EQUAL(search("TA", "TATATAT"), 3);
  
  BOOST_CHECK_EQUAL(search("A", "AAAAA"), 5);
  
  BOOST_CHECK_EQUAL(search("R", "AGTAGAC"), 5);		    


}

BOOST_AUTO_TEST_CASE(test_count_empty)
{
  std::string text = "AGT";
  suffix_array sa(text, false);
  BOOST_CHECK_EQUAL(sa.count(""), 3);		    
  BOOST_CHECK_EQUAL(sa.count(""), 3);		    
  //BOOST_CHECK_EQUAL(sa.count_iupac(""), 3);		    
}

BOOST_AUTO_TEST_SUITE_END()

