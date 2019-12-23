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

#include "dinucleotide.hpp"
#include "test/helper.hpp"


dinuc_model<double>
adm_from_background(const std::vector<double>& bg, int k)
{
  assert(k >= 1);
  dmatrix result(16, k);
  dmatrix matrix_bg(4,1);
  for (int b=0; b < 4; ++b)
    matrix_bg(b, 0) = bg[b];
      
  for (int j=0; j < k; ++j) {
    int amax = j==0 ? 1 : 4;
    for (int a=0; a < amax; ++a) {
      result.inject(matrix_bg, 4*a, j);
    }
  }

  dinuc_model<double> adm;
  adm.init(result);
  return adm;
}

dinuc_model<double>
random_adm(int k)
{
  assert(k >= 1);
  dmatrix result(16, k);

  int seed = 0;
  //std::random_device rd;  //Will be used to obtain a seed for the random number engine
  //  std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
  std::mt19937 gen(seed); //Standard mersenne_twister_engine seeded with rd()
  std::uniform_real_distribution<> dis(0.0, 1.0);
  dmatrix m(4,1);
  
  for (int j=0; j < k; ++j) {
    int amax = j==0 ? 1 : 4;
    for (int a=0; a < amax; ++a) {
      for (int b=0; b < 4; ++b)
	m(b, 0) = dis(gen);
      normalize_matrix_columns(m);
      result.inject(m, 4*a, j);
    }
  }

  dinuc_model<double> adm(result);
  return adm;
}

// Does not work
// template <typename A, typename B>
// std::ostream&
// operator<<(std::ostream& str, std::pair<A, B> p)
// {
//   str << "<" << p.first << "," << p.second << ">";

//   return str;
// }

//typedef double FloatType;

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

struct F {
  F() : hoxb13(get_model_filename())       { BOOST_TEST_MESSAGE( "setup fixture" ); }
  ~F()         { BOOST_TEST_MESSAGE( "teardown fixture" ); }

  dinuc_model<double> hoxb13;
  
private:

  std::string
  get_model_filename() const
  {
    std::string temp(boost::unit_test::framework::master_test_suite().argv[0]);
    std::string dir = dirname((char*)temp.c_str());
    //printf("argv[0]: %s\n", temp.c_str());
    //printf("dirname(argv[0]): %s\n", dir.c_str());
    return dir + "/HOXB13-CYMRTAAAA.adm";
  }
  
};

//____________________________________________________________________________//


BOOST_FIXTURE_TEST_SUITE( dinucleotide, F )

BOOST_AUTO_TEST_CASE(test_sub_adm)
{
  int k = 9;
  dinuc_model<double> adm = hoxb13;
  BOOST_CHECK_EQUAL(adm.dim(), std::make_pair(16, k));

  dinuc_model<double> another_adm = sub_adm(adm, 0, k);
  BOOST_CHECK_EQUAL(adm, another_adm);
}


BOOST_AUTO_TEST_CASE(test_left_extend_adm)
{
  int k = 9;
  dinuc_model<double> adm = hoxb13;
  std::vector<double> bg(4, 0.25);
  
  for (int extension=0; extension < 10; ++extension) {

    dinuc_model<double> another_adm = left_extend_adm(adm, bg, extension);
    BOOST_CHECK_EQUAL(another_adm.dim(), std::make_pair(16, k+extension));
    BOOST_CHECK_EQUAL(adm, sub_adm(another_adm, extension, k));
  }
}

BOOST_AUTO_TEST_CASE(test_right_extend_adm)
{
  int k = 9;
  dinuc_model<double> adm = hoxb13;
  std::vector<double> bg(4, 0.25);
  
  for (int extension=0; extension < 10; ++extension) {

    dinuc_model<double> another_adm = right_extend_adm(adm, bg, extension);
    BOOST_CHECK_EQUAL(another_adm.dim(), std::make_pair(16, k+extension));
    BOOST_CHECK_EQUAL(adm, sub_adm(another_adm, 0, k));
  }
}


BOOST_AUTO_TEST_CASE(test_reverse_complement)
{
  dinuc_model<double> adm = hoxb13;

  dinuc_model<double> another_adm = reverse_complement(adm);
  
  BOOST_CHECK_EQUAL(adm, reverse_complement(another_adm));

}

BOOST_AUTO_TEST_CASE(test_combination)
{
  int k = 9;
  int extension = 5;
  std::vector<double> bg(4, 0.25);
  
  dinuc_model<double> adm = hoxb13;

  dinuc_model<double> another_adm = left_extend_adm(adm, bg, extension);
  another_adm = reverse_complement(another_adm);
  another_adm = sub_adm(another_adm, 0, k);
  another_adm = reverse_complement(another_adm);
  BOOST_CHECK_EQUAL(adm, another_adm);

}


BOOST_AUTO_TEST_CASE(test_force_adms_equal)
{
  int k = 9;
  std::vector<double> bg(4, 0.25);
  
  dinuc_model<double> adm = hoxb13;

  dinuc_model<double> background_adm = adm_from_background(bg, k);
  dinuc_model<double> another_adm = force_adms_equal(adm, background_adm);
  BOOST_CHECK_EQUAL(adm, another_adm);
  
}

BOOST_AUTO_TEST_CASE(test_representation)
{
  
  dinuc_model<double> adm = hoxb13;
  int k = adm.get_length();
  
  dmatrix rep = adm.representation();
  BOOST_CHECK_EQUAL(rep.dim(), std::make_pair(16, k));
  for (int a=4; a < 16; ++a)
    BOOST_CHECK_EQUAL(rep(a,0), 0.0);
  dinuc_model<double> adm2 = dinuc_model<double>(rep);
  BOOST_CHECK_EQUAL(adm, adm2);
}

BOOST_AUTO_TEST_CASE(test_probability)
{
  
  dinuc_model<double> adm = hoxb13;
  dinuc_model<double> rev = reverse_complement(adm);
  int k = adm.get_length();

  for (int code=0; code < pow(4, k); ++code) {
    std::string s = number_to_dna(code, k);
    double p1 = rev.probability(s);
    double p2 = adm.probability(reverse_complement(s));
    BOOST_CHECK_CLOSE(p1, p2, 0.0001);
  }
}

BOOST_AUTO_TEST_CASE(test_cut_and_join)
{
  
  dinuc_model<double> adm = hoxb13;
  int k = adm.get_length();

  for (int j=1; j < k; ++j) {
    dinuc_model<double> left = sub_adm(adm, 0, j);
    dinuc_model<double> right = sub_adm(adm, j-1, k-j+1);
    BOOST_CHECK_EQUAL(left.get_length(), j);
    BOOST_CHECK_EQUAL(right.get_length(), k-j+1);
    dinuc_model<double> join = join_adms(left, right);
    BOOST_CHECK_EQUAL(adm, join);
  }
  
}

BOOST_AUTO_TEST_CASE(test_dinucleotide_counts)
{
  bool use_rna = false;
  use_two_strands = false;
  std::string nucs ="ACGT";
  random_sequence rand(nucs, 0);
  int n = 1000;
  int hd=2;
  model_type model_type = adm_unfixed;
  std::string seed = "ACGTAAG";
  int k = seed.length();
  std::vector<std::string> sequences;
  for (int i=0; i < n; ++i)
    sequences.push_back(rand.get_sequence(40));
  
  std::string str=join(sequences, '#');
  if (use_two_strands) {
    str.append("#");
    str += join_rev(sequences, '#');
  }
  suffix_array sa(str, use_rna);
  
  std::vector<dmatrix> v;
  std::vector<dmatrix> v2;
  v = dinucleotide_counts_scan_better<myuint128>(seed, sequences, hd, model_type);
  v2 = dinucleotide_counts_suffix_array(seed, sequences, sa, hd, model_type);
				      
  for (int row=0; row < 16; ++row) {
    for (int col=0; col < k; ++col) {
      BOOST_CHECK_EQUAL(v[0](row, col), v2[0](row, col));
    }
  }
    
}

//____________________________________________________________________________//

BOOST_AUTO_TEST_SUITE_END()
