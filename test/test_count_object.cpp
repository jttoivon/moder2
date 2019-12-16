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
#include "count_object.hpp"
#include "iupac.hpp"
#include "combinatorics.hpp"
#include "test/helper.hpp"

int hamming_radius=3;

double
lowest_dinucleotide_count(const std::string& w, const dmatrix& all_counts);

std::vector<double>
lowest_dinucleotide_count_iupac(const std::string& w, const dmatrix& all_counts);

std::tuple<std::vector<double>, std::vector<double> >
compute_bias_and_low_counts(const std::string& s, int t, int j, const dmatrix& all_counts, const dmatrix& corrected);

std::tuple<std::vector<double>, std::vector<double> >
compute_bias_and_low_counts_iupac(const std::string& s, int t, int j, const dmatrix& all_counts, const dmatrix& corrected);

bool extra_debug = false;

// see https://stackoverflow.com/questions/10976130/boost-check-equal-with-pairint-int-and-custom-operator
namespace boost
{
    namespace test_tools
    {
#if BOOST_VERSION >= 105900
      namespace tt_detail
      {
#endif
	template<typename T>
	struct print_log_value<matrix<T> >
	{
	  void operator()(std::ostream& os, matrix<T> const& m)
	  {
	    write_matrix(stdout, m, "");
	  }
	};
#if BOOST_VERSION >= 105900
      }
#endif 
    }
}

void
normalize_adm(dmatrix& adm)
{
  int k = adm.get_columns();
  for (int c=0; c < k; ++c) {
    for (int a=0; a < 4; ++a) {
      double s = 0.0;
      for (int b=0; b < 4; ++b) {
	s += adm(4*a + b, c);
      }
      if (s > 0.0) {
	for (int b=0; b < 4; ++b) {
	  adm(4*a + b, c) /= s;
	}
      }
    }
  }
}

bool
almost_equal_matrices(const dmatrix& a, const dmatrix& b)
{
  int rows = a.get_rows();
  int cols = a.get_columns();
  for (int row=0; row<rows; ++row) {
    for (int col=0; col<cols; ++col) {
      if (fabs(a(row, col) - b(row, col)) > 0.0000000001)
	return false;
    }
  }
  return true;
}

BOOST_AUTO_TEST_SUITE(test_count_object)

BOOST_AUTO_TEST_CASE(test_lowest_dinucleotide_count)
{
  std::string nucs = "ACGT";
  //  random_sequence rand(iupac_chars, 0);
  random_sequence rand(nucs, 0);
  random_matrix rand2(0, 1000, 0);
  for (int i=0; i < 10; ++i) {
    std::string s = rand.get_sequence(10);
    //printf("%s\n", s.c_str());
    for (int j=0; j < 10; ++j) {
      dmatrix all_counts = rand2.get_matrix(16, 15);
      //write_matrix(stdout, all_counts, "", "%.0f");
      std::vector<double> v = lowest_dinucleotide_count_iupac(s, all_counts);
      double result = *std::max_element(v.begin(), v.end());
      BOOST_CHECK_EQUAL(lowest_dinucleotide_count(s, all_counts),
			result);
    }
    // int hd = 3;
    // std::vector<std::string> expanded_neighbourhood = get_n_neighbourhood(s, hd, true);
    // std::vector<std::string> unexpanded_neighbourhood = get_n_neighbourhood(s, hd, false);

    // int size = 0;
    // for (std::string t : unexpanded_neighbourhood) {
    //   size += number_of_iupac_combinations(t);
    // }
    // BOOST_CHECK_EQUAL(expanded_neighbourhood.size(), size);

  }
}


BOOST_AUTO_TEST_CASE(test_lowest_dinucleotide_count2)
{
  std::string nucs = "ACGT";
  random_sequence rand(iupac_chars, 0);
  //random_sequence rand(nucs, 0);
  random_matrix rand2(0, 1000, 0);
  for (int j=0; j < 10; ++j) {
    dmatrix all_counts = rand2.get_matrix(16, 15);
    for (int i=0; i < 10; ++i) {
      std::string s = rand.get_sequence(10);
      //printf("%s\n\n", s.c_str());
      //int hd = 3;
      //std::vector<std::string> expanded_neighbourhood = get_n_neighbourhood(s, hd, true);
      //std::vector<std::string> unexpanded_neighbourhood = get_n_neighbourhood(s, hd, false);
      //write_matrix(stdout, all_counts, "", "%.0f");
      double count1 = std::numeric_limits<double>::lowest();
      sequences_of_iupac_string ss(s);
      int number_of_seqs = ss.number_of_strings();
      for (int a=0; a < number_of_seqs; ++a) {
	std::string t = ss.next();
	//printf("%s\n", t.c_str());
	count1 = std::max(count1, lowest_dinucleotide_count(t, all_counts));
      }
      std::vector<double> v = lowest_dinucleotide_count_iupac(s, all_counts);
      double count2 = *std::max_element(v.begin(), v.end());
      BOOST_CHECK_EQUAL(count1, count2);
    }

    // int size = 0;
    // for (std::string t : unexpanded_neighbourhood) {
    //   size += number_of_iupac_combinations(t);
    // }
    // BOOST_CHECK_EQUAL(expanded_neighbourhood.size(), size);

  }
}

BOOST_AUTO_TEST_CASE(test_compute_bias_and_low_counts)
{
  std::string nucs = "ACGT";
  random_sequence rand(iupac_chars, 0);
  random_matrix rand2(0, 1000, 0);
  //random_sequence rand(nucs, 0);
  int k = 10;
  int hamming_radius = 3;
  for (int i=0; i < 10; ++i) {
    std::string seed = rand.get_sequence(k);
    for (int i2=0; i2 < 10; ++i2) {
    
      dmatrix all_counts = rand2.get_matrix(16, k);
  
      //  write_matrix(stdout, motif, to_string("%s dinucleotide motif matrix counts:\n", name.c_str()), "%.0f");
      dmatrix corrected = all_counts;
      normalize_adm(corrected);
      for (int j=k-1; j >= 0; --j) {
    
	int tmax = std::min(hamming_radius-1, k-j-1);          // Compute the correction factors tau[b][t]
	if (j == k-1) {

	}
	else {
	  for (int t=0; t <= tmax; ++t) {
	    int suffix_len = k-j-1;
	    std::vector<double> my_tau;
	    std::vector<double> my_low_counts;
	    std::vector<double> my_tau2;
	    std::vector<double> my_low_counts2;
	    std::string sub = seed.substr(j+1, suffix_len);
	    std::tie(my_tau, my_low_counts)   = compute_bias_and_low_counts(sub, t, j, all_counts, corrected);
	    std::tie(my_tau2, my_low_counts2) = compute_bias_and_low_counts_iupac(sub, t, j, all_counts, corrected);
	    for (int b=0; b < 4; ++b) {
	      BOOST_CHECK_CLOSE(my_tau[b], my_tau2[b], 1.0e-10);
	      BOOST_CHECK_EQUAL(my_low_counts[b], my_low_counts2[b]);
	      //BOOST_CHECK_MESSAGE(my_low_counts[b] == my_low_counts2[b], to_string("j=%i, t=%i, b=%i, sub=%s\n", j, t, b, sub.c_str()));
	    }
	  }
	}
      }
    }
  }
}

BOOST_AUTO_TEST_CASE(test_create_pseudo_count_tables)
{
  std::string nucs = "ACGT";
  int k=10;
  std::vector<double> bg(4, 0.25);
  count_object co(adm_fixed, k);
  random_sequence rand(iupac_chars, 0);
  //random_sequence rand(nucs, 0);
  for (int i=0; i < 10; ++i) {
    std::string seed = rand.get_sequence(k);
    std::vector<dmatrix> pc1 = co.create_pseudo_count_tables<myuint128>(seed, bg);
    std::vector<dmatrix> pc2 = co.create_pseudo_count_tables_slow<myuint128>(seed, bg);
    BOOST_CHECK_EQUAL(pc1.size(), pc2.size());
    for (int r=0; r<pc1.size(); ++r) {
      //BOOST_CHECK_EQUAL(pc1[r], pc2[r]);
      BOOST_CHECK(almost_equal_matrices(pc1[r], pc2[r]));
    }
  }
}


//____________________________________________________________________________//

BOOST_AUTO_TEST_SUITE_END()
