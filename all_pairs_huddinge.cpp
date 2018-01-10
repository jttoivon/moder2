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

///////////////////////////////
//
// (c) Kimmo Palin
//
// email: kimmo.palin@helsinki.fi
//
///////////////////////////////

#define TIMING

#ifndef PACKAGE_VERSION
#define PACKAGE_VERSION "unknown"
#endif

#include "iupac.hpp"
#include "matrix.hpp"
#include "vectors.hpp"
#include "timing.hpp"
#include "matrix_tools.hpp"
#include "common.hpp"
#include "probabilities.hpp"
#include "parameters.hpp"
#include "orientation.hpp"
#include "kmer_tools.hpp"
#include "multinomial_helper.hpp"
#include "suffix_array_wrapper.hpp"
#include "huddinge.hpp"

#include <sys/stat.h>
#include <cstdio>
#include <cmath>
#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif
#include <fenv.h>

#include <libgen.h>   // for dirname and basename
#include <omp.h>

#include <string>
#include <iostream>
#include <fstream>
#include <numeric>
#include <queue>
#include <cfloat>

#include <boost/algorithm/string.hpp>
#include <boost/foreach.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/program_options.hpp>
#include <boost/numeric/ublas/symmetric.hpp>

namespace po = boost::program_options;

typedef long double FloatType;


typedef std::vector<boost::tuple<int, int, int, int, std::string> > overlapping_dimer_cases_t;
typedef std::vector<boost::tuple<int, int, int, int> > spaced_dimer_cases_t;

typedef boost::multi_array<dmatrix, 2> cob_of_matrices;




bool use_palindromic_correction=false;
bool use_multimer=true;
bool use_meme_init=false;
bool get_full_flanks=false;  // make motifs for each model that have width 2*L-motif_width

int max_iter = 50;  // was 300
int minimum_distance_for_learning = 4;
int global_dmax = 10;
//int global_max_dist_for_deviation = -1;
//int global_max_dist_for_deviation = 4;
int global_max_dist_for_deviation = 1000;
//int min_flank = 3;  // On each side of overlapping area of a dimer at least min_flank number of positions
                    // must be reserved. That is the minimum non-overlapping part on each side.
double ic_threshold = 0.40;
double learning_fraction = 0.02;
//int hamming_threshold = 4;

int character_count=0;
int digram_count=0;

std::vector<double> background_frequencies(4);
std::vector<double> background_probabilities(4);
dmatrix background_frequency_matrix(4,4);   // for Markov background
dmatrix background_probability_matrix(4,4); // for Markov background


prior<double> pseudo_counts;

double cob_cutoff = 0.001;  // if an element in a cob table is smaller than this constant,
                           // the element is excluded. Note! This works for spaced dimers as well

bool adjust_seeds = true;
bool use_multinomial=true;
bool local_debug = true;
bool extra_debug = false;   // Even more printing
bool allow_extension = false;
bool use_dimers = true;
bool seeds_given = false;
bool no_unique = false;
bool use_output = false; // whether to write model parameters to files
bool maximize_overlapping_seeds=true;
bool require_directional_seed = false;
bool avoid_palindromes = false; // If tf1==tf2, the orientation is HH or TT and the PPM \tau_ht1,ht2,o,d, the
                                // probability of a sequence is the same in both directions.
                                // If we recognize this situation, we can try to escape from this palindromicity
                                // by using temporarily only a single strand.
int hamming_distance_overlapping_seeds_N=0;
int hamming_distance_overlapping_seeds_OR=2;
int hamming_radius = 1; // For learning the model from sequences in the hamming neighbourhood of the seed

std::string unbound; // Filename where to store the sequences that had greatest probability under background model
std::vector<std::string> names;
std::string outputfile = "huddinge.dists";


int main(int argc, char* argv[])
{
  TIME_START(t);
  WALL_TIME_START(t2);
#ifdef FE_NOMASK_ENV
  feenableexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW | FE_UNDERFLOW);   // These exceptions cause trap to occur
#endif

  if (argc > 1)
    print_command_line(argc, argv);

  use_two_strands = true;
  int unique_param = -1; // either -1, 0, 1, ...
  using namespace boost;
  using std::string;
  using std::cout;

  std::vector<std::string> sequences;
  std::vector<bool> keep_monomer_fixed;
  

  int number_of_threads=1;
#if defined(_OPENMP)
  {
    char* p = getenv("OMP_NUM_THREADS"); 
    if (p)
      number_of_threads = atoi(p);
  }
#endif

  
  // Declare the supported options.
  po::options_description hidden("Allowed options");
  hidden.add_options()
    ("sequencefile",                           "Name of the sequence file containing sequences")
    ;

  po::options_description nonhidden("Allowed options");
  nonhidden.add_options()
    ("help",                           "Produce help message")
    ("output", m("Write model parameters to files", use_output).c_str())
    ("names", po::value<string>(),     "Names for the monomer models. Comma separated list. Default: TF0,TF1, ...")
    ("outputfile", po::value<string>(), m("Output file for distances", outputfile).c_str())
    ("number-of-threads", po::value<int>(), m("Number of threads", number_of_threads).c_str())
    ("single-strand", m("Only consider binding sites in forward strand", not use_two_strands).c_str())
    ("unique", po::value<std::string>(), "Uniqueness of sequences. Either off, unique, 1, 2, 3, ..., default: off")
    ("quiet", m("Don't print intermediate results", not local_debug).c_str())
    ;

  po::options_description desc("Allowed options");
  desc.add(hidden).add(nonhidden);

  po::positional_options_description popt;
  popt.add("sequencefile", 1);

  
  string seqsfile;

  string prior_parameter;
  po::variables_map vm;
  // int overlapping_dimer_p = 0; // number of overlapping models
  // int spaced_dimer_p = 0; // number of spaced models
  boost::multi_array<std::string, 2> dimer_seeds;

   ////////////////////////////////
  // parse command line parameters
  ////////////////////////////////

  try {
    po::store(po::command_line_parser(argc, argv).
	      options(desc).positional(popt).run(), vm);
    po::notify(vm);

    // are all positional parameters present
    bool all_positional = vm.count("sequencefile");

    if (vm.count("help") || not all_positional) {
      cout << basename(argv[0]) << " version " << PACKAGE_VERSION << std::endl;
      cout << "Usage:\n" 
	   << basename(argv[0]) << " [options] sequencefile\n";
      cout << nonhidden << "\n";
      return 1;
    }

    if (vm.count("single-strand")) 
      use_two_strands = false;

    if (vm.count("unique")) {
      std::string arg = vm["unique"].as< string >();
      if (arg == "off")
	unique_param = -1;
      else if (arg == "unique")
	unique_param = 0;
      else {
	unique_param = atoi(arg);
	error(unique_param <= 0, "Argument to --unique must be one of the following: off, unique, 1, 2, 3, ...");
      }
    }
    

    if (vm.count("quiet")) 
      local_debug = false;

    if (vm.count("flanks")) 
      get_full_flanks = true;

    if (vm.count("output")) 
      use_output = true;

    if (vm.count("outputfile")) {
      outputfile = vm["outputfile"].as< string >();
    } else {
      error(vm.count("outputfile"), "Must give output file.");
    }


#if defined(_OPENMP)
    if (vm.count("number-of-threads")) {
      number_of_threads = vm["number-of-threads"].as< int >();
    }
    omp_set_num_threads(number_of_threads);
    printf("Using %i openmp threads.\n", number_of_threads);
#endif

    seqsfile   = vm["sequencefile"].as< string >();


   
  }
  catch (std::exception& e) {
    std::cerr << e.what() << std::endl;
    std::cerr << desc << "\n";
    exit(1);
  }

  
  int lines, bad_lines;


  boost::tie(lines, bad_lines) = read_sequences(seqsfile, sequences);
  check_data(sequences);
  printf("Read %zu good lines from file %s\n", sequences.size(), seqsfile.c_str());
  printf("Discarded %i bad lines\n", bad_lines);
  if (unique_param >= 0) {
    TIME_START(s1);
    sequences = remove_duplicate_reads_faster(sequences, unique_param);
    TIME_PRINT("\tRemoving Hamming duplicates took %.2f seconds.\n", s1);
  }
  printf("Using %zu sequences\n", sequences.size());

  const unsigned int N= sequences.size();

  //boost::numeric::ublas::symmetric_matrix<int> huddinge_dists(N,N);
  unsigned char *huddinge_dists;
  size_t dist_size = (size_t)(N*((N-1.0)/2.0));

  printf("%lu - %g = %g == 0??\n",dist_size,(double)N*(N-1.0)/2.0, (double)dist_size - (double)N*(N-1.0)/2.0);
  assert(fabs( (double)dist_size - (double)N*(N-1.0)/2.0) <0.01) ;

  huddinge_dists = (unsigned char*)malloc(dist_size);

  if(huddinge_dists==NULL) {
  std::cerr << "Couldn't allocate "<<dist_size <<" bytes ("<<(dist_size/(1024.*1024.0*1024.0))
	    <<" GB) for pairwise distances of "<< N <<" sequences"<<std::endl;
  exit(2);
  }

#pragma omp parallel for schedule(guided)
  for(long unsigned int i = 1; i<N; i++) {
    long unsigned int idx;
    if(i%2==0) {  // idx=i*(i-1)/2 + j without huge multiplication
      idx = (i>>1)*(i-1);
    } else {
      idx = i*((i-1)>>1);
    }

    if((i%1000)==0 ) {
      printf("Sequence %lu '%s' at idx %lu. Still %lu to go. Error %g.\n",i,sequences[i].c_str(),
	     idx, dist_size -idx,
     	     (double)idx-(i*(i-1.0)/2.0));
    }


    for(long unsigned int j = 0; j<i; j++,idx++) {
      int d = huddinge_distance(sequences[i],sequences[j]);
        huddinge_dists[idx] = (unsigned char)d;
    }
  }

  printf("Writing %g distances\n",(double)N*(N-1.0)/2.0);
  // Output format:
  // Number N of sequences
  // List of N sequences, each on it's own row
  // N*(N-1)/2 rows with format i<tab>j<tab>hudding distance between sequences i and j
  std::ofstream fout;
  fout.open(outputfile);

  fout << N <<'\n';
  for(int i = 0; i<N; i++) {
    fout << sequences[i] <<'\n';
  }
  fout.write((const char*)huddinge_dists,dist_size);
  fout.close();
  free(huddinge_dists);
  
  TIME_PRINT("Whole program took %.2f seconds cpu-time\n", t);
  WALL_TIME_PRINT("Whole program took %.2f seconds wall-time\n", t2);

  return 0;
}
