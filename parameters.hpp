#ifndef PARAMETERS_HPP
#define PARAMETERS_HPP
extern bool use_markov_background;
extern bool use_positional_background;
extern bool use_two_strands;
extern bool use_pseudo_counts;
extern bool use_submotif;
extern bool use_em;
extern bool reestimate_mixing_parameters;
extern bool reestimate_error_rate;
extern bool use_logodds_score;
extern bool use_permutation_test;
extern bool use_bernoulli_read_method;
extern bool print_alignment;

const int min_matrix_len=3;
const int max_matrix_len=20;

#endif // PARAMETERS_HPP
