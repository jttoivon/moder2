/*

    MODER is a program to learn DNA binding motifs from SELEX datasets.
    Copyright (C) 2016  Jarkko Toivonen

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
