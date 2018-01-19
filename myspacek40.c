#include <time.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <inttypes.h>
#define MAX(a,b) (((a)>(b))?(a):(b))
#define MIN(a,b) (((a)>(b))?(b):(a))

char *VERSION = "spacek40 v0.173 FEB 22 2015";
/* multinomial 2 bug fixed */
/* from 0.144 lower memory use for kmer counting */

char *COMMAND;

/* GLOBAL VARIABLES */
float pseudocount = 0.0001;
short int Nlength = 20;
short int max_Nlength = 40;
short int max_width_of_pwm = 84;	/* Nlength * 2 + 2 */
long int max_number_of_sequences = 10000000;
short int number_of_files = 2;
double pvalue_cache[1001][1001];
short int head_to_tail = 0;
short int head_to_head = 1;
short int tail_to_tail = 2;
short int print_local_max = 0;
short int align_matches = 0;
short int pwm_align = 0;
short int contacts = 0;
short int rna = 0;
double local_max_min_percent = 0.1;
short int nocall = 0;
short int methylCGcompare = 0;
char *seed_story;
int core_length1=0;                   // Length of the boxed area in nucleotides
int core_length2=0;                   // Length of the boxed area in nucleotides
int core_distance=0;                  // Distance between two boxes, can be negative
int core_line_thickness = 4;
int core_line_offset = 4;            // Offset in y-direction in the top and bottom of the box

__uint128_t mask_ULL[42][42];	/* PRIMARY mask_ULL FOR EACH SEPARATE NUCLEOTIDE, VALUES GIVEN IN MAIN PROGRAM */
__uint128_t lowmask_ULL[42];	/* LOW MASK FOR EACH KMER, VALUES GIVEN IN MAIN PROGRAM */
__uint128_t highmask_ULL[42];	/* HIGH MASK FOR EACH KMER, VALUES GIVEN IN MAIN PROGRAM */

char *tf_kmers[] =
    { "AATCAA", "ACATGT", "ACCGCA", "ACCGGA", "AGATAA", "AGGTCA", "ATCGAT", "CAACAG", "CACCTG", "CACGCA", "CACGTC", "CACGTG", "CACTTA", "CAGCTG",
"CATAAA", "CATATG", "CATTCC", "CCATAT", "CCATTA", "CCCGCC", "CCCGGA", "CCGGAT", "CCGTTA", "CGAAAC", "CGTAAA", "CTAGTG", "CTGTCA", "GAAACC", "GAACAA", "GACCAC",
"GACGTC", "GAGGAA", "GCCACG", "GGCAAC", "GGCGCC", "GGGGAA", "GGTACA", "GGTGTG", "GTCACG", "TAAACA", "TAATTA", "TACGTA", "TATGCA", "TGACAG", "TGCATA", "TGCCAA",
"TGCGCA", "TGCGGG", "TGCTGA", "TGGAAA", "TTCTAG", "ATGCCC", "GTCCGC", "GTGAAA", "CCGCCA", "TCGCGA", "CGACCA", "CGCTGT", "ACCCAC", "ACCGGT", "CCATGG", "ATCAAA",
"AACGAT", "TTCGAA", "AAATAG", "TCTAGA", "TGCCCT", "CACGCC", "GATGCA", "TGACTC", "TGAGTC", "TGACAC", "TCCCCA", "TAAACG", "TAATTG", "CAATAA", "ACATGA", "CAAGGT",
"GTCCAA", "AAGTCA", "AGTTCA", "END" };
char *tf_names[] =
    { "PBX", "P53", "RUNX", "ETSI", "GATA", "NucRes1", "CUT", "SCRT", "TCF4", "EGR", "CREB3", "EboxHLH", "NKX", "NHLH", "HOX9to12", "bHLHcat", "TEAD",
"SRF", "PITX", "HighE2F", "ETSII", "SPDEF", "MYB", "IRF", "HOX13", "CTCF", "MEIS", "IRF2", "SOX", "GLI", "cre_bZIP", "ETSIII", "SP_KLF", "RFX", "LowE2F",
"NFKB1", "AR", "Tbox", "PAX", "FOX", "Homeo", "TEF_GMEB", "POU", "MEIS", "CUT2", "NFI", "CEBP", "GCM", "MAF", "NFAT", "HSF", "HIC", "HINFP", "PRD_IRF", "YY",
"ZBED", "ZBTB7", "ZIC", "ZNF143", "GRHL", "EBF", "TCF7_LEF", "HOMEZ", "HSFY", "MEF", "SMAD3", "TFAP", "SREBF", "CRE_CEBP", "prebZIP1", "prebZIP2", "oddFox_a",
"FOXO", "BARHLa", "BARHLb", "HOX9to12b", "IRX", "ESRR", "HNF4", "NRE", "VDR", "END" };
char *pwm_row_ids[] =
    { "Top_Strand_Base_Count_A", "Top_Strand_Base_Count_C", "Top_Strand_Base_Count_G", "Top_Strand_Base_Count_T", "Top_Strand_Maj_Groove_Phosphate",
"Top_Strand_Min_Groove_Phosphate", "Top_Strand_Maj_Groove_Sugar", "Top_Strand_Min_Groove_Sugar", "Top_Strand_Maj_Groove_Base", "Top_Strand_Min_Groove_Base",
"Bottom_Strand_Min_Groove_Base", "Bottom_Strand_Maj_Groove_Base", "Bottom_Strand_Min_Groove_Sugar", "Bottom_Strand_Maj_Groove_Sugar",
"Bottom_Strand_Min_Groove_Phosphate", "Bottom_Strand_Maj_Groove_Phosphate" };
long int *tf_kmer_values;
short int number_of_tf_kmers;

char *yesnocall[] = { "No", "Yes", "ND" };

char *dna_iupac = "ACGTRYMKWSBDHVN";
char *dna_rc_iupac = "TGCAYRKMWSVHDBN";
char *dna_bitiupac = "?ACMGRSVTWYHKDBN";
char *rna_iupac = "ACGURYMKWSBDHVN";
char *rna_rc_iupac = "UGCAYRKMWSVHDBN";
char *rna_bitiupac = "?ACMGRSVUWYHKDBN";

char *dnaforward = "ACGTN";
char *rnaforward = "ACGUN";
char *dnareverse = "TGCAN";
char *rnareverse = "UGCAN";
char *forward_lc = "acgtn";
short int iupac_length = 15;
long int **flank_kmer_count;

/* GLOBAL FLAGS */
short int count_both_instances_of_palindromic_hit = 0;
short int count_unequal_hits_only = 0;
short int count_only_forward_instance_of_palindromic_hit = 0;
short int count_only_reverse_instance_of_palindromic_hit = 0;
short int prefer_forward_strand = 0;
short int prefer_reverse_strand = 0;
short int last_count_was_forward = 0;
short int noname = 0;
short int paths = 0;
short int barcodelogo = 0;
short int gray_bars = 0;
short int barcodelogolabels = 0;
double scale_bars = 0;
short int max_scale_bar = 0;
short int circles = 1;
signed short int flank_kmer_pos = -100;

/* GLOBAL STRUCTURES */
struct dinucleotide_properties {
  short int number_of_dinucleotide_properties;
  char **dinucleotide_property_string;
  double **dinucleotide_property;
};
struct dinucleotide_properties di;

/* STRUCTURE DECLARATIONS AND INITIALIZATION SUBROUTINES */
struct flags {
  short int extendedoutput;
  short int remove_non_unique;
  short int print_counts;
  short int print_p_values;
  short int print_input_sequences;
  short int output_all_gap_lengths_in_one_line;
  short int print_nucleotides;
  short int print_frequencies;
  short int count_also_spaced_kmers;
  short int only_palindromes;
  short int dinucleotide_properties;
  short int kmer_table;
  short int information_content_output;
  short int even_background;
  short int complex_background;
  short int flank_with_pwm;
  short int flank;
  short int kmer_count;
};

/* RGB COLOR */
struct rgb_color {
  short int red;
  short int green;
  short int blue;
};

/* ALIGNMENT SCORE */
struct alignscore {
  long int **score;
  long int **count;
  long int **direct_repeat;
  long int **inverted_repeat;
};

/* ORIENTED MATCH */
struct oriented_match {
  short int position;
  signed short int strand;
  double score;
  short int id;
};

/* BIT REPRESENTATION OF IUPAC */
struct bitiupac_structure {
  short int length;
  __uint128_t *base;
  __uint128_t sequence_value_ULL;
};

/* STATISTICS FOR KMER-MONO CORRELATIONS */
struct similarity_stats {
  double correlation;
  long int sum_total_deviation;
  long int sum_relative_deviation;
};

struct sumtable {
  __uint128_t *sums;
  long int *counts;
  long int *max_min_counts;
  long int *max_max_counts;
  __uint128_t *max_max_kmer;
  __uint128_t *max_min_kmer;
  long int ***cloud_counts;
};

/* ADJACENT DINUC MARKOV MODEL */
struct adjacent_dinucleotide_model {
  char *name;
  short int width;
  double **fraction;
  double **mononuc_fraction;
};

/* BASE DEPENDENCY MATRIX */
struct base_dependency_matrix {
  char *name;
  short int width;
  long int ***incidence;
  double **total_relative_deviation;
  double **count_statistic_expected_total_relative_deviation;
  double **information_content;
  double **eo_correlation;
  double **permutated_correlation;
  double max_expected_relative_deviation;
  double max_relative_deviation;
};



/* DINUCLEOTIDE PROPERTIES */

/* COUNT PAIR */
struct countpair {
  long int sequence_value;
  short int contains_CpG;
  short int add_logo;
  short int gap_position;
  short int gap_width;
  long int x_count;
  long int y_count;
  short int x_local_max;
  short int y_local_max;
  long int sum;
};

/* SEQUENCE INCIDENCE TABLE */
struct sequence_incidence_table {
  __uint128_t sequence_value_ULL;
  long int incidence;
};

/* KMER INCIDENCE TABLE */
struct kmer_incidence_table {
  short int kmer_length;
  long int kmer;
  long int incidence;
  float max_enrichment;
  float local_max_enrichment;
  long int local_max_max_incidence;
  short int max_gap_length;
  short int local_max_max_gap_length;
  short int local_max;
  short int any_local_max;
  short int preferred;
};

/* KMER LIST */
struct kmer_list {
  char **kmer;
  __uint128_t *kmer_value;
  short int kmer_length;
  short int number_of_kmers;
};

/* KMER MATCHES */
struct kmer_matches {
  char **kmer;
  __uint128_t *kmer_value;
  short int kmer_length;
  short int number_of_kmers;
  long int number_of_sequences;
  short int **number_of_matches;
};

/* COUNT PWM */
struct count_pwm {
  char *name;
  short int width;
  long int max_counts;
  double **incidence;
};


short int count_pwm_init(struct count_pwm *i, char *name, short int width, double initial_value)
{
  short int maximum_width = max_width_of_pwm;
  short int counter;
  short int counter2;
  (*i).name = malloc(1000);
  strcpy((*i).name, name);
  (*i).width = width;
  (*i).max_counts = initial_value;
  (*i).incidence = malloc(sizeof(double *) * (5 + contacts * 12) + 5);
  for (counter = 0; counter < 5 + contacts * 12; counter++) {
    (*i).incidence[counter] = malloc(sizeof(double) * maximum_width + 5);
    for (counter2 = 0; counter2 < maximum_width; counter2++)
      (*i).incidence[counter][counter2] = initial_value;
  }
  return (0);
}


struct pairwise_correlation {
  short int first_base;
  short int second_base;
  double delta_ic;
  short int max_dinucleotide;
  short int min_dinucleotide;
  double min_fold_change;
  double max_fold_change;
};

/* NORMALIZED PWM */
struct normalized_pwm {
  char *name;
  char *seed;
  short int width;
  long int max_counts;
  double *information_content;
  short int *original_position;
  double *position_score;
  long int *total_counts_for_column;
  double **fraction;
  struct pairwise_correlation *pairwise_correlation;
  short int negative_values_allowed;
  struct oriented_match match;
};
short int normalized_pwm_init(struct normalized_pwm *i, char *name, short int width, double initial_value)
{
  short int maximum_width = max_width_of_pwm;
  short int counter;
  short int counter2;
  (*i).negative_values_allowed = 0;
  (*i).pairwise_correlation = malloc(sizeof(struct pairwise_correlation) * 10 + 5);
  for (counter = 0; counter < 10; counter++) {
    (*i).pairwise_correlation[counter].first_base = 0;
    (*i).pairwise_correlation[counter].delta_ic = 0;
    (*i).pairwise_correlation[counter].second_base = 0;
    (*i).pairwise_correlation[counter].max_dinucleotide = 0;
    (*i).pairwise_correlation[counter].min_dinucleotide = 0;
    (*i).pairwise_correlation[counter].max_fold_change = 0;
    (*i).pairwise_correlation[counter].min_fold_change = 0;
  }
  (*i).name = malloc(100);
  strcpy((*i).name, name);
  (*i).seed = malloc(1000);
  strcpy((*i).seed, "UNKNOWN");
  (*i).width = width;
  (*i).max_counts = initial_value;
  (*i).fraction = malloc(sizeof(double *) * (5 + contacts * 12) + 5);
  (*i).information_content = malloc(sizeof(double) * maximum_width + 5);
  (*i).position_score = malloc(sizeof(double) * maximum_width + 5);
  (*i).original_position = malloc(sizeof(short int) * maximum_width + 5);
  (*i).total_counts_for_column = malloc(sizeof(long int) * maximum_width + 5);

  for (counter = 0; counter < 5 + contacts * 12; counter++) {
    (*i).fraction[counter] = malloc(sizeof(double) * maximum_width + 5);
    for (counter2 = 0; counter2 < maximum_width; counter2++)
      (*i).fraction[counter][counter2] = initial_value;
  }
  for (counter2 = 0; counter2 < maximum_width; counter2++) {
    (*i).information_content[counter2] = 0;
    (*i).position_score[counter2] = 0;
    (*i).original_position[counter2] = counter2;
    (*i).total_counts_for_column[counter2] = 0;
  }
  (*i).match.position = 0;
  (*i).match.strand = 0;
  (*i).match.score = 0;
  return (0);
}


/* DINUCLEOTIDE PROPERTIES MATRIX */
struct dinucleotide_properties_matrix {
  char *name;
  short int width;
  short int number_of_dinucleotide_properties;
  long int max_counts;
  long int **count;
  double **score;
  char **dinucleotide_property_string;
  short int query_sequence_length;
};



/* COUNT CONNECTING MATRIX */
struct count_connecting_matrix {
  char *name;
  short int width;
  short int height;
  long int number_of_total_matches;
  long int one_hit_matches;
  long int two_hit_matches;
  long int **incidence;
};


short int count_connecting_matrix_init(struct count_connecting_matrix *i, char *name, short int width, double initial_value)
{
  short int maximum_width = max_width_of_pwm;
  short int counter;
  short int counter2;
  (*i).number_of_total_matches = 0;
  (*i).one_hit_matches = 0;
  (*i).two_hit_matches = 0;
  (*i).name = malloc(100);
  strcpy((*i).name, name);
  (*i).width = width;
  (*i).incidence = malloc(sizeof(double *) * 5 + 5);
  for (counter = 0; counter < 5; counter++) {
    (*i).incidence[counter] = malloc(sizeof(double) * maximum_width + 5);
    for (counter2 = 0; counter2 < maximum_width; counter2++)
      (*i).incidence[counter][counter2] = initial_value;
  }
  return (0);
}

/* NORMALIZED CONNECTING MATRIX */
struct normalized_connecting_matrix {
  char *name;
  short int width;
  short int height;
  long int *orientation_count;
  double *orientation_fraction;
  double **fraction;
};
short int normalized_connecting_matrix_init(struct normalized_connecting_matrix *i, char *name, short int width, double initial_value)
{
  short int maximum_width = max_width_of_pwm;
  short int counter;
  short int counter2;
  (*i).name = malloc(100);
  strcpy((*i).name, name);
  (*i).width = width;
  (*i).fraction = malloc(sizeof(double *) * 5 + 5);
  (*i).orientation_fraction = malloc(sizeof(double) * 5 + 5);
  (*i).orientation_count = malloc(sizeof(long int) * 5 + 5);
  for (counter = 0; counter < 4; counter++) {
    (*i).fraction[counter] = malloc(sizeof(double) * maximum_width + 5);
    for (counter2 = 0; counter2 < maximum_width; counter2++)
      (*i).fraction[counter][counter2] = initial_value;
  }
  for (counter2 = 0; counter2 < 4; counter2++) {
    (*i).orientation_fraction[counter2] = 0;
    (*i).orientation_count[counter2] = 0;
  }
  return (0);
}

struct hit_position_matrix {
  char *name;
  short int width;
  long int **incidence;
  double **fraction;
};

short int hit_position_matrix_init(struct hit_position_matrix *i, char *name, short int width, double initial_value)
{
  short int maximum_width = max_Nlength;
  short int counter;
  short int counter2;
  (*i).name = malloc(100);
  strcpy((*i).name, name);
  (*i).width = width;
  (*i).fraction = malloc(sizeof(double *) * 3 + 5);
  (*i).incidence = malloc(sizeof(long int *) * 3 + 5);
  for (counter = 0; counter < 2; counter++) {
    (*i).fraction[counter] = malloc(sizeof(double) * maximum_width + 5);
    (*i).incidence[counter] = malloc(sizeof(long int) * maximum_width + 5);
    for (counter2 = 0; counter2 < maximum_width; counter2++) {
      (*i).fraction[counter][counter2] = initial_value;
      (*i).incidence[counter][counter2] = initial_value;
    }
  }
  return (0);
}

struct match {
  short int *position;
  double *score;
};
short int match_init(struct match *i, short int width)
{
  short int maximum_width = Nlength + 10;
  short int counter;
  short int counter2;
  (*i).position = malloc(sizeof(short int) * maximum_width + 5);
  (*i).score = malloc(sizeof(double) * maximum_width + 5);
  for (counter = 0; counter < maximum_width; counter++) {
    (*i).position[counter] = 0;
    (*i).score[counter] = 0;
  }
  return (0);
}

/* SVG_TILE */
struct svg_tile {
  short int x;
  short int y;
  short int height;
  short int width;
  short int red;
  short int green;
  short int blue;
  short int stroke;
};

/* SUBROUTINES */

/* SUBROUTINE THAT ESCAPES OFFENDING CHARACTERS FROM STRINGS FOR INCLUSION IN SVG COMMENTS */
char *svgsafe(char *string)
{
  signed long int counter;
  signed long int counter2;

  for (counter = strlen(string); counter >= 0; counter--) {
    if (string[counter] == '-' && string[counter + 1] == '-')
      string[counter] = '_';
  }
  return (string);
}


/* SUBROUTINE THAT ADDS NUCLEOTIDE PATHS */
short int Add_nucleotide_paths(FILE * outfile)
{
  fprintf(outfile,
	  "<g display=\"none\"> <path style=\"stroke:none\" fill=\"green\" id=\"A\" d=\"m 6.8910692,-2.720526 -3.7670144,0 -0.3782159,1.030503 0.3782159,0 q 0.5597597,0 0.7942476,0.230832 0.2344959,0.230832 0.2420559,0.618304 0,0.370975 -0.2344879,0.601807 Q 3.6913745,-0.008248 3.1240548,0 L 1.0363116,0 Q 0.47655181,0 0.24205591,-0.230832 0.00756801,-0.461664 0,-0.849135 0,-1.228359 0.24962391,-1.467439 0.49923981,-1.706511 1.0741276,-1.690023 l 2.435703,-6.611698 -1.0136156,0 q -0.5597598,0 -0.7942477,-0.230832 -0.2344959,-0.23084 -0.2420559,-0.618304 0,-0.379223 0.2344879,-0.610055 0.2344959,-0.230832 0.8018157,-0.23908 l 3.3509827,0.0082 3.0862308,8.301721 q 0.544632,0 0.718608,0.131904 0.34796,0.272056 0.34796,0.717232 0,0.370975 -0.234496,0.601807 -0.234488,0.230832 -0.794248,0.23908 l -2.0877433,0 q -0.5597597,0 -0.7942556,-0.230832 -0.2344879,-0.230832 -0.2420559,-0.618303 0,-0.370976 0.2344959,-0.601816 0.2344879,-0.230832 0.8018156,-0.239072 l 0.3782159,0 -0.3706559,-1.030503 z M 6.2481095,-4.410549 4.999998,-7.757618 l -1.2556715,3.347069 2.503783,0 z\" />\n");
  fprintf(outfile,
	  "<path style=\"stroke:none\" id=\"C\" fill=\"blue\" d=\"m 7.9562048,-9.347881 q 0.1733519,-0.204239 0.3740799,-0.306359 0.2007359,-0.10212 0.4470878,-0.10212 0.410584,0 0.666056,0.243512 0.255472,0.243519 0.2646,0.824815 l 0,1.366839 q 0,0.581295 -0.25548,0.824815 -0.255472,0.243512 -0.675176,0.251368 -0.3740879,0 -0.6021918,-0.180672 Q 7.9470768,-6.606355 7.8375889,-7.10125 7.7737489,-7.43117 7.582117,-7.611842 7.2080291,-7.965338 6.5419654,-8.177434 5.8759096,-8.389529 5.1916059,-8.389529 q -0.8394156,0 -1.5419754,0.314215 -0.7025517,0.314216 -1.2408715,1.0212 -0.5383197,0.706983 -0.5383197,1.681046 l 0,1.044767 q 0,1.162599 0.9762716,1.940279 0.9762796,0.777679 2.7281029,0.777679 1.0401436,0 1.7609513,-0.243512 0.4197038,-0.141399 0.8941596,-0.557735 0.2919679,-0.251368 0.4561998,-0.322072 0.16424,-0.07072 0.374088,-0.07856 0.374088,0 0.656936,0.24352 0.282848,0.243512 0.282848,0.57344 0,0.329927 -0.383208,0.706983 -0.556576,0.54988 -1.4324878,0.864095 Q 7.0072932,0 5.5839418,0 3.9233584,0 2.591239,-0.589152 1.5145994,-1.060471 0.75729571,-2.073814 0,-3.087158 0,-4.296885 L 0,-5.38878 Q 0,-6.645643 0.67518371,-7.729682 1.3503675,-8.813729 2.563871,-9.410736 q 1.2135035,-0.597008 2.554743,-0.589152 0.8120396,0 1.5145913,0.157112 0.7025598,0.157104 1.3229995,0.494887 z\" />\n");
  fprintf(outfile,
	  "<path style=\"stroke:none\" id=\"G\" fill=\"orange\" d=\"m 9.2454565,-3.181462 0,2.160255 Q 8.0399288,-0.432048 7.1986567,-0.219952 6.3573851,-0.007856 5.3686815,0 4.007034,0 2.8448665,-0.369208 1.9342188,-0.652 1.413843,-1.084055 0.89346721,-1.516103 0.44248341,-2.356638 -0.00850839,-3.197173 1.6360563e-4,-4.249797 l 0,-1.162607 q 0,-1.610374 1.08411159437,-2.875101 1.4570394,-1.712487 3.9635024,-1.712487 0.7285277,0 1.3789914,0.133544 0.6504638,0.133544 1.2488955,0.400624 0.3555843,-0.290648 0.6938322,-0.290648 0.3902798,0 0.6331198,0.24352 0.24284,0.243519 0.24284,0.824823 l 0,1.044775 q 0,0.581304 -0.24284,0.824823 -0.24284,0.24352 -0.6417918,0.251368 -0.3122239,0 -0.5377202,-0.172816 Q 7.6496485,-6.857811 7.5368966,-7.258442 7.4241526,-7.659066 7.2680407,-7.808322 7.0338728,-8.03613 6.444113,-8.208945 5.8543613,-8.381769 5.1084896,-8.389625 q -1.0580876,0 -1.8386473,0.392775 -0.5550638,0.290648 -1.0233996,1.013352 -0.4683358,0.722703 -0.4683358,1.571094 l 0,1.162607 q 0,1.296143 0.8239197,1.963863 0.8239276,0.667711 2.6365589,0.675567 1.2228715,0 2.2289271,-0.416343 l 0,-1.154752 -1.8126313,0 q -0.6417917,0 -0.9106476,-0.219951 -0.2688639,-0.219952 -0.2775359,-0.58916 0,-0.361352 0.2688559,-0.581303 0.2688639,-0.219952 0.9193276,-0.227808 l 3.1569272,0.008 q 0.6418,0 0.910656,0.219952 0.268856,0.219951 0.2775365,0.581303 0,0.2828 -0.1821365,0.494896 -0.182128,0.212096 -0.572408,0.314215 z\" />\n");
  if (rna == 0)
    fprintf(outfile,
	    "<path style=\"stroke:none\" id=\"T\" fill=\"red\" d=\"m 5.9510296,-8.301721 0,6.611698 1.2994395,0 q 0.6967997,0 0.9886956,0.230832 0.2919039,0.230832 0.3013199,0.618304 0,0.370975 -0.2919039,0.601807 Q 7.9566848,-0.008248 7.2504691,0 L 2.7306949,0 Q 2.0338952,0 1.7419993,-0.230832 1.4500954,-0.461664 1.4406794,-0.849135 q 0,-0.370976 0.2919039,-0.601816 0.2918959,-0.230832 0.9981116,-0.239072 l 1.2900155,0 0,-6.611698 -2.0903912,0 0,1.599334 q 0,0.610056 -0.2636559,0.865624 -0.2636479,0.255567 -0.70621569,0.263807 -0.4237278,0 -0.6873757,-0.255567 Q 0.00941601,-6.084083 0,-6.702387 l 0,-3.297605 9.9999965,0.0082 0,3.289357 q 0,0.610056 -0.263656,0.865624 -0.263648,0.255567 -0.706216,0.263807 -0.4237279,0 -0.6873758,-0.255567 -0.2636559,-0.25556 -0.2730719,-0.873864 l 0,-1.599334 -2.1186472,0 z\" />\n</g>\n");
  else
    fprintf(outfile,
	    "<path style=\"stroke:none\" id=\"U\" fill=\"red\" d=\"m 1.0966549,-8.8225952 0,5.1428531 c 0,2.3025166 1.5799292,3.68066910008 3.9219361,3.68066910008 2.6394046,0 3.9405131,-1.31092600008 3.9405131,-3.98318990008 l 0,-4.8403323 c 0.037178,0 0.092936,0 0.130111,0 0.6691451,0 0.9107811,-0.1344542 0.9107811,-0.5882348 0,-0.386554 -0.204461,-0.588235 -0.650557,-0.588235 -0.1487001,0 -0.278812,0 -0.42751,0 l -2.0074259,0 c -0.6691449,0 -0.985131,0.06722 -0.985131,0.588235 0,0.4873935 0.353161,0.5882348 1.1895921,0.5882348 0.1486989,0 0.3159855,0 0.4832717,0 l 0,4.6554585 c 0,1.9999961 -0.6691472,2.9915941 -2.4907082,2.9915941 -0.9107805,0 -1.710039,-0.3865551 -2.1375483,-1.0588229 C 2.5650578,-2.8898264 2.5650574,-3.5452893 2.5650574,-4.318397 l 0,-4.5041982 c 0.1486989,0 0.2788108,0 0.3903349,0 0.7992565,0 1.1338299,-0.1008413 1.1338299,-0.5882348 0,-0.571428 -0.3717481,-0.588235 -1.1338299,-0.588235 l -1.8587374,0 C 0.37174808,-9.999065 0,-9.965445 0,-9.41083 c 0,0.4537806 0.2602238,0.5882348 0.87360628,0.5882348 0.07435,0 0.14869912,0 0.22304862,0 z\" />\n</g>\n");

}

struct plot_point {
  double x;
  double y;
};


/* SUBROUTINE THAT GENERATES FORWARD SEQUENCE VALUE as long int */
long int Generate_sequence_value(char *searchstring)
{
  short int query_sequence_length = strlen(searchstring);
  long int query_sequence_value;
  short int position;
  short int nucleotide_value;
  long int position_value;

  for (query_sequence_value = 0, position = 0, position_value = pow(4, query_sequence_length - 1); position < query_sequence_length;
       position++, position_value /= 4) {
    for (nucleotide_value = 0; nucleotide_value < 4 && searchstring[position] != dnaforward[nucleotide_value]; nucleotide_value++) ;
    if (nucleotide_value == 4) {
      printf("\nERROR IN QUERY SEQUENCE\n");
      exit(1);
    }
    query_sequence_value += position_value * nucleotide_value;
  }
  return (query_sequence_value);
}


/* SUBROUTINE THAT CHECKS IF CpG DINUCLEOTIDE FREQUENCY IS DIFFERENT */

void
draw_box(FILE* outfile, int offset, int start, int end)
{
  char style[1000];
  double y1 = 0 - core_line_offset;
  double y2 = 100 + core_line_offset;
  snprintf(style, 999, "stroke:rgb(0,0,0); stroke-linecap:square; stroke-width:%i", core_line_thickness);
  /* Top line */
  fprintf(outfile, "<line x1=\"%i\" y1=\"%.2f\" x2=\"%i\" y2=\"%.2f\" style=\"%s\" />\n",
	  offset + 20*start, y1, offset + 20*end, y1, style);
  /* Bottom line */
  fprintf(outfile, "<line x1=\"%i\" y1=\"%.2f\" x2=\"%i\" y2=\"%.2f\" style=\"%s\" />\n",
	  offset + 20*start, y2, offset + 20*end, y2, style);
  /* Left line */
  fprintf(outfile, "<line x1=\"%i\" y1=\"%.2f\" x2=\"%i\" y2=\"%.2f\" style=\"%s\" />\n",
	  offset + 20*start, y1, offset + 20*start, y2, style);
  /* Right line */
  fprintf(outfile, "<line x1=\"%i\" y1=\"%.2f\" x2=\"%i\" y2=\"%.2f\" style=\"%s\" />\n",
	  offset + 20*end, y1, offset + 20*end, y2, style);
}

// draw vertical line in dashed stroke
// draw vertical line in dashed stroke
void
draw_vertical_line(FILE* outfile, int offset, int pos)
{
  char style[1000];
  double y1 = 0 - core_line_offset;
  double y2 = 100 + core_line_offset;
  snprintf(style, 999, "stroke:rgb(0,0,0); stroke-linecap:square; stroke-width:%i; stroke-dasharray: 10;", core_line_thickness);
  fprintf(outfile, "<line x1=\"%i\" y1=\"%.2f\" x2=\"%i\" y2=\"%.2f\" style=\"%s\" />\n",
	  offset + 20*pos, y1, offset + 20*pos, y2, style);
}

/* SUBROUTINE THAT GENERATES AN SVG LOGO FILE */
short int Svg_logo(char *filename, short int number_of_pwms, struct normalized_pwm **n, struct count_connecting_matrix **cm,
		   struct dinucleotide_properties_matrix *d, struct base_dependency_matrix *bd_background, struct base_dependency_matrix *bd_signal,
		   struct base_dependency_matrix *bd_expected, struct alignscore *all_hits_align_scores, double lambda, short int warning)
{
  FILE *outfile;
  outfile = fopen(filename, "w");
  if (warning > 2 || warning < 0)
    warning = 2;

  short int circles = 0;
  short int top_backbone_contacts = 0;
  short int top_base_contacts = 0;
  short int bottom_backbone_contacts = 0;
  short int bottom_base_contacts = 0;
  short int current_contact_position;
  short int contacts_from_pwm;
  short int contacts_defined;
  short int color_start;
  signed short int dot_offset;
  short int CG_difference;
  short int current_pwm = 0;
  short int top_position;
  short int counter;
  short int nucleotide_value;
  short int height_scale;
  signed short int start_rectangle;
  double rectangle_cutoff = 0.25;
  short int pwm_position = 0;
  short int offset = 0;

  short int **colorcode;
  colorcode = malloc(sizeof(short int *) * 2 + 5);
  colorcode[0] = malloc(sizeof(short int) * 4 + 5);
  colorcode[1] = malloc(sizeof(short int) * 4 + 5);
  for (counter = 0; counter < 4; counter++) {
    colorcode[0][counter] = 0;
    colorcode[1][counter] = 0;
  }

  double **order;
  order = malloc(sizeof(double *) * 4 + 5);
  for (counter = 0; counter < 3; counter++)
    order[counter] = malloc(sizeof(double) * 6 + 5);
  double swap;
  double font_position;
  char *nucleotide_char;
  char *nucleotide_iupac;
  char *nucleotide_bitiupac;

  if (rna == 1) {
    nucleotide_char = rnaforward;
    nucleotide_iupac = dna_iupac;
    nucleotide_bitiupac = dna_bitiupac;
  } else {
    nucleotide_char = dnaforward;
    nucleotide_iupac = rna_iupac;
    nucleotide_bitiupac = rna_bitiupac;
  }

  short int iupac_bits = 0;
  short int rgbcolors[5][4];
  rgbcolors[0][0] = 0;
  rgbcolors[0][1] = 128;
  rgbcolors[0][2] = 0;
  rgbcolors[1][0] = 0;
  rgbcolors[1][1] = 0;
  rgbcolors[1][2] = 255;
  rgbcolors[2][0] = 255;
  rgbcolors[2][1] = 165;
  rgbcolors[2][2] = 0;
  rgbcolors[3][0] = 255;
  rgbcolors[3][1] = 0;
  rgbcolors[3][2] = 0;
  rgbcolors[4][0] = 211;
  rgbcolors[4][1] = 211;
  rgbcolors[4][2] = 211;

  char **colors;
  colors = malloc(sizeof(char *) * 16 + 5);
  for (counter = 0; counter < 15; counter++)
    colors[counter] = malloc(sizeof(char) * 20 + 5);
  strcpy(colors[0], "green");
  strcpy(colors[1], "blue");
  strcpy(colors[2], "orange");
  strcpy(colors[3], "red");
  strcpy(colors[4], "black");
  strcpy(colors[5], "lightgray");
  strcpy(colors[6], "white");
  strcpy(colors[7], "cornflowerblue");
  strcpy(colors[8], "none");
  strcpy(colors[9], "orchid");
  strcpy(colors[10], "midnightblue");
  strcpy(colors[11], "aliceblue");
  strcpy(colors[12], "saddlebrown");
  strcpy(colors[13], "moccasin");

/* INITIALIZES TF KMER SEQUENCE VALUES */
  for (number_of_tf_kmers = 0; strcmp(tf_kmers[number_of_tf_kmers], "END") != 0; number_of_tf_kmers++) ;
  tf_kmer_values = malloc(sizeof(char *) * number_of_tf_kmers + 5);
  for (counter = 0; counter < number_of_tf_kmers; counter++) {
    tf_kmer_values[counter] = Generate_sequence_value(tf_kmers[counter]);
/*  printf("\n%s %li", tf_kmers[counter], tf_kmer_values[counter]); */
  }

  char *font = "Courier";
  fprintf(outfile,
	  "<?xml version=\"1.0\" standalone=\"no\"?>\n<!DOCTYPE svg PUBLIC \"-//W3C//DTD SVG 1.1//EN\"\n\"http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd\">\n");
  fprintf(outfile, "<!--%s : command %s -->\n", svgsafe(VERSION), svgsafe(COMMAND));

  if (cm == '\0') {
    if (noname == 1 && number_of_pwms == 1) {
      int width=(*n)[0].width * 20;
      fprintf(outfile, "<svg width=\"");
      fprintf(outfile, "%i", width);
      if (core_length1 > 0) {
	double extra = core_line_thickness/2.0 + core_line_offset;
	//double extra = 0;
	int height=100 + 2*extra + 130 * contacts;
	fprintf(outfile,
		"\" height=\"%i\" viewBox=\"0 %.2f %i %i\" version=\"1.1\" xmlns=\"http://www.w3.org/2000/svg\" xmlns:xlink=\"http://www.w3.org/1999/xlink\">\n",
		height, -extra + 300 * current_pwm, width, height);
      }
      else
	fprintf(outfile,
		"\" height=\"%i\" x=\"0\" y=\"%i\" version=\"1.1\" xmlns=\"http://www.w3.org/2000/svg\" xmlns:xlink=\"http://www.w3.org/1999/xlink\">\n",
		100 + 130 * contacts, 300 * current_pwm);
	
    }

    else
      fprintf(outfile,
	      "<svg width=\"2500\" height=\"%i\" x=\"0\" y=\"%i\" version=\"1.1\" xmlns=\"http://www.w3.org/2000/svg\" xmlns:xlink=\"http://www.w3.org/1999/xlink\">\n",
	      number_of_pwms * 100 + 130 * contacts, 300 * current_pwm);
  } else
    fprintf(outfile,
	    "<svg width=\"2500\" height=\"5000\" x=\"0\" y=\"%i\" version=\"1.1\" xmlns=\"http://www.w3.org/2000/svg\" xmlns:xlink=\"http://www.w3.org/1999/xlink\">\n",
	    300 * current_pwm);

  fprintf(outfile, "<title>%s %s</title>\n", svgsafe(VERSION), svgsafe(COMMAND));
  if (paths == 1)
    Add_nucleotide_paths(outfile);	/* Adds nucleotide paths */

/* GENERATES LOGO */
  for (top_position = 20; current_pwm < number_of_pwms; current_pwm++) {

/* GENERATES SHAPE HEATMAP, DINUCLEOTIDE LINE AND BASE DEPENDENCY HEATMAP */
/*     if (current_pwm == 3) */
/*       if (d != '\0') { */
/* 	offset = Svg_dinucleotide_heatmap(outfile, 0, current_pwm * 120 - 40, d, bd_signal, bd_expected) - 350; */
/* 	/\* Dinucleotide_line(outfile, n[current_pwm], 0, current_pwm * 120 - 120,1); *\/ */
/* 	offset = Svg_base_dependency_heatmap(outfile, (*n[current_pwm]).seed, 0, current_pwm * 120 + offset, bd_signal, bd_expected, n) - 500; */
/*       } */

/* /\* DRAWS PIECHART *\/ */
/*     if (current_pwm == 0) */
/*       if (cm != '\0') { */
/* 	double *hitcounts; */
/* 	hitcounts = malloc(sizeof(double) * 4 + 5); */
/* 	hitcounts[0] = (*cm[1]).one_hit_matches; */
/* 	hitcounts[1] = (*cm[1]).two_hit_matches; */
/* 	hitcounts[2] = (*cm[1]).number_of_total_matches - (*cm[1]).one_hit_matches - (*cm[1]).two_hit_matches; */
/* 	char *hitnames[] = { "One hit", "Two hits", "More hits" }; */
/* 	char *hitcolors[] = { "lightgreen", "lightsalmon", "red" }; */
/* 	Svg_piechart(outfile, "hit_piechart", 100, 60 + current_pwm + offset, 50, hitcounts, hitnames, hitcolors, 0, 2); */
/* 	free(hitcounts); */
/*       } */

/* /\* PRINTS SPACING HEATMAPS *\/ */
/*     if (current_pwm == 6) */
/*       if (cm != '\0') { */
/* 	Svg_heatmap(outfile, 180, current_pwm * 120 + offset, -1, -1, -1, cm[0]); */
/* 	Svg_heatmap(outfile, 180, current_pwm * 120 + 150 + offset, -1, -1, -1, cm[1]); */
/* 	Svg_heatmap(outfile, 180, current_pwm * 120 + 300 + offset, 0, 1, -1, cm[2]); */
/* 	offset += 470; */
/*       } */

    if (cm != '\0')
      fprintf(outfile, "<g id=\"group%i\" transform=\"translate(%i, %i)\" >", current_pwm, 0 + contacts * 20, current_pwm * 120 + offset);
    else
      fprintf(outfile, "<g id=\"group%i\" transform=\"translate(%i, %i)\" >", current_pwm, ((*(n[current_pwm])).match.position + contacts) * 20,
	      current_pwm * 95 + offset);

/* DRAWS KMER LINE */
    if (current_pwm > 0)
      ;//Kmerlines(outfile, (*(n[current_pwm])).name, (*n[current_pwm]).fraction, (*n[current_pwm]).width, 0, 0, 6, 2, 2.0);
    else
      printf("\n");

    for (pwm_position = 0; pwm_position < (*(n[current_pwm])).width; pwm_position++) {

/* DETERMINES ORDER BY BUBBLE SORT*/
      double positive_sum = 0;
      double negative_sum = 0;
      for (counter = 0; counter < 4; counter++) {
	order[0][counter] = counter;
	order[1][counter] = (*(n[current_pwm])).fraction[counter][pwm_position];
	if (order[1][counter] > 0)
	  positive_sum = positive_sum + order[1][counter];
	else
	  negative_sum = negative_sum + order[1][counter];
      }
      for (counter = 0; counter < 3; counter++) {
	for (nucleotide_value = counter; nucleotide_value < 4; nucleotide_value++) {
	  if (order[1][counter] == (0.0 / 0.0) || order[1][counter] == (1.0 / 0.0))
	    order[1][counter] = 0;
	  if (order[1][counter] < order[1][nucleotide_value]) {
	    swap = order[0][counter];
	    order[0][counter] = order[0][nucleotide_value];
	    order[0][nucleotide_value] = swap;
	    swap = order[1][counter];
	    order[1][counter] = order[1][nucleotide_value];
	    order[1][nucleotide_value] = swap;
	  }
	}
      }

/* CHECKS IF THERE IS DIFFERENCE IN CG OR GC FREQUENCY BETWEEN CONTROL AND SAMPLE */
      /* if (methylCGcompare == 1) */
      /* 	CG_difference = CG_difference_call(bd_background, bd_signal, pwm_position, colorcode); */

      for (font_position = 0, nucleotide_value = 0, iupac_bits = 0; nucleotide_value < 4; nucleotide_value++) {
	if (order[1][nucleotide_value] > 0 || (order[1][nucleotide_value] < 0 && (*(n[current_pwm])).negative_values_allowed == 1)) {
	  if (paths == 0 && barcodelogo == 0) {
// if (current_pwm == 3) CG_difference = 1;

/* PRINTS OUT SCALED FONT NUCLEOTIDES */
	    fprintf(outfile, " <text  x=\"%i\" y=\"%f", pwm_position * 20,
		    font_position / (order[1][nucleotide_value] * 4.5) + top_position - order[1][0] * 0.9);
/* else fprintf(outfile, "%f", font_position); */
	    fprintf(outfile, "\" fill = \"%s\" stroke=\"", colors[(int)order[0][nucleotide_value]]);
	    if (colorcode[0][(int)order[0][nucleotide_value]] == 0 || methylCGcompare == 0 || current_pwm > 5)
	      fprintf(outfile, "%s", colors[(int)order[0][nucleotide_value]]);
	    else {
	      if (colorcode[0][(int)order[0][nucleotide_value]] == 1)
		fprintf(outfile, "%s", "black\" stroke-width = \"1.5");
	      if (colorcode[0][(int)order[0][nucleotide_value]] == -1)
		fprintf(outfile, "%s", "black\" opacity = \"0.25");
	    }
	    fprintf(outfile, "\" font-size=\"30\" font-family = \"");
	    fprintf(outfile, "%s\" ", font);
	    if (contacts == 1 && pwm_align == 1 && current_pwm < 2 && (*(n[current_pwm])).position_score[pwm_position] < rectangle_cutoff)
	      fprintf(outfile, "opacity = \"0.1\" ");
	    fprintf(outfile, "transform = \"scale(1, ");
	    fprintf(outfile, "%f", order[1][nucleotide_value] * 4.5);
	    fprintf(outfile, ")\" >");
	    fprintf(outfile, "%c", nucleotide_char[(int)order[0][nucleotide_value]]);
	    fprintf(outfile, "</text>\n");
	    font_position += (order[1][nucleotide_value] * 90);
	  } else if (barcodelogo == 1) {
	    /* PRINTS OUT BARCODE LOGO RECTANGLE */
	    if (max_scale_bar == 1)
	      height_scale = 0;
	    else
	      height_scale = nucleotide_value;
	    fprintf(outfile, " <rect x=\"%.3f\" y=\"%.2f\" width=\"%.3f", pwm_position * 20 + font_position,
		    (double)top_position + 25 * (1 - scale_bars * order[1][height_scale]), order[1][nucleotide_value] * 15);
	    fprintf(outfile, "\" height=\"");
	    fprintf(outfile, "%.2f", 50 * (scale_bars == 0) + 50 * (scale_bars * order[1][height_scale]));
	    fprintf(outfile, "\" style=\"");
	    if (gray_bars == 0)
	      fprintf(outfile, "fill:%s;stroke-width:%i;stroke:%s\"/>\n", colors[(int)order[0][nucleotide_value]], 1,
		      colors[(int)order[0][nucleotide_value]]);
	    else {
	      fprintf(outfile, "fill:rgb(%.0f,%.0f,%.0f);stroke-width:%i;stroke:rgb(0,0,0)\"/>\n",
		      order[1][nucleotide_value] * rgbcolors[(int)order[0][nucleotide_value]][0] + (1 - order[1][nucleotide_value]) * rgbcolors[4][0],
		      order[1][nucleotide_value] * rgbcolors[(int)order[0][nucleotide_value]][1] + (1 - order[1][nucleotide_value]) * rgbcolors[4][1],
		      order[1][nucleotide_value] * rgbcolors[(int)order[0][nucleotide_value]][2] + (1 - order[1][nucleotide_value]) * rgbcolors[4][2],
		      0);
	    }
	    if (order[1][nucleotide_value] > order[1][0] / 2)
	      iupac_bits |= 1 << (int)order[0][nucleotide_value];
	    if ((nucleotide_value == 3 || (nucleotide_value < 4 && order[1][nucleotide_value + 1] == 0)) && barcodelogo == 1
		&& barcodelogolabels == 1)
	      fprintf(outfile, " <text  x=\"%i\" y=\"%i\" fill = \"white\" stroke =\"white\" font-size=\"15\" >%c</text>\n", pwm_position * 20 + 2,
		      top_position + 30, nucleotide_bitiupac[iupac_bits]);
	    font_position += (order[1][nucleotide_value] * 15);

	  } else if (paths == 1) {
	    /* PRINTS OUT SCALED PATH NUCLEOTIDES */
	    fprintf(outfile, "<use xlink:href=\"#%c\" ", nucleotide_char[(int)order[0][nucleotide_value]]);
	    if ((*(n[current_pwm])).negative_values_allowed == 0) {
	      fprintf(outfile, " transform=\"translate(%i,%f) scale(2,%f)\" visibility=\"visible\" />\n", pwm_position * 20 + offset,
		      font_position + (order[1][nucleotide_value] * 100), order[1][nucleotide_value] * 10);
	      font_position += (order[1][nucleotide_value] * 100);
	    } else {
	      /* printf("\n%f",order[1][nucleotide_value]); */
	      if (order[1][nucleotide_value] < 0)
		order[1][nucleotide_value] = -order[1][nucleotide_value];
	      fprintf(outfile, " transform=\"translate(%i,%f) scale(2,%f)\" visibility=\"visible\" />\n", pwm_position * 20 + offset,
		      (1 - positive_sum) * 50 + font_position + (order[1][nucleotide_value] * 50), order[1][nucleotide_value] * 5);
	      font_position += (order[1][nucleotide_value] * 50);
	    }

	  }

	}
      }

/* PRINTS BASE CONTACTS */
      if (contacts == 1 && ((pwm_align == 0) || (pwm_align == 1 && current_pwm == 2))) {
	current_contact_position = pwm_position;
	dot_offset = 0;
	if (pwm_align == 0)
	  contacts_from_pwm = current_pwm;
	else {
	  contacts_from_pwm = 0;
	  dot_offset--;
	}

	for (color_start = 7; contacts_from_pwm <= current_pwm - pwm_align; contacts_from_pwm++, color_start += 2, dot_offset += 2) {
	  if (pwm_align == 1)
	    current_contact_position = pwm_position + MAX((*(n[0])).width, (*(n[1])).width) - (*(n[contacts_from_pwm])).match.position;
	  /* printf("\nPosition : %i", current_contact_position); */
	  contacts_defined = 1;
	  top_backbone_contacts = 0;
	  top_base_contacts = 0;
	  bottom_backbone_contacts = 0;
	  bottom_base_contacts = 0;
	  if (current_contact_position < 0 || current_contact_position > (*(n[contacts_from_pwm])).width)
	    contacts_defined = 0;
	  if (contacts_defined == 1) {
	    top_backbone_contacts =
		(*(n[contacts_from_pwm])).fraction[4][current_contact_position] + (*(n[contacts_from_pwm])).fraction[5][current_contact_position] +
		(*(n[contacts_from_pwm])).fraction[6][current_contact_position] + (*(n[contacts_from_pwm])).fraction[7][current_contact_position];
	    top_base_contacts =
		(*(n[contacts_from_pwm])).fraction[8][current_contact_position] + (*(n[contacts_from_pwm])).fraction[9][current_contact_position];
	  }
	  /* TOP BACKBONE ELLIPSE */
	  fprintf(outfile, " <ellipse cx=\"%i\" cy=\"%i\" rx=\"%i\" ry=\"%i\" style=\"opacity:0.5;fill:%s;stroke-width:%i;stroke:%s\"/>\n",
		  pwm_position * 20, top_position + 90, 8, 12, colors[color_start - (top_backbone_contacts == 0)], 1, "saddlebrown");
	  /* CONNECTOR LINES */
	  fprintf(outfile, " <polyline points =\"%i,%i %i,%i %i,%i\" style=\"opacity:1;fill:%s;stroke-width:%i;stroke:%s\"/>\n",
		  pwm_position * 20 - 5 + 10 * (pwm_position == 0), top_position + 110, pwm_position * 20, top_position + 103, pwm_position * 20 + 5,
		  top_position + 110, "none", 1, "gray");
	  /* TOP BACKBONE DOTS */
	  for (circles = 0; circles < top_backbone_contacts; circles++) {
	    fprintf(outfile, " <circle cx=\"%i\" cy=\"%.2f\" r=\"%i\" style=\"fill:%s;stroke-width:%.2f;stroke:%s\"/>\n",
		    pwm_position * 20 + dot_offset, top_position + 84 + (double)circles / top_backbone_contacts * 20 + dot_offset, 2,
		    colors[color_start + 4 -
			   (circles <
			    (*(n[contacts_from_pwm])).fraction[4][current_contact_position] +
			    (*(n[contacts_from_pwm])).fraction[6][current_contact_position])], 0.5, "black");
	  }
	  /* TOP BASE RECTANGLE */
	  fprintf(outfile, " <rect x=\"%i\" y=\"%i\" width=\"%i\" height=\"%i\" style=\"opacity:0.5;fill:%s;stroke-width:%i;stroke:%s\"/>\n",
		  pwm_position * 20 + 2, top_position + 110, 16, 30, colors[color_start - (top_base_contacts == 0)], 1, "black");
	  /* TOP BASE DOTS */
	  for (circles = 0; circles < top_base_contacts; circles++) {
	    fprintf(outfile, " <circle cx=\"%i\" cy=\"%.2f\" r=\"%i\" style=\"fill:%s;stroke-width:%.2f;stroke:%s\"/>\n",
		    pwm_position * 20 + 10 + dot_offset, top_position + 118 + (double)circles / top_base_contacts * 20 + dot_offset, 2,
		    colors[color_start + 4 - (circles < (*(n[contacts_from_pwm])).fraction[8][current_contact_position])], 0.5, "black");
	  }

	  if (contacts_defined == 1) {
	    bottom_backbone_contacts =
		(*(n[contacts_from_pwm])).fraction[12][current_contact_position] + (*(n[contacts_from_pwm])).fraction[13][current_contact_position] +
		(*(n[contacts_from_pwm])).fraction[14][current_contact_position] + (*(n[contacts_from_pwm])).fraction[15][current_contact_position];
	    bottom_base_contacts =
		(*(n[contacts_from_pwm])).fraction[10][current_contact_position] + (*(n[contacts_from_pwm])).fraction[11][current_contact_position];
	  }
	  /* BOTTOM BACKBONE ELLIPSE */
	  fprintf(outfile, " <ellipse cx=\"%i\" cy=\"%i\" rx=\"%i\" ry=\"%i\" style=\"opacity:0.5;fill:%s;stroke-width:%i;stroke:%s\"/>\n",
		  pwm_position * 20 + 20, top_position + 195, 8, 12, colors[color_start - (bottom_backbone_contacts == 0)], 1, "saddlebrown");
	  /* CONNECTOR LINES */
	  fprintf(outfile, " <polyline points =\"%i,%i %i,%i %i,%i\" style=\"opacity:1;fill:%s;stroke-width:%i;stroke:%s\"/>\n",
		  pwm_position * 20 + 15, top_position + 175, pwm_position * 20 + 20, top_position + 182,
		  pwm_position * 20 + 25 - 10 * (pwm_position == (*(n[current_pwm])).width - 1), top_position + 175, "none", 1, "gray");
	  /* BOTTOM BACKBONE DOTS */
	  for (circles = 0; circles < bottom_backbone_contacts; circles++) {
	    fprintf(outfile, " <circle cx=\"%i\" cy=\"%.2f\" r=\"%i\" style=\"fill:%s;stroke-width:%.2f;stroke:%s\"/>\n",
		    pwm_position * 20 + 20 + dot_offset, top_position + 201 - (double)circles / bottom_backbone_contacts * 20 - dot_offset, 2,
		    colors[color_start + 4 -
			   (circles <
			    (*(n[contacts_from_pwm])).fraction[15][current_contact_position] +
			    (*(n[contacts_from_pwm])).fraction[13][current_contact_position])], 0.5, "black");
	  }
	  /* BOTTOM BASE RECTANGLE */
	  fprintf(outfile, " <rect x=\"%i\" y=\"%i\" width=\"%i\" height=\"%i\" style=\"opacity:0.5;fill:%s;stroke-width:%i;stroke:%s\"/>\n",
		  pwm_position * 20 + 2, top_position + 145, 16, 30, colors[color_start - (bottom_base_contacts == 0)], 1, "black");
	  /* BOTTOM BASE DOTS */
	  for (circles = 0; circles < bottom_base_contacts; circles++) {
	    fprintf(outfile, " <circle cx=\"%i\" cy=\"%.2f\" r=\"%i\" style=\"fill:%s;stroke-width:%.2f;stroke:%s\"/>\n",
		    pwm_position * 20 + 10 + dot_offset, top_position + 167 - (double)circles / bottom_base_contacts * 20 - dot_offset, 2,
		    colors[color_start + 4 - (circles < (*(n[contacts_from_pwm])).fraction[11][current_contact_position])], 0.5, "black");
	  }
	}
      }

/* PRINTS REPEAT TILES */
      if (align_matches == 1 && current_pwm == 2) {
/* REPEAT */
	fprintf(outfile,
		"<rect x=\"%i\" y=\"%.2f\" width=\"%i\" height=\"%.2f\" style=\"fill:cornflowerblue;stroke:none;stroke-width:2;opacity:0.5\" />",
		(pwm_position + 5) * 20,
		top_position + 70 -
		0.5 * ((double)all_hits_align_scores[1].score[pwm_position][0]) / ((double)all_hits_align_scores[1].count[pwm_position][0] + 0.01),
		20, ((double)all_hits_align_scores[1].score[pwm_position][0]) / ((double)all_hits_align_scores[1].count[pwm_position][0] + 0.01));
/* INVERTED REPEAT */
	fprintf(outfile,
		"<rect x=\"%i\" y=\"%.2f\" width=\"%i\" height=\"%.2f\" style=\"fill:midnightblue;stroke:none;stroke-width:2;opacity:0.5\" />",
		(pwm_position + 5) * 20,
		top_position + 80 -
		0.5 * ((double)all_hits_align_scores[1].score[pwm_position][1]) / ((double)all_hits_align_scores[1].count[pwm_position][1] + 0.01),
		20, ((double)all_hits_align_scores[1].score[pwm_position][1]) / ((double)all_hits_align_scores[1].count[pwm_position][1] + 0.01));

      }
    }
    /* CORE LINE BY JARKKO TOIVONEN */
    int w = (*(n[current_pwm])).width;       // width of the complete pwm (with possible flanks)
    int total_core_length = core_length1 + core_length2 + core_distance;
    int L = (w + total_core_length)/2;             // width of the SELEX sequence
    if (core_length1 > 0) {
      if (core_length2 == 0)
	draw_box(outfile, offset, L-core_length1, L);
      else {
	int start = L-total_core_length;
	if (core_distance < 0) {
	  draw_box(outfile, offset, start, L);
	  draw_vertical_line(outfile, offset, start + core_length1 + core_distance);
	  draw_vertical_line(outfile, offset, start + core_length1);
	}
	else {
	  draw_box(outfile, offset, start, start + core_length1);            // first motif box
	  draw_box(outfile, offset, start + core_length1 + core_distance, L);  // second motif box
	}
      }
    }

    if ((*(n[current_pwm])).negative_values_allowed == 1) {	/* PRINTS OUT POLYLINE */
      fprintf(outfile, "<polyline points =\"");
      fprintf(outfile, "%i,%i %i,%i \" ", offset, 50, pwm_position * 20 + offset, 50);
      fprintf(outfile, "fill = \"none\" stroke = \"black\" stroke-width = \"0.5\"/>");
    }

/* PRINTS OUT NAME OF LOGO AND OTHER DATA */
    if (noname == 0) {
      fprintf(outfile, "<text  x=\"%i\" y=\"%i\" fill = \"black\" stroke = \"%s\" font-size=\"30\" font-family = \"%s\" >%s",
	      pwm_position * 20 + 20 + Nlength * 10 * (current_pwm == 0), top_position + 30, colors[4 - 2 * warning + 3 * (warning == 2)], font,
	      (*(n[current_pwm])).name);
      if (current_pwm == 0 || current_pwm == 4)
	fprintf(outfile, " lambda=%.3f", lambda);
      fprintf(outfile, "</text>\n");
    }

/* PRINTS ALIGNMENT BOX IF TWO LOGOS ARE ALIGNED */
    if (cm == '\0' && contacts == 0) {
      for (start_rectangle = -1, pwm_position = 0; pwm_position < (*(n[current_pwm])).width; pwm_position++) {
/* printf("\nBase %i score = %.2f", pwm_position, (*(n[current_pwm])).position_score[pwm_position]); */
	if ((*(n[current_pwm])).position_score[pwm_position] > rectangle_cutoff && start_rectangle == -1) {
	  start_rectangle = pwm_position;	/* printf("\nStarts rectangle at %i", pwm_position); */
	  continue;
	}
	if (start_rectangle != -1
	    && ((*(n[current_pwm])).position_score[pwm_position] < rectangle_cutoff || pwm_position == (*(n[current_pwm])).width)) {
/* printf("\nEnds rectangle at %i", pwm_position); */
	  fprintf(outfile,
		  "<rect x=\"%i\" y=\"%i\" width=\"%i\" height=\"%i\" style=\"fill:cornflowerblue;stroke:none;stroke-width:2;opacity:0.1\" />",
		  start_rectangle * 20 - 1, top_position - 18 - 95 * (current_pwm == 2), (pwm_position - start_rectangle) * 20, 187);
	  fprintf(outfile, "<rect x=\"%i\" y=\"%i\" width=\"%i\" height=\"%i\" style=\"fill:none;stroke:midnightblue;stroke-width:2;opacity:1\" />",
		  start_rectangle * 20 - 1, top_position - 18 - 95 * (current_pwm == 2), (pwm_position - start_rectangle) * 20, 187);
	  start_rectangle = -1;
	}
      }
      if (start_rectangle != -1) {
	fprintf(outfile, "<rect x=\"%i\" y=\"%i\" width=\"%i\" height=\"%i\" style=\"fill:cornflowerblue;stroke:none;stroke-width:2;opacity:0.1\" />",
		start_rectangle * 20 - 1, top_position - 18 - 95 * (current_pwm == 2), (pwm_position - start_rectangle) * 20, 187);
	fprintf(outfile, "<rect x=\"%i\" y=\"%i\" width=\"%i\" height=\"%i\" style=\"fill:none;stroke:midnightblue;stroke-width:2;opacity:1\" />",
		start_rectangle * 20 - 1, top_position - 18 - 95 * (current_pwm == 2), (pwm_position - start_rectangle) * 20, 187);
      }
    }

    fprintf(outfile, "</g>\n");
  }

  free(tf_kmer_values);
  free(colors);
  free(order);

  fprintf(outfile, "</svg>\n");
  fclose(outfile);
  return (0);
}

/* SUBROUTINE THAT SUBTRACTS NORMALIZED c FROM NORMALIZED PWM n */
short int Subtract_normalized_pwm(struct normalized_pwm *n, struct normalized_pwm *c, short int offset)
{
  short int counter;
  short int position;
  double total_nucleotides = 0;
  if ((*c).width - offset < (*n).width)
    (*n).width = (*c).width - offset;
  for (position = 0; position < (*n).width; position++) {
    for (counter = 0; counter < 4; counter++)
      (*n).fraction[counter][position] -= (*c).fraction[counter][position + offset];
  }
  return (0);
}


/* SUBROUTINE THAT RENORMALIZES NORMALIZED PWM (ROWS IN EACH COLUMN ADD TO 1) */
short int Normalize_pwm(struct normalized_pwm *n)
{
  short int counter;
  short int position;
  double total_nucleotides = 0;
  double normalized_value = 0;
  for (position = 0; position < (*n).width; position++) {
    for (counter = 0, total_nucleotides = 0; counter < 4; counter++) {
      if ((*n).fraction[counter][position] > 0)
	total_nucleotides += (*n).fraction[counter][position];
      else if ((*n).negative_values_allowed == 1)
	total_nucleotides += -(*n).fraction[counter][position];
    }

    for (counter = 0; counter < 4; counter++) {
      normalized_value = ((double)(*n).fraction[counter][position]) / total_nucleotides;
      if ((normalized_value < 0) && ((*n).negative_values_allowed == 0))
	normalized_value = 0;
      (*n).fraction[counter][position] = normalized_value;
    }
  }
  return (0);
}

/* SUBROUTINE THAT LOADS A PWM AND NORMALIZES IT */
short int Load_pwm(struct normalized_pwm *p, char *filename, short int normalize)
{
  long int counter;
  char text1;
  short int line = 0;
  short int pwm_position = 0;
  char *current_string;
  current_string = malloc(200);
  FILE *pwmfile;
  if ((pwmfile = fopen(filename, "r")) == '\0') {
    printf("\nNo File: %s", filename);
    exit(2);
  }
  for (line = 0; line <= 3 + contacts * 12;) {
    for (counter = 0; counter < 30; counter++) {
      text1 = getc(pwmfile);
      if (text1 == EOF || text1 == '\n' || text1 == '\t') {
	current_string[counter] = '\0';
	if (counter > 0
	    && (current_string[0] == '0' || current_string[0] == '1' || current_string[0] == '2' || current_string[0] == '3'
		|| current_string[0] == '4' || current_string[0] == '5' || current_string[0] == '6' || current_string[0] == '7'
		|| current_string[0] == '8' || current_string[0] == '9' || current_string[0] == ' ' || current_string[0] == '-')) {
	  (*p).fraction[line][pwm_position] = atof(current_string);
	  /* printf("\n%f", (*p).fraction[line][pwm_position]); */
	  pwm_position++;
	}
	if (text1 == '\n' || text1 == EOF) {
	  (*p).width = pwm_position;
	  line++;
	  pwm_position = 0;
	}
	break;
      }
      current_string[counter] = text1;
      /* printf ("%c", text1); */
    }
  }
  free(current_string);
  if (normalize == 1)
    Normalize_pwm(p);
// for (line = 0; line < 4 + contacts * 12; line++) {printf("\n"); for (pwm_position = 0; pwm_position < (*p).width; pwm_position++) printf("\t%f", (*p).fraction[line][pwm_position]);}
  if (text1 == EOF && line != 3)
    return (1);
  else
    return (0);
}


struct alignment {
  short int strand;
  short int position;
};


/* ########################## MAIN PROGRAM ############################# */
/* ########################## MAIN PROGRAM ############################# */
/* ########################## MAIN PROGRAM ############################# */
/* ########################## MAIN PROGRAM ############################# */
/* ########################## MAIN PROGRAM ############################# */
int main(int argc, char *argv[])
{

  Nlength++;

  time_t t0, t1;
  t0 = time('\0');

  double *error_values;
  error_values = malloc(sizeof(double) * 10 + 5);
  error_values[0] = 0;

  /* STRINGS FOR DIFFERENT ORIENTATIONS */
  char **orientation_string;
  orientation_string = malloc(sizeof(char *) * 3 + 5);
  orientation_string[0] = "HT";
  orientation_string[1] = "HH";
  orientation_string[2] = "TT";

  char **repeat_report;
  repeat_report = malloc(sizeof(char *) * 6 + 5);
  repeat_report[0] = "--";
  repeat_report[1] = "Odd_Repeat";
  repeat_report[2] = "Even_Repeat";
  repeat_report[3] = "Odd_Palindrome";
  repeat_report[4] = "Even_Palindrome";

  long int *last_kmer;
  last_kmer = malloc(sizeof(long int) * 100 + 5);

  double pvalue;
  short int shortcounter;
  long int counter;
  long int counter2;
  long int counter3;

  /* clears pvalue cache */
  for (counter = 0; counter < 1000; counter++)
    for (counter2 = 0; counter2 < 1000; counter2++)
      pvalue_cache[counter][counter2] = 2;

  FILE *open_file;
  char **file_name;
  file_name = malloc(sizeof(char *) * (number_of_files + 1) + 5);
  for (counter = 0; counter < number_of_files; counter++) {
    file_name[counter] = malloc(1000);
    strcpy(file_name[counter], "no file");
  }
  short int file_number = 0;
  double *number_of_sequences_analyzed;
  number_of_sequences_analyzed = malloc(sizeof(double) * (number_of_files + 1) + 5);
  double *number_of_sequences_with_hits;
  number_of_sequences_with_hits = malloc(sizeof(double) * (number_of_files + 1) + 5);
  double *number_of_sequences_with_no_hits;
  number_of_sequences_with_no_hits = malloc(sizeof(double) * (number_of_files + 1) + 5);
  for (counter = 0; counter < number_of_files; counter++) {
    number_of_sequences_analyzed[counter] = 0;
    number_of_sequences_with_hits[counter] = 0;
    number_of_sequences_with_no_hits[counter] = 0;
  }
  short int firstnoncommandposition = 0;
  short int print_counts = 0;
  short int extendedoutput = 0;
  short int remove_non_unique = 0;
  short int print_input_sequences = 0;
  short int print_p_values = 0;
  short int output_all_gap_lengths_in_one_line = 0;
  char *linecommand;
  short int too_long_kmer = 8;
  short int shortest_kmer = 8;

  short int kmer_length;
  long int line = 0;
  long int current_kmer;
  short int current_kmer_length;
  short int current_gap_position;
  short int current_gap_length;

  short int *kmermatch_position;
  kmermatch_position = malloc(sizeof(short int) * 2 + 5);
  long int *observed_count_p;
  char *plotfilename;
  plotfilename = malloc(1000);
  seed_story = malloc(10000 * sizeof(char) + 5);
  strcpy(seed_story, "\n");
  char *current_sequence = malloc(1000);
  strcpy(current_sequence, "INITIALVALUE");
  char *tempstring = malloc(10000);
  __uint128_t current_sequence_value_ULL;
  __uint128_t forward_sequence_value_ULL;
  __uint128_t deleted_sequence_value_ULL;
  __uint128_t position_value;
  char *forward_lc = "acgt";
  char *dnareverse = "TGCA";
  signed short int position;
  short int start_position;
  short int end_position;
  short int nucleotide_value;
  short int kmer_count = 1;
  long int kmer_incidence;
  long int kmer2_incidence;
  double max_incidence;
  long int total_incidence;
  float current_fold_expected;
  short int nucleotide;
  short int max_nucleotide;
  long int max_kmer_count;
  long int test_kmer;
  long int test_kmer_count;
  short int iupac_bits;
  short int represents_n_nucleotides;
  long int represents_n_kmers;
  long int top_normalized_count;
  double kmer_length_difference_cutoff = 0.35;
  double iupac_cutoff = 0.25;
  short int position2;
  short int localmaxes = 0;

  __uint128_t query_sequence_value = 2 * 64 + 3 * 4;	/* SEQUENCE: GATA */
  short int query_sequence_length = 4;

  /* FORWARD AND REVERSE NUCLEOTIDE COUNT PWM STRUCTURES */
  struct count_pwm *nc;
  nc = malloc(sizeof(struct count_pwm) * 3 + 5);
  count_pwm_init(&nc[0], "EMPTY", max_Nlength, 0);
  count_pwm_init(&nc[1], "EMPTY", max_Nlength, 0);
  short int print_nucleotides = 0;
  short int print_frequencies = 0;
  char text1;
  long int charcounter;
  short int deletion_size = 1;
  short int Nmer_position;
  short int too_long_nmer_position;

  double minimum_kmer_count = 1;
  short int strand;

  short int non_unique_flag = 0;
  short int current_sorted_list = 0;
  short number_of_unordered_sequences = 0;
  struct sequence_incidence_table *index_position;
  struct sequence_incidence_table **sorted_list;
  sorted_list = malloc(sizeof(struct sequence_incidence_table *) * 2 + 5);
  struct sequence_incidence_table *unordered_list;
  unordered_list = malloc(4000);
  struct sequence_incidence_table **sorted_index;
  sorted_index = malloc(sizeof(struct sequence_incidence_table *) * 4098 + 10);
  short int current_sequence_contains_match = 0;

  long int current_number_of_different_sequences = 0;

  long int matches_to_filter = 0;

  short int eof_reached;

  double **flank_kmer_expected_count;
  long int *palindromic_hits;
  palindromic_hits = malloc(number_of_files * sizeof(long int) + 5);
  flank_kmer_count = malloc(number_of_files * sizeof(long int **) + 5);
  flank_kmer_expected_count = malloc(number_of_files * sizeof(double **) + 5);
  for (file_number = 0; file_number < number_of_files; file_number++) {
    palindromic_hits[file_number] = 0;
    flank_kmer_count[file_number] = malloc(256 * sizeof(long int) + 5);
    flank_kmer_expected_count[file_number] = malloc(256 * sizeof(double) + 5);
    for (counter = 0; counter < 256; counter++) {
      flank_kmer_count[file_number][counter] = 0;
      flank_kmer_expected_count[file_number][counter] = 0;
    }
  }
  long int *match_orientations = '\0';
  short int match_length = 0;

  /* VARIABLES AND STRUCTURES FOR HIT POSITION DATA GENERATION */
  struct hit_position_matrix hit_position;
  hit_position_matrix_init(&hit_position, "hit positions", max_Nlength, 0);
  short int correct_position = 1;
  char limit_hits_to_strand = 'N';
  short int limit_hits_to_position = 0;
  char *exclude_strand;
  exclude_strand = malloc(sizeof(char *) * max_Nlength + 5);
  short int *exclude_position;
  exclude_position = malloc(sizeof(short int *) * max_Nlength + 5);
  for (counter = 0; counter < max_Nlength; counter++) {
    exclude_position[counter] = 0;
    exclude_strand[counter] = 'N';
  }

  /* FLAGS */
  short int xyplot = 0;
  long int max_counts = 0;
  short int max_spacing = 0;
  short int max_orientation = 0;
  short int flank_with_pwm = 0;
  short int even_background = 0;
  short int only_palindromes = 0;
  short int palindrome_correction = 1;
  short int information_content_output = 0;
  short int count_also_spaced_kmers = 1;
  short int dinucleotide_properties = 0;
  short int searchstring_is_palindrome = 0;
  short int flank = 0;
  short int kmer_table = 0;
  short int multinomial = 0;
  short int complex_background = 1;
  short int center_spaced_kmer_table = 0;
  short int svg_only = 0;
  short int match_filter = 0;
  short int modeldistance = 0;
  short int difference_logo = 0;
  short int negative_values_allowed = 0;
  short int offset = 0;
  short int lowpaircount = 0;
  short int expected_observed_plot = 0;
  short int number_of_heatmap_rows = 20;
  short int iterate_seed = 0;
  short int remember_iterate_seed = 0;
  short int end_trim_seed = 0;
  short int same_seed_size = 0;
  short int max_seed_size = 0;
  short int auto_seed = 1;
  short int seed_from_local_max_number = 1;

  long int *signal_kmer_count_p;
  long int *background_kmer_count_p;
  long int *kmer_count_p_memory;

  char *user_specified_output_file;
  user_specified_output_file = malloc(1000);
  user_specified_output_file[0] = '\0';
  short int **max_align_score;
  max_align_score = malloc(sizeof(short int *) * 3 + 5);
  max_align_score[0] = malloc(sizeof(short int) * 3 + 5);
  max_align_score[1] = malloc(sizeof(short int) * 3 + 5);
  max_align_score[0][0] = 0;
  max_align_score[0][1] = 0;
  max_align_score[1][0] = 0;
  max_align_score[1][1] = 0;
  char *forward;
  forward = dnaforward;
  char *nucleotide_bitiupac = dna_bitiupac;
  char *nucleotide_iupac = dna_iupac;

  COMMAND = malloc(5000);
  strcpy(COMMAND, "");
  for (counter = 0; counter < argc; counter++) {
    strcat(COMMAND, argv[counter]);
    strcat(COMMAND, " ");
  }

/* COMMAND PARSER */
/* PARSES OPTIONS */
  if (argc >= 2 + firstnoncommandposition) {
    for (;; firstnoncommandposition++) {
      linecommand = argv[1 + firstnoncommandposition];

      if (linecommand[0] == '-' && linecommand[1] != '\0') {
	if (strcmp(linecommand, "--pwm") == 0) {
	  flank_with_pwm = 1;
	  flank = 1;
	}
	if (strcmp(linecommand, "--f") == 0) {
	  flank = 1;
	  kmer_count = 0;
	}
	if (strcmp(linecommand, "--logo") == 0)
	  svg_only = 1;
	if (strcmp(linecommand, "--difflogo") == 0) {
	  svg_only = 1;
	  difference_logo = 1;
	  paths = 1;
	}
	if (strcmp(linecommand, "--dist") == 0)
	  modeldistance = 1;
	if (strcmp(linecommand, "--pwmalign") == 0)
	  pwm_align = 1;
	if (strcmp(linecommand, "-editdist") == 0) {
	  lowpaircount = 1;
	  count_also_spaced_kmers = 1;
	};
	if (strcmp(linecommand, "-paths") == 0)
	  paths = 1;
	if (strcmp(linecommand, "-barcodelogo") == 0)
	  barcodelogo = 1;
	if (strcmp(linecommand, "-heightscaledbars") == 0)
	  scale_bars = 1;
	if (strcmp(linecommand, "-maxheightscaledbars") == 0) {
	  scale_bars = 1;
	  max_scale_bar = 1;
	}
	if (strcmp(linecommand, "-colorscaledbars") == 0)
	  gray_bars = 1;
	if (strcmp(linecommand, "-label") == 0)
	  barcodelogolabels = 1;
	if (strcmp(linecommand, "-neg") == 0) {
	  negative_values_allowed = 1;
	  paths = 1;
	}

	if (strncmp(linecommand, "-core=", 6) == 0) {
	  char * param = linecommand + 6;
	  if (sscanf(param, "%i,%i,%i", &core_length1, &core_length2, &core_distance) != 3) {
	    core_length1 = atoi(param);
	    core_length2 = 0;
	    core_distance = 0;
	  }
	}
	if (strcmp(linecommand, "-noname") == 0)
	  noname = 1;

	if (strcmp(linecommand, "-rna") == 0) {
	  rna = 1;
	  nucleotide_bitiupac = rna_bitiupac;
	  nucleotide_iupac = rna_iupac;
	  forward = rnaforward;
	}
	if (strcmp(linecommand, "-x") == 0)
	  extendedoutput = 1;
	if (strcmp(linecommand, "-u") == 0)
	  remove_non_unique = 1;
	if (strcmp(linecommand, "-c") == 0)
	  print_counts = 1;
	if (strcmp(linecommand, "-p") == 0)
	  print_p_values = 1;
	if (strcmp(linecommand, "-i") == 0)
	  print_input_sequences = 1;
	if (strcmp(linecommand, "-s") == 0)
	  output_all_gap_lengths_in_one_line = 1;
	if (strcmp(linecommand, "-n") == 0)
	  print_nucleotides = 1;
	if (strcmp(linecommand, "-q") == 0)
	  print_frequencies = 1;
	if (strcmp(linecommand, "-a") == 0)
	  align_matches = 1;
	if (strcmp(linecommand, "-xyplot") == 0)
	  xyplot = 1;
	if (strcmp(linecommand, "-eoplot") == 0)
	  expected_observed_plot = 1;
	if (strcmp(linecommand, "-contacts") == 0)
	  contacts = 1;
	if (strcmp(linecommand, "-printlocalmax") == 0)
	  print_local_max = 1;
	if (strcmp(linecommand, "-nogaps") == 0)
	  count_also_spaced_kmers = 0;
	if (strcmp(linecommand, "-allgaps") == 0)
	  count_also_spaced_kmers = 2;
	if (strcmp(linecommand, "-both") == 0)
	  count_both_instances_of_palindromic_hit = 1;
	if (strcmp(linecommand, "-nocall") == 0)
	  nocall = 1;
	if (strcmp(linecommand, "-CpG") == 0) {
	  methylCGcompare = 1;
	  dinucleotide_properties = 1;
	}
	if (strcmp(linecommand, "-bothifnotequal") == 0)
	  count_unequal_hits_only = 1;
	if (strcmp(linecommand, "-forwardonly") == 0)
	  count_only_forward_instance_of_palindromic_hit = 1;
	if (strcmp(linecommand, "-reverseonly") == 0)
	  count_only_reverse_instance_of_palindromic_hit = 1;
	if (strcmp(linecommand, "-forwardifequal") == 0)
	  prefer_forward_strand = 1;
	if (strcmp(linecommand, "-reverseifequal") == 0)
	  prefer_reverse_strand = 1;
	if (strcmp(linecommand, "-dimer") == 0)
	  only_palindromes = 1;
	if (strcmp(linecommand, "-dinuc") == 0)
	  dinucleotide_properties = 1;
	if (strcmp(linecommand, "-kmer") == 0)
	  kmer_table = 1;
	if (strcmp(linecommand, "-ic") == 0)
	  information_content_output = 1;
	if (strcmp(linecommand, "-e") == 0) {
	  even_background = 1;
	  file_number++;
	}
	if (strcmp(linecommand, "-14N") == 0)
	  Nlength = 15;
	if (strcmp(linecommand, "-26N") == 0)
	  Nlength = 27;
	if (strcmp(linecommand, "-30N") == 0)
	  Nlength = 31;
	if (strcmp(linecommand, "-40N") == 0)
	  Nlength = 41;
	if (strcmp(linecommand, "-mono") == 0)
	  complex_background = 0;
	if (strcmp(linecommand, "-iterate") == 0) {
	  iterate_seed = 1;
	}
	if (strcmp(linecommand, "-iterate-fast") == 0) {
	  iterate_seed = 1;
	  end_trim_seed = 1;
	}
	if (strcmp(linecommand, "-iterate-samesize") == 0) {
	  iterate_seed = 1;
	  same_seed_size = 1;
	}
	if (linecommand[0] == '-' && linecommand[1] == 'o' && linecommand[2] == '=')
	  strcpy(user_specified_output_file, linecommand + 3);

	if (linecommand[0] == '-' && linecommand[1] == 'm' && linecommand[2] == '=')
	  multinomial = atoi(linecommand + 3);

	if (linecommand[0] == '-' && linecommand[1] == 'f' && linecommand[2] == 'k' && linecommand[3] == '=')
	  flank_kmer_pos = atoi(linecommand + 4);

	if (linecommand[0] == '-' && linecommand[1] == 'm' && linecommand[2] == 'f' && linecommand[3] == '=')
	  match_filter = atoi(linecommand + 4);

	if (linecommand[0] == '-' && linecommand[1] == 'k' && linecommand[2] == 'l' && linecommand[3] == '=') {
	  shortest_kmer = atoi(linecommand + 4);
	  too_long_kmer = shortest_kmer + 1;
	}

	if (linecommand[0] == '-' && linecommand[1] == 'l' && linecommand[2] == 'o' && linecommand[7] == '=') {
	  kmer_length_difference_cutoff = atof(linecommand + 8);
	}

	if (linecommand[0] == '-' && linecommand[1] == 'i' && linecommand[2] == 'u' && linecommand[6] == '=') {
	  iupac_cutoff = atof(linecommand + 7);
	}

	if (linecommand[0] == '-' && linecommand[1] == 'h' && linecommand[2] == 'r' && linecommand[6] == '=') {
	  number_of_heatmap_rows = atof(linecommand + 7);
	}

	if (linecommand[0] == '-' && linecommand[1] == 'e' && linecommand[2] == 'd' && linecommand[9] == '=') {
	  lowpaircount = 1;
	  local_max_min_percent = ((double)atoi(linecommand + 10)) / 100;
	}

	if (linecommand[0] == '-' && linecommand[1] == 'l' && linecommand[2] == 'i' && linecommand[3] == 'm' && linecommand[4] == '=') {
	  limit_hits_to_position = atoi(linecommand + 5);
	  limit_hits_to_strand = linecommand[strlen(linecommand) - 1];
	}

	if (linecommand[0] == '-' && linecommand[1] == 'e' && linecommand[2] == 'x' && linecommand[3] == 'c' && linecommand[4] == '=') {
	  for (counter = 0, linecommand += 5; counter < max_Nlength; counter++, linecommand++) {
	    exclude_position[counter] = atoi(linecommand);
	    for (; *linecommand != ',' && *linecommand != '\0'; linecommand++) ;
	    exclude_strand[counter] = *(linecommand - 1);
	    if (*linecommand == '\0')
	      break;
	  }
	}

	/* if (linecommand[0] == '-' && linecommand[1] == 'm' && linecommand[2] == 'a' && linecommand[3] == 't' && linecommand[6] == '=') { */
	/*   char *halfsite1; */
	/*   char *halfsite2; */
	/*   halfsite1 = malloc(1000); */
	/*   halfsite2 = malloc(1000); */
	/*   for (counter = 7; counter < strlen(linecommand), linecommand[counter] != ','; counter++) ; */
	/*   strcpy(halfsite2, linecommand + counter + 1); */
	/*   strcpy(halfsite1, linecommand + 7); */
	/*   halfsite1[counter - 7] = '\0'; */
	/*   //printf("\nHalfsites %s,%s, %i\n", halfsite1, halfsite2, counter); */
	/*   for (query_sequence_value = 0, position = 0, position_value = pow(4, strlen(halfsite1) - 1); position < strlen(halfsite1); */
	/*        position++, position_value /= 4) { */
	/*     for (nucleotide_value = 0; nucleotide_value < 4 && halfsite1[position] != forward[nucleotide_value]; nucleotide_value++) ; */
	/*     if (nucleotide_value == 4) { */
	/*       printf("\nERROR IN MATCH SEQUENCE\n"); */
	/*       exit(1); */
	/*     } */
	/*     query_sequence_value += position_value * nucleotide_value; */
	/*   } */
	/*   long int halfsite1_value = query_sequence_value; */
	/*   for (query_sequence_value = 0, position = 0, position_value = pow(4, strlen(halfsite2) - 1); position < strlen(halfsite2); */
	/*        position++, position_value /= 4) { */
	/*     for (nucleotide_value = 0; nucleotide_value < 4 && halfsite2[position] != forward[nucleotide_value]; nucleotide_value++) ; */
	/*     if (nucleotide_value == 4) { */
	/*       printf("\nERROR IN MATCH SEQUENCE\n"); */
	/*       exit(1); */
	/*     } */
	/*     query_sequence_value += position_value * nucleotide_value; */
	/*   } */
	/*   long int halfsite2_value = query_sequence_value; */

	/*   match_orientations = malloc(sizeof(long int) * 10 + 5); */
	/*   short int shift1 = strlen(halfsite1) * 2; */
	/*   short int shift2 = strlen(halfsite2) * 2; */

	/*   /\* HETERODIMER *\/ */
	/*   match_orientations[0] = halfsite1_value << shift1 | halfsite2_value; */
	/*   match_orientations[1] = halfsite2_value << shift2 | halfsite1_value; */
	/*   match_orientations[2] = Reverse_complement_sequence_value_li(halfsite1_value, strlen(halfsite1)) << shift1 | halfsite2_value; */
	/*   match_orientations[3] = halfsite1_value << shift1 | Reverse_complement_sequence_value(halfsite2_value, strlen(halfsite2)); */

	/*   /\* HOMODIMER1 *\/ */
	/*   match_orientations[4] = halfsite1_value << shift1 | halfsite1_value; */
	/*   match_orientations[5] = Reverse_complement_sequence_value(halfsite1_value, strlen(halfsite1)) << shift1 | halfsite1_value; */
	/*   match_orientations[6] = halfsite1_value << shift1 | Reverse_complement_sequence_value(halfsite1_value, strlen(halfsite1)); */

	/*   /\* HOMODIMER2 *\/ */
	/*   match_orientations[7] = halfsite2_value << shift2 | halfsite2_value; */
	/*   match_orientations[8] = Reverse_complement_sequence_value(halfsite2_value, strlen(halfsite2)) << shift2 | halfsite2_value; */
	/*   match_orientations[9] = halfsite2_value << shift2 | Reverse_complement_sequence_value(halfsite2_value, strlen(halfsite2)); */
	/*   match_length = strlen(halfsite1) + strlen(halfsite2); */
	/* } */

      } else
	break;
    }
  }

/* INITIALIZES VARIABLES WHOSE SIZE OR USE IS DEPENDENT ON ARGUMENTS */
  __uint128_t left_position_value = 1;
  left_position_value <<= ((Nlength - 2) * 2);

  long int ***align_score_histogram;
  align_score_histogram = malloc(sizeof(long int *) * 3 + 5);
  align_score_histogram[0] = malloc(sizeof(long int *) * 3 + 5);
  align_score_histogram[1] = malloc(sizeof(long int *) * 3 + 5);
  align_score_histogram[0][0] = malloc(sizeof(long int) * Nlength + 5);
  align_score_histogram[0][1] = malloc(sizeof(long int) * Nlength + 5);
  align_score_histogram[1][0] = malloc(sizeof(long int) * Nlength + 5);
  align_score_histogram[1][1] = malloc(sizeof(long int) * Nlength + 5);
  for (counter = 0; counter < Nlength; counter++) {
    align_score_histogram[0][0][counter] = 0;
    align_score_histogram[1][0][counter] = 0;
    align_score_histogram[0][1][counter] = 0;
    align_score_histogram[1][1][counter] = 0;
  }

/* PRIMARY mask_ULL FOR EACH SEPARATE NUCLEOTIDE */
  for (mask_ULL[1][0] = 3, counter = 0; counter < Nlength; counter++)
    mask_ULL[1][counter + 1] = mask_ULL[1][counter] << 2;

  /* FLANK TOOL VARIABLES */
  char *searchstring;
  searchstring = malloc(1000);
  strcpy(searchstring, "-");
  char *valuestring;
  valuestring = malloc(1000);

  short int *number_of_matches;
  number_of_matches = malloc(sizeof(short int) * 200 + 5);
  number_of_matches[0] = 0;
  number_of_matches[1] = 0;
  struct match *match;
  match = malloc(sizeof(struct match) * 3 + 5);
  match_init(&match[0], Nlength);
  match_init(&match[1], Nlength);
  match_init(&match[2], Nlength);
  struct match *dinucmatch;
  dinucmatch = malloc(sizeof(struct match) * 3 + 5);
  match_init(&dinucmatch[0], Nlength);
  match_init(&dinucmatch[1], Nlength);
  match_init(&dinucmatch[2], Nlength);
  struct count_connecting_matrix cm;
  count_connecting_matrix_init(&cm, "all hits connecting matrix", Nlength, 0);
  struct count_connecting_matrix **cms_to_heatmap;
  cms_to_heatmap = malloc(sizeof(struct count_connecting_matrix *) * 4 + 5);
  cms_to_heatmap[0] = &cm;

  short int iupac_query = 0;

  struct count_connecting_matrix two_hits_connecting_matrix;
  count_connecting_matrix_init(&two_hits_connecting_matrix, "two hit connecting matrix", Nlength, 0);
  struct count_connecting_matrix two_hits_fold_connecting_matrix;
  count_connecting_matrix_init(&two_hits_fold_connecting_matrix, "fold two hit connecting matrix", Nlength, 0);

  struct count_pwm *all_hits_pwm;
  all_hits_pwm = malloc(sizeof(struct count_pwm) * 3 + 5);
  struct count_pwm *one_hit_pwm;
  one_hit_pwm = malloc(sizeof(struct count_pwm) * 3 + 5);
  struct count_pwm *one_hit_exponential;
  one_hit_exponential = malloc(sizeof(struct count_pwm) * 3 + 5);
  struct count_pwm *all_hits_exponential;
  all_hits_exponential = malloc(sizeof(struct count_pwm) * 3 + 5);
  short int reverse_strand_position;
  struct count_pwm ***two_hits_pwm;
  two_hits_pwm = malloc(sizeof(struct count_pwm *) * 3 + 5);
  two_hits_pwm[0] = malloc(sizeof(struct count_pwm *) * 5 + 5);
  two_hits_pwm[1] = malloc(sizeof(struct count_pwm *) * 5 + 5);
  for (counter = 0; counter < 5; counter++) {
    two_hits_pwm[0][counter] = malloc(sizeof(struct count_pwm) * (Nlength * 2) + 5);
    two_hits_pwm[1][counter] = malloc(sizeof(struct count_pwm) * (Nlength * 2) + 5);
  }

  for (file_number = 0; file_number < number_of_files; file_number++) {
    count_pwm_init(&all_hits_pwm[file_number], "empty", Nlength * 2, 0);
    count_pwm_init(&one_hit_pwm[file_number], "empty", Nlength * 2, 0);
    count_pwm_init(&one_hit_exponential[file_number], "empty", Nlength * 2, 0);
    count_pwm_init(&all_hits_exponential[file_number], "empty", Nlength * 2, 0);
    for (counter2 = 0; counter2 < 4; counter2++)
      for (counter3 = 0; counter3 < Nlength * 2; counter3++)
	count_pwm_init(&two_hits_pwm[file_number][counter2][counter3], "empty", Nlength * 2, 0);
  }

  /* QUERY PWM STRUCTURE */
  struct normalized_pwm qp;
  normalized_pwm_init(&qp, "empty", Nlength * 2, 0);

  /* MATCH FILTER PWM STRUCTURE */
  struct normalized_pwm mfp;
  normalized_pwm_init(&mfp, "empty", Nlength * 2, 0);

  /* DINUC MULTINOMIAL 2 PWM STRUCTURE */
  struct normalized_pwm mmp;
  normalized_pwm_init(&mmp, "empty", Nlength * 2, 0);

  /* BACKGROUND QUERY PWM STRUCTURE */
  struct normalized_pwm bqp;
  normalized_pwm_init(&bqp, "empty", Nlength * 2, 0);

  /* PRINT PWM STRUCTURE */
  struct normalized_pwm p;
  normalized_pwm_init(&p, "empty", Nlength * 2, 0);

  /* NORMALIZED BACKGROUND PWM STRUCTURES */
  struct normalized_pwm background_pwm[2];
  normalized_pwm_init(&background_pwm[0], "empty", Nlength * 2, 0.25);
  normalized_pwm_init(&background_pwm[1], "empty", Nlength * 2, 0.25);

  /* NORMALIZED CONNECTING MATRIX */
  struct normalized_connecting_matrix cp;
  normalized_connecting_matrix_init(&cp, "normalized connecting matrix", Nlength * 2, 0);

  /* ADM STRUCTURES */
  struct adjacent_dinucleotide_model unflanked_adm;
  struct adjacent_dinucleotide_model flanked_adm;
  struct adjacent_dinucleotide_model *flanked_adm_p;

  short int first;
  short int last;
  short int first_match_position;
  short int spacing;
  short int orientation;
  __uint128_t first_sequence_value_ULL;
  char *pwm_name;
  double cut_off = 9;
  double filter_cut_off = 9;
  double multi_2_cut_off = 9;

  /* VARIABLES AND STRUCTURES FOR LOGO GENERATION */
  short int current_logo = 0;
  struct normalized_pwm *np;
  np = malloc(sizeof(struct normalized_pwm) * 200 + 5);
  struct normalized_pwm **np_p;
  np_p = malloc(sizeof(struct normalized_pwm *) * 200 + 5);
  for (counter = 0; counter < 200; counter++) {
    normalized_pwm_init(&np[counter], "empty", Nlength * 2, 0);
    np_p[counter] = &np[counter];
  }
  long int total_number_of_two_hits;
  long int total_possible_spacings;

/* VARIABLES FOR MULTINOMIAL */
  double swap = 0;
  file_number = 0;
  double lambda = 1;
  short int no_insertions = 0;
  short int multinomial_2_dinuc = 0;
  long int background_kmer_count;
  long int signal_kmer_count;
  long int number_of_kmers = 0;

/* VARIABLES FOR INFORMATION CONTENT */
  double background_info = 16;
  double signal_info = 16;
  double background_nonhit_info = -100;
  double signal_nonhit_info = -100;
  double fraction;
  double total_exp;
  double total_obs;

  double current_ic_score;
  double max_ic_score = 0;
  double ic_score = 0;
  double delta_ic = 0;
  short int max_first = 0;
  short int max_second = 0;
  short int end_trim = 15;
  short int max_dinucleotide;
  short int min_dinucleotide;
  short int warning = 0;
  double forward_pos_ic;
  double reverse_pos_ic;
  double min_fold;
  double max_fold;
  double current_fold;

  short int load_pwm = 0;
  short int load_adm = 0;
  short int loaded_pwm_width = 0;

/* double **total_relative_deviation;
total_relative_deviation = malloc(sizeof(double *) * Nlength * 2 + 5);
for(counter = 0; counter < Nlength * 2; counter++) total_relative_deviation[counter] = malloc(sizeof(double) * Nlength * 2 + 5); */
  double **uncentered_correlation;
  uncentered_correlation = malloc(sizeof(double *) * Nlength * 2 + 5);
  for (counter = 0; counter < Nlength * 2; counter++)
    uncentered_correlation[counter] = malloc(sizeof(double) * Nlength * 2 + 5);
  long double sum_x_squared;
  long double sum_y_squared;
  long double sum_xy;

/* PWM ALIGN OPTION ENTIRE MAIN PROGRAM CODE */
/*   if (pwm_align == 1 && argc >= 3 + firstnoncommandposition) { */
/*     short int number_of_pwms = 2; */
/*     short int width = 0; */

/*     strcpy(searchstring, argv[1 + firstnoncommandposition]); */
/*     printf("\n\nMATRIX 1:\t%s", searchstring); */
/*     Load_pwm(&qp, searchstring, 0); */
/*     strcpy(qp.name, searchstring); */
/*     Normalize_pwm(&qp); */
/* // Log_ratio_pwm(&qp); */

/*     strcpy(tempstring, argv[2 + firstnoncommandposition]); */
/*     printf("\n\nMATRIX 2:\t%s", tempstring); */
/*     Load_pwm(&mfp, tempstring, 0); */
/*     strcpy(mfp.name, tempstring); */
/*     Normalize_pwm(&mfp); */
/* // Log_ratio_pwm(&mfp); */

/*     if (argc >= 4 + firstnoncommandposition) { */
/*       number_of_pwms = 3; */
/*       strcpy(valuestring, argv[3 + firstnoncommandposition]); */
/*       printf("\n\nMATRIX 3:\t%s", valuestring); */
/*       Load_pwm(&p, valuestring, 0); */
/*       strcpy(p.name, valuestring); */
/*       Normalize_pwm(&p); */
/* // Log_ratio_pwm(&p);         */
/*     } */

/*     width = qp.width; */
/*     if (number_of_pwms == 3 && p.width > qp.width) */
/*       width = p.width; */

/*         Add_flanks(&mfp, width, 0.25); */

/*     struct normalized_pwm rc_qp; */
/*     normalized_pwm_init(&rc_qp, qp.name, qp.width, 0); */
/*     Reverse_complement_normalized_pwm(&rc_qp, &qp); */

/*     struct normalized_pwm *qp_p; */
/*     qp_p = Align_pwm(&qp, &rc_qp, &mfp); */
/*     mfp.match.position = mfp.match.position + width; */

/*     struct normalized_pwm *p_p; */
/*     if (number_of_pwms == 3) { */
/*       struct normalized_pwm rc_p; */
/*       normalized_pwm_init(&rc_p, p.name, p.width, 0); */
/*       Reverse_complement_normalized_pwm(&rc_p, &p); */

/*       p_p = Align_pwm(&p, &rc_p, &mfp); */

/*     } */

/*     strcat(searchstring, "_"); */
/*     strcat(searchstring, tempstring); */
/*     if (number_of_pwms == 3) { */
/*       strcat(searchstring, "_"); */
/*       strcat(searchstring, valuestring); */
/*     } */
/*     strcat(searchstring, ".svg"); */

/*     Load_pwm(&mfp, tempstring, 1); */

/*     np_p[0] = qp_p; */
/*     np_p[1] = &mfp; */
/*     if (number_of_pwms == 3) { */
/*       if (contacts == 0) */
/* 	np_p[2] = p_p; */
/*       else { */
/* 	np_p[1] = p_p; */
/* 	np_p[2] = &mfp; */
/*       } */
/*       nocall = 2; */
/*     } else */
/*       nocall = 1; */

/*     noname = 1; */

/*     Svg_logo(searchstring, number_of_pwms, np_p, '\0', '\0', '\0', '\0', '\0', '\0', 0, 0); */

/*     printf("\n\nPWM1 best match score to PWM2 is %.2f position is %i strand %i", (*qp_p).match.score, (*qp_p).match.position - width + 1, */
/* 	   (*qp_p).match.strand); */
/*     if (number_of_pwms == 3) { */
/*       printf("\nPWM3 best match score to PWM2 is %.2f position is %i strand %i\n", (*p_p).match.score, (*p_p).match.position - width + 1, */
/* 	     (*p_p).match.strand); */

/*       short int pwm_position; */
/*       struct normalized_pwm *first_match_pwm_p; */
/*       struct normalized_pwm *second_match_pwm_p; */
/*       if ((*p_p).match.position < (*qp_p).match.position) { */
/* 	first_match_pwm_p = p_p; */
/* 	second_match_pwm_p = qp_p; */
/*       } else { */
/* 	first_match_pwm_p = qp_p; */
/* 	second_match_pwm_p = p_p; */
/*       } */

/* /\* PRINTS ALIGNMENT SCORES BY BASE *\/ */
/*       char *overlapcall[] = { "none", "OVERLAP" }; */
/*       char *strengthcall[] = { "strong", "weak" }; */
/*       printf("\nSIDE\tNAME\tINPUT_PWM_POS\tSCORE\tOVERLAP\tALIGNMENT"); */
/*       for (counter = 0; counter < (*first_match_pwm_p).width; counter++) */
/* 	printf("\nleft\t%s\t%i\t%.2f\t%s\t%s\talign_score", (*first_match_pwm_p).name, (*first_match_pwm_p).original_position[counter] + 1, */
/* 	       (*first_match_pwm_p).position_score[counter], */
/* 	       overlapcall[counter + (*first_match_pwm_p).match.position >= (*second_match_pwm_p).match.position], */
/* 	       strengthcall[(*first_match_pwm_p).position_score[counter] < 0.2]); */
/*       for (counter = 0; counter < (*second_match_pwm_p).width; counter++) */
/* 	printf("\nright\t%s\t%i\t%.2f\t%s\t%s\talign_score", (*second_match_pwm_p).name, (*second_match_pwm_p).original_position[counter] + 1, */
/* 	       (*second_match_pwm_p).position_score[counter], */
/* 	       overlapcall[(*second_match_pwm_p).match.position + counter < (*first_match_pwm_p).match.position + (*first_match_pwm_p).width], */
/* 	       strengthcall[(*first_match_pwm_p).position_score[counter] < 0.25]); */

/* /\* PRINTS TOTAL GAP WIDTH *\/ */
/*       printf("\n\nComposite site total gap width: \t%3i", abs((*p_p).match.position - (*qp_p).match.position) - (*first_match_pwm_p).width); */

/*       if (contacts == 1) { */
/* /\* PRINTS CORE GAP WIDTH *\/ */
/* 	for (counter = 0; counter < (*second_match_pwm_p).width; counter++) */
/* 	  if ((*second_match_pwm_p).fraction[8][counter] + (*second_match_pwm_p).fraction[9][counter] + (*second_match_pwm_p).fraction[10][counter] + */
/* 	      (*second_match_pwm_p).fraction[11][counter] != 0) */
/* 	    break; */
/* 	for (counter2 = (*first_match_pwm_p).width - 1; counter2 > 0; counter2--, counter++) */
/* 	  if ((*first_match_pwm_p).fraction[8][counter2] + (*first_match_pwm_p).fraction[9][counter2] + (*first_match_pwm_p).fraction[10][counter2] + */
/* 	      (*first_match_pwm_p).fraction[11][counter2] != 0) */
/* 	    break; */
/* 	printf("\nComposite site core gap width:  \t%3li", */
/* 	       abs((*p_p).match.position - (*qp_p).match.position) - (*first_match_pwm_p).width + counter); */

/* /\* PRINTS BACKBONE GAP WIDTH *\/ */
/* 	for (counter = 0; counter < (*second_match_pwm_p).width; counter++) */
/* 	  if ((*second_match_pwm_p).fraction[4][counter] + (*second_match_pwm_p).fraction[5][counter] + (*second_match_pwm_p).fraction[6][counter] + */
/* 	      (*second_match_pwm_p).fraction[7][counter] + (*second_match_pwm_p).fraction[12][counter] + (*second_match_pwm_p).fraction[13][counter] + */
/* 	      (*second_match_pwm_p).fraction[14][counter] + (*second_match_pwm_p).fraction[15][counter] != 0) */
/* 	    break; */
/* 	for (counter2 = (*first_match_pwm_p).width - 1; counter2 > 0; counter2--, counter++) */
/* 	  if ((*first_match_pwm_p).fraction[4][counter2] + (*first_match_pwm_p).fraction[5][counter2] + (*first_match_pwm_p).fraction[6][counter2] + */
/* 	      (*first_match_pwm_p).fraction[7][counter2] + (*first_match_pwm_p).fraction[12][counter2] + (*first_match_pwm_p).fraction[13][counter2] + */
/* 	      (*first_match_pwm_p).fraction[14][counter2] + (*first_match_pwm_p).fraction[15][counter2] != 0) */
/* 	    break; */
/* 	printf("\nComposite site backbone gap width:  \t%3li", */
/* 	       abs((*p_p).match.position - (*qp_p).match.position) - (*first_match_pwm_p).width + counter); */

/* /\* RELOADS ORIGINAL MIDDLE (COMPOSITE) PWM *\/ */
/* 	Load_pwm(&mfp, tempstring, 0); */
/* /\* ADDS UP CONTACTS *\/ */
/* 	short int position_on_left_pwm; */
/* 	short int position_on_right_pwm; */
/* 	for (nucleotide = 4; nucleotide < 16; nucleotide++) */
/* 	  for (counter = 0; counter < mfp.width; counter++) { */
/* 	    mfp.fraction[nucleotide][counter] = 0; */
/* 	    position_on_left_pwm = counter - (*first_match_pwm_p).match.position + width + 1; */
/* 	    position_on_right_pwm = counter - (*second_match_pwm_p).match.position + width + 1; */
/* 	    if (position_on_left_pwm >= 0 && position_on_left_pwm < (*first_match_pwm_p).width) */
/* 	      mfp.fraction[nucleotide][counter] += (*first_match_pwm_p).fraction[nucleotide][position_on_left_pwm]; */
/* 	    if (position_on_right_pwm >= 0 && position_on_right_pwm < (*second_match_pwm_p).width) */
/* 	      mfp.fraction[nucleotide][counter] += (*second_match_pwm_p).fraction[nucleotide][position_on_right_pwm]; */
/* 	  } */
/* /\* PRINTS COMPOSITE PWM WITH TOTAL CONTACTS *\/ */
/* 	printf("\n\nOriginal composite PWM (2) with contacts from the individual PWMs (1,3)"); */
/* 	for (nucleotide = 0; nucleotide < 16; nucleotide++) { */
/* 	  printf("\n%s_Composite_PWM\t%s", mfp.name, pwm_row_ids[nucleotide]); */
/* 	  for (pwm_position = 0; pwm_position < mfp.width; pwm_position++) */
/* 	    printf("\t%.0f", mfp.fraction[nucleotide][pwm_position]); */
/* 	} */

/* /\* GENERATES CORE-TO-CORE PWM *\/ */
/* 	for (counter = 0; counter < (*first_match_pwm_p).width; counter++) */
/* 	  if ((*first_match_pwm_p).fraction[8][counter] + (*first_match_pwm_p).fraction[9][counter] + (*first_match_pwm_p).fraction[10][counter] + */
/* 	      (*first_match_pwm_p).fraction[11][counter] != 0) */
/* 	    break; */
/* 	for (counter2 = (*second_match_pwm_p).width - 1; counter2 > 0; counter2--) */
/* 	  if ((*second_match_pwm_p).fraction[8][counter2] + (*second_match_pwm_p).fraction[9][counter2] + */
/* 	      (*second_match_pwm_p).fraction[10][counter2] + (*second_match_pwm_p).fraction[11][counter2] != 0) */
/* 	    break; */
/* 	printf("\n\nCore-to-core PWM start position:\t%3li\n", (*first_match_pwm_p).match.position + counter - width + 1); */
/* 	printf("Core-to-core PWM end position:   \t%3li\n", */
/* 	       (*second_match_pwm_p).match.position + (*second_match_pwm_p).width - 1 - counter2 - width); */
/* 	printf("\nCore-to-core PWM"); */

/* 	for (nucleotide = 0; nucleotide < 16; nucleotide++) { */
/* 	  printf("\n%s_CTC_PWM\t%s", mfp.name, pwm_row_ids[nucleotide]); */
/* 	  for (pwm_position = (*first_match_pwm_p).match.position + counter - width; */
/* 	       pwm_position < (*second_match_pwm_p).match.position + (*second_match_pwm_p).width - 1 - counter2 - width; pwm_position++) */
/* 	    printf("\t%.0f", mfp.fraction[nucleotide][pwm_position]); */
/* 	} */
/* 	printf("\n"); */

/*       } */
/*     } */
/*     exit(0); */
/*   } */

/* LOGO OPTION ENTIRE MAIN PROGRAM CODE */
  if (svg_only == 1 && argc >= 2 + firstnoncommandposition) {
    strcpy(searchstring, argv[1 + firstnoncommandposition]);
    qp.negative_values_allowed = negative_values_allowed;

    if (strstr(searchstring + strlen(searchstring) - 4, ".adm") == 0) {
      Load_pwm(&qp, searchstring, 1);
      strcpy(qp.name, searchstring);

      if (difference_logo == 1 && argc >= 3 + firstnoncommandposition) {
	offset = 0;
	if (argc >= 4 + firstnoncommandposition)
	  offset = atoi(argv[3 + firstnoncommandposition]);
	strcpy(tempstring, argv[2 + firstnoncommandposition]);
	Load_pwm(&mfp, tempstring, 1);
	Subtract_normalized_pwm(&qp, &mfp, offset);
	qp.negative_values_allowed = 1;
	strcat(qp.name, "_minus_");
	strcat(qp.name, tempstring);
      }

      np_p[0] = &qp;
      strcpy(tempstring, qp.name);
      strcpy(searchstring, qp.name);
      strcat(tempstring, ".svg");
      if (argc >= 3 + firstnoncommandposition && difference_logo == 0) {
	strcpy(tempstring, argv[2 + firstnoncommandposition]);
	strcpy(searchstring, argv[2 + firstnoncommandposition]);
      }
      Svg_logo(tempstring, 1, np_p, '\0', '\0', '\0', '\0', '\0', '\0', 0, 0);
    }
    /* else { */
    /*   printf("\nloads ADM"); */
    /*   adjacent_dinucleotide_model_init(&unflanked_adm, "unflanked_adm", Nlength); */
    /*   Load_ADM(&unflanked_adm, searchstring); */
    /*   strcpy(unflanked_adm.name, searchstring); */
    /*   strcpy(tempstring, unflanked_adm.name); */
    /*   strcpy(searchstring, unflanked_adm.name); */
    /*   strcat(tempstring, ".svg"); */
    /*   if (argc >= 3 + firstnoncommandposition) {   // JARKKO 28.4.2017 */
    /* 	strcpy(tempstring, argv[2 + firstnoncommandposition]); */
    /* 	strcpy(searchstring, argv[2 + firstnoncommandposition]); */
    /*   } */
    /*   Svg_riverlake_logo(tempstring, 0, 0, &unflanked_adm, 0.05, 0.1, 1000, 0.1); */
    /* } */

    char *system_command;
    system_command = malloc(2000);
    strcpy(system_command, "convert ");
    strcat(system_command, tempstring);
    strcat(system_command, " ");
    strcat(system_command, searchstring);
    strcat(system_command, ".png");
//printf("\n%s\n", system_command);
//system(system_command);
    exit(0);
  }
  /* SOMETHING WRONG, PRINT INSTRUCTIONS */
  else {
    printf
      ("Usage: \n"
       "./spacek --tool -options [background file name] [sample file name] [shortest kmer] [longest kmer] [minimum count or fold (with option -s) for printing] [PWM file (only with option -kmer)] \n"
       "\tCount number of occurrences of subsequences (kmers) of the indicated size. Indicates if count is higher than sequences within 1 Huddinge distance.\n"
       "\tWrites also a graphical summary to an svg file. \n"
       "\n"
       "  Other tools:\n"
       "  --f [background file name or - for none] [sample file name] [seed kmer sequence (consensus or IUPAC) or number for automatic (number specifies rank of local max start seed)] [minimum incidence] or\n"
       "  --pwm [background file name or - for none] [sample file name] [pwm file name] [score cut-off e.g. 0.1 or 10%%] [minimum incidence] \n"
       "\tFind occurrences of input sequence / pwm and print flanking preference, spacing and orientation information. \n"
       "\tGraphical summary of data is also written to an svg file\n"
       "\n"
       "  --logo [PWM file name] [output file name (opt)]\n"
       "\tGenerate SVG logo file from input PWM, also generates .png if convert is installed in path\n"
       "  --difflogo [PWM1 file name] [PWM2 file name] [offset (opt)] \n"
       "\tGenerate difference logo (PWM1-PWM2)\n"
       "  --dist [PWM1 file name] [PWM2 file name] ['gapped' or kmer1 length (opt)] [length of gapped kmer or kmer2 (opt)] \n"
       "\tCalculate uncentered correlation between maximum scores from PWM1 and PWM2 for all kmers. Default kmers are ungapped, default lengths are PWM widths (8 for gapped)\n"
       "  --pwmalign [PWM1 file name] [PWM2 file name] [PWM3 file name (optional)]\n"
       "\tAlign PWMs and generate a svg with logos aligned to PWM2, with boxes indicating similar base positions.\n"
       "\tOption -contacts shows also base contacts of PWM1 and PWM3 on PWM2.\n"
       "\n"
       "Options:\n"
       "\n"
       "Logo and sequence formatting:\n"
       "  -contacts\tShow hydrogen bond contacts (PWMs where contact information is included must be used)\n"
       "  -paths\tGenerate PWM logos using paths and not fonts\n"
       "  -neg\tAllow negative values in PWM (using paths)\n"
       "  -nocall\t do not include kmer-based prediction calls into logos\n"
       "  -core=[number]\tHilight this number of positions in the middle of the logo\n"
       "  -noname\tDo not include filename in the logo\n"
       "  -barcodelogo\tGenerate PWM barcode logos using colored rectangles\n"
       "  -heightscaledbars\tScale height of barcode logos based on nucleotide frequency\n"
       "  -maxheightscaledbars\tScale height of barcode logos based on maximum nucleotide frequency at position\n"
       "  -colorscaledbars\tScale color of barcode logos based on nucleotide frequency\n"
       "  -label\tInclude white IUPAC label to barcodelogo (cutoff = more than 50%% of max)\n"
       "  -rna\tInput sequence is from RNA, use uracil (U) instead of thymidine (T) in sequences and logos\n"
       "  -CpG\tHighlight CpGs in logos if they occur with lower (opaque) or higher (black outline) frequency in sample file than in control file \n"
       "\n"
       "Input sequences:\n"
       "  -14N\tSet sequence length to 14N (default is 20N)\n"
       "  -30N\tSet sequence length to 30N\n"
       "  -40N\tSet sequence length to 40N\n"
       "  -u\tUse only unique input sequences\n"
       "\n"
       "Kmer counting:\n"
       "  -nogaps\tCount only kmers without gaps\n"
       "  -allgaps\tCount kmers with gaps in any position (default is only middle 1 or 2 positions)\n"
       "  -longer≈[kmer_length _difference_cutoff for local max] (default = 0.4)]\n"
       "  -iupac≈[cutoff for making a base to an iupac] (default = 0.25 of maximum base)]\n"
       "\n"
       "Matching and PWM generation:\n"
       "  -m=[number]\tGenerate PWM using multinomial [number] distribution\n"
       "  -iterate\titerate input seed based on generated one hit pwm allowing up to 2 bp longer seed per round\n"
       "  -iterate-samesize\titerate seed, limit seed size to specified length\n"
       "  -iterate-fast\t iterate seed using free seed length\n"
       "  -lim=[position followed by strand (one or both must be given)]\t Show only hits at indicated strand and/or position (e.g. -lim=F, -lim=5 or -lim=4F)\n"
       "  -exc=[position][strand],[position][strand]‚..,[position][strand]\tExclude indicated positions and strands (e.g. -exc=4,3F,6R)\n"
       "  -both\tCount both instances of palindromic hits\n"
       "  -bothifnotequal\tCount both instances of palindromic hits if hit scores are not equal\n"
       "  -forwardonly\tCount only forward instance of palindromic hits (default counts hit from strand with better score, if score equal alternates between strands)\n"
       "  -reverseonly\tCount only reverse instance of palindromic hits\n"
       "  -forwardifequal\tCount only forward instance of palindromic hits if scores equal (default alternates between strands)\n"
       "  -reverseifequal\tCount only reverse instance of palindromic hits if scores equal\n"
       "  -e\tUse even background instead of background from file\n"
       "  -mono\tUse mononucleotide background instead of multinomial background\n"
       "  -fk=[position] count 4-mers at given position relative to match (e.g. -fk=-4 counts 4-mers that precede search sequence)\n"
       "  -mf=[multinomial] match filter: consider only sequences with one hit to indicated multinomial as true one hits\n"
       "  -kl=[kmer length] kmer length used to calculate lambda and automatic start seed (default 8)\n"
       "\n"
       "Text output:\n"
       "  -q\tPrint frequencies instead of counts\n"
       "  -n\tPrint nucleotide counts for full sequences\n"
       "  -s\tPrint values for different gap lengths on same line\n"
       "  -c\tPrint raw counts from both input files\n"
       "  -p\tPrint p-values (Winflat program needs to be installed in path)\n"
       "  -i\tPrint incidence of all input sequences\n"
       "  -kmer\tPrint kmer count table with scores for each kmer against the PWM indicated\n"
       "  -editdist[≈min_local_max_percent_cutoff (default = 10)] (for --f option)\tPrint table with max count kmer pairs, all local maxima, and all cloud counts for each edit distance\n"
       "  -dimer\tPrint counts for dimeric sequences only\n"
       "  -ic\tOutput information content for all spacings\n"
       "  -dinuc\tPrint dinucleotide data to output and svg (counting dinucleotides uses multinomial 2)\n"
       "  -x\tExtended output for debugging\n"
       "\n"
       "Svg output:\n"
       "  -o=[output file name] base name of svg output file(s)\n"
       "  -match=[kmer1],[kmer2]\tInclude spacing and orientation heatmap for these kmers in kmer summary svg\n"
       "  -hrows≈[number of rows for kmer svg summary heatmaps (default = 20)]\n"
       "  -xyplot\tGenerate scatterplot of kmer counts for the shortest kmer, with CpG-containing kmers indicated in red\n"
       "  -eoplot\tGenerate scatterplot of kmer counts observed and expected from the one hit PWM\n"
       "\n"
       " ");

    exit(0);
  }


}
/* END OF LOGO OPTION */

  
/* /\* /\\* DIST OPTION ENTIRE MAIN PROGRAM CODE *\\/ *\/ */
/* /\*   if (modeldistance == 1 && argc >= 3 + firstnoncommandposition) { *\/ */
/* /\*     short int gapped = 0; *\/ */
/* /\*     short int gappedkmerlength = 4; *\/ */
/* /\*     short int gap; *\/ */
/* /\*     long double onetwo; *\/ */
/* /\*     long double twoone; *\/ */
/* /\*     long double gapped_correlation; *\/ */
/* /\*     long double firstscore; *\/ */
/* /\*     long double secondscore; *\/ */
/* /\*     double **firstpwm; *\/ */
/* /\*     double **secondpwm; *\/ */
/* /\*     signed short int firstlength; *\/ */
/* /\*     signed short int secondlength; *\/ */
/* /\*     signed short int secondpwmwidth; *\/ */
/* /\*     signed short int firstpwmwidth; *\/ */
/* /\*     signed short int max_pwm_width; *\/ */
/* /\*     long long int kmervalue; *\/ */
/* /\*     long long int kmervalue2; *\/ */

/* /\*     strcpy(searchstring, argv[1 + firstnoncommandposition]); *\/ */
/* /\*     printf("\n\nNORMALIZED MATRIX 1:\t%s", searchstring); *\/ */
/* /\*     Load_pwm(&qp, searchstring, 0); *\/ */
/* /\*     firstlength = qp.width; *\/ */
/* /\*     Normalize_pwm(&qp); *\/ */
/* /\*     printf("\n\nLOG MATRIX 1:\t%s", searchstring); *\/ */
/* /\*     Log_ratio_pwm(&qp); *\/ */

/* /\*     strcpy(tempstring, argv[2 + firstnoncommandposition]); *\/ */
/* /\*     printf("\n\nNORMALIZED MATRIX 2:\t%s", tempstring); *\/ */
/* /\*     Load_pwm(&mfp, tempstring, 0); *\/ */
/* /\*     secondlength = mfp.width; *\/ */
/* /\*     Normalize_pwm(&mfp); *\/ */
/* /\*     printf("\n\nLOG MATRIX 2:\t%s", tempstring); *\/ */
/* /\*     Log_ratio_pwm(&mfp); *\/ */

/* /\*     if (firstlength > secondlength) *\/ */
/* /\*       too_long_kmer = firstlength; *\/ */
/* /\*     else *\/ */
/* /\*       too_long_kmer = secondlength; *\/ */
/* /\*     max_pwm_width = too_long_kmer; *\/ */

/* /\*     if (argc >= 4 + firstnoncommandposition) { *\/ */
/* /\*       if (strcmp(argv[3 + firstnoncommandposition], "gapped") == 0) *\/ */
/* /\* 	gapped = 1; *\/ */
/* /\*       else *\/ */
/* /\* 	firstlength = atoi(argv[3 + firstnoncommandposition]); *\/ */
/* /\*     } *\/ */
/* /\*     if (argc >= 5 + firstnoncommandposition) { *\/ */
/* /\*       secondlength = atoi(argv[4 + firstnoncommandposition]); *\/ */
/* /\*       secondlength = secondlength - secondlength % 2; *\/ */
/* /\*       gappedkmerlength = secondlength / 2; *\/ */
/* /\*     } *\/ */
/* /\*     if (gapped == 1) { *\/ */
/* /\*       firstlength = secondlength / 2; *\/ */
/* /\*       if (max_pwm_width < gappedkmerlength * 2) *\/ */
/* /\* 	gappedkmerlength = max_pwm_width / 2; *\/ */
/* /\*     } *\/ */

/* /\*     if (firstlength > secondlength) *\/ */
/* /\*       too_long_kmer = firstlength; *\/ */
/* /\*     else *\/ */
/* /\*       too_long_kmer = secondlength; *\/ */

/* /\*     Add_flanks(&qp, too_long_kmer - 2, -0.60206); *\/ */
/* /\*     Add_flanks(&mfp, too_long_kmer - 2, -0.60206); *\/ */

/* /\*     firstpwm = qp.fraction; *\/ */
/* /\*     firstpwmwidth = qp.width; *\/ */

/* /\*     secondpwm = mfp.fraction; *\/ */
/* /\*     secondpwmwidth = mfp.width; *\/ */

/* /\*     long long int last_kmer_value; *\/ */

/* /\*     sum_x_squared = 0; *\/ */
/* /\*     sum_y_squared = 0; *\/ */
/* /\*     sum_xy = 0; *\/ */
/* /\*     long int different = 0; *\/ */
/* /\*     long int same = 0; *\/ */

/* /\*     last_kmer_value = pow(4, firstlength); *\/ */
/* /\* #pragma omp parallel for schedule(static) reduction(+:sum_xy, sum_x_squared, sum_y_squared) default(shared) private(kmervalue,firstscore,secondscore) *\/ */
/* /\*     for (kmervalue = last_kmer_value; kmervalue > 0; kmervalue--) { *\/ */
/* /\*       firstscore = fastKmerscore(firstpwm, firstpwmwidth, kmervalue, firstlength); *\/ */
/* /\*       secondscore = fastKmerscore(secondpwm, secondpwmwidth, kmervalue, firstlength); *\/ */
/* /\*       sum_x_squared += pow(10, firstscore + firstscore); *\/ */
/* /\*       sum_y_squared += pow(10, secondscore + secondscore); *\/ */
/* /\*       sum_xy += pow(10, firstscore + secondscore); *\/ */
/* /\*       /\\*     if(firstscore != secondscore) different++; *\/ */
/* /\*          else same++; *\\/ *\/ */
/* /\*     } *\/ */
/* /\*     onetwo = sum_xy / (sqrt(sum_y_squared) * sqrt(sum_x_squared)); *\/ */

/* /\*     /\\* printf("\ndifferent %li, same %li\n", different, same); *\\/ *\/ */

/* /\*     if (firstlength != secondlength) { *\/ */
/* /\*       sum_x_squared = 0; *\/ */
/* /\*       sum_y_squared = 0; *\/ */
/* /\*       sum_xy = 0; *\/ */

/* /\*       last_kmer_value = pow(4, secondlength); *\/ */
/* /\* #pragma omp parallel for schedule(static) reduction(+:sum_xy, sum_x_squared, sum_y_squared) default(shared) private(kmervalue,firstscore,secondscore) *\/ */
/* /\*       for (kmervalue = 0; kmervalue < last_kmer_value; kmervalue++) { *\/ */
/* /\* 	firstscore = fastKmerscore(firstpwm, firstpwmwidth, kmervalue, secondlength); *\/ */
/* /\* 	secondscore = fastKmerscore(secondpwm, secondpwmwidth, kmervalue, secondlength); *\/ */
/* /\* 	sum_x_squared += pow(10, firstscore + firstscore); *\/ */
/* /\* 	sum_y_squared += pow(10, secondscore + secondscore); *\/ */
/* /\* 	sum_xy += pow(10, firstscore + secondscore); *\/ */
/* /\*       } *\/ */
/* /\*       twoone = sum_xy / (sqrt(sum_y_squared) * sqrt(sum_x_squared)); *\/ */
/* /\*     } else *\/ */
/* /\*       twoone = onetwo; *\/ */

/* /\*     /\\* GAPPED kmers *\\/ *\/ */
/* /\*     if (gapped == 1) { *\/ */

/* /\*       sum_x_squared = 0; *\/ */
/* /\*       sum_y_squared = 0; *\/ */
/* /\*       sum_xy = 0; *\/ */

/* /\*       last_kmer_value = pow(4, gappedkmerlength); *\/ */

/* /\* #pragma omp parallel for schedule(static) reduction(+:sum_xy, sum_x_squared, sum_y_squared) default(shared) private(kmervalue,kmervalue2, gap, firstscore,secondscore) *\/ */
/* /\*       for (gap = max_pwm_width - 2 * gappedkmerlength; gap >= 0; gap--) { *\/ */
/* /\* 	for (kmervalue = last_kmer_value; kmervalue > 0; kmervalue--) { *\/ */
/* /\* 	  for (kmervalue2 = last_kmer_value; kmervalue2 > 0; kmervalue2--) { *\/ */

/* /\* 	    firstscore = fastgappedKmerscore(firstpwm, firstpwmwidth, kmervalue, kmervalue2, gappedkmerlength, gap); *\/ */
/* /\* 	    secondscore = fastgappedKmerscore(secondpwm, secondpwmwidth, kmervalue, kmervalue2, gappedkmerlength, gap); *\/ */
/* /\* 	    sum_x_squared += pow(10, firstscore + firstscore); *\/ */
/* /\* 	    sum_y_squared += pow(10, secondscore + secondscore); *\/ */
/* /\* 	    sum_xy += pow(10, firstscore + secondscore); *\/ */
/* /\* 	  } *\/ */
/* /\* 	} *\/ */
/* /\*       } *\/ */
/* /\*       gapped_correlation = sum_xy / (sqrt(sum_y_squared) * sqrt(sum_x_squared)); *\/ */

/* /\*       printf("\n\nSimilarity (cosine angle correlation) \n%s to %s (%i-mer)\t%Lg\n%s to %s (%i-mer)\t%Lg\n%s to %s (gapped %i-mer)\t%Lg\n", *\/ */
/* /\* 	     searchstring, tempstring, firstlength, onetwo, tempstring, searchstring, secondlength, twoone, searchstring, tempstring, *\/ */
/* /\* 	     gappedkmerlength * 2, gapped_correlation); *\/ */
/* /\*     } else { *\/ */
/* /\*       printf("\n\nSimilarity (cosine angle correlation) \n%s to %s (%i-mer)\t%Lg\n%s to %s (%i-mer)\t%Lg\n", searchstring, tempstring, firstlength, *\/ */
/* /\* 	     onetwo, tempstring, searchstring, secondlength, twoone); *\/ */
/* /\*     } *\/ */

/* /\*     t1 = time('\0'); *\/ */
/* /\*     printf("\nTime: %ld seconds\n", (long)(t1 - t0)); *\/ */
/* /\*     exit(0); *\/ */
/* /\*   } *\/ */

/* /\* /\\* END OF DIST OPTION *\\/ *\/ */

/* /\* /\\* CHECKS IF WINFLAT IS INSTALLED *\\/ *\/ */
/* /\*   if (print_p_values == 1) *\/ */
/* /\*     if (Winflat(10, 10, 100, 100) == 2) { *\/ */
/* /\*       printf("\nERROR: Winflat not installed in path\n"); *\/ */
/* /\*       argc = 1; *\/ */
/* /\*       fflush(stdout); *\/ */
/* /\*     } *\/ */

/* /\* /\\* PARSES ARGUMENTS *\\/ *\/ */
/* /\*   if ((argc >= 5 + firstnoncommandposition || (flank == 1 && argc >= 3 + firstnoncommandposition)) && strcmp("HELP", argv[1]) != 0 *\/ */
/* /\*       && strcmp("help", argv[1]) != 0 && strcmp("-h", argv[1]) != 0 && strcmp("-H", argv[1]) != 0) { *\/ */
/* /\* /\\* FLANK TOOL *\\/ *\/ */
/* /\*     if (flank == 1) { *\/ */
/* /\*       strcpy(file_name[0], argv[1 + firstnoncommandposition]); *\/ */
/* /\*       if (strcmp(file_name[0], "-") == 0) { *\/ */
/* /\* 	even_background = 1; *\/ */
/* /\* 	file_number++; *\/ */
/* /\*       } *\/ */
/* /\*       strcpy(file_name[1], argv[2 + firstnoncommandposition]); *\/ */
/* /\*       strcpy(searchstring, argv[3 + firstnoncommandposition]); *\/ */
/* /\* /\\* CHECKS IF SEARCHSTRING IS A NUCLEOTIDE SEQUENCE, IF NOT, ASSUMES IT IS A FILENAME FOR PWM *\\/ *\/ */
/* /\*       for (counter = 0; counter < strlen(searchstring); counter++) *\/ */
/* /\* 	if (searchstring[counter] != 'A' & searchstring[counter] != 'C' & searchstring[counter] != 'G' & searchstring[counter] != 'T') *\/ */
/* /\* 	  flank_with_pwm = 1; *\/ */
/* /\*       if (flank_with_pwm == 1) { *\/ */
/* /\* 	iupac_query = 1; *\/ */
/* /\* 	flank_with_pwm = 1; *\/ */
/* /\* 	for (counter = 0; counter < strlen(searchstring); counter++) { *\/ */
/* /\* 	  for (counter2 = 0; (counter2 < strlen(nucleotide_iupac)) & (searchstring[counter] != nucleotide_iupac[counter2]); counter2++) ; *\/ */
/* /\* 	  if (counter2 == strlen(nucleotide_iupac)) { *\/ */
/* /\* 	    iupac_query = 0; *\/ */
/* /\* 	    flank_with_pwm = 1; *\/ */
/* /\* 	    break; *\/ */
/* /\* 	  } *\/ */
/* /\* 	} *\/ */
/* /\*       } *\/ */

/* /\*       if (argc >= 5 + firstnoncommandposition) { *\/ */
/* /\* 	strcpy(valuestring, argv[4 + firstnoncommandposition]); *\/ */
/* /\* 	if (valuestring[strlen(valuestring) - 1] == '%') { *\/ */
/* /\* 	  valuestring[strlen(valuestring) - 1] = '\0'; *\/ */
/* /\* 	  minimum_kmer_count = atof(valuestring) / 100; *\/ */
/* /\* 	} else *\/ */
/* /\* 	  minimum_kmer_count = atof(valuestring); *\/ */
/* /\*       } *\/ */
/* /\*       if (cut_off != 0) *\/ */
/* /\* 	cut_off = log10(minimum_kmer_count); *\/ */
/* /\*       else *\/ */
/* /\* 	cut_off = log10(pseudocount); *\/ */
/* /\*       if (argc >= 6 + firstnoncommandposition && flank_with_pwm == 1) *\/ */
/* /\* 	minimum_kmer_count = atof(argv[5 + firstnoncommandposition]); *\/ */
/* /\*       pwm_name = searchstring; *\/ */
/* /\*     } *\/ */
/* /\* /\\* NMER COUNT TOOL *\\/ *\/ */
/* /\*     else { *\/ */
/* /\*       strcpy(file_name[0], argv[1 + firstnoncommandposition]); *\/ */
/* /\*       strcpy(file_name[1], argv[2 + firstnoncommandposition]); *\/ */
/* /\*       shortest_kmer = atoi(argv[3 + firstnoncommandposition]); *\/ */
/* /\*       too_long_kmer = atoi(argv[4 + firstnoncommandposition]) + 1; *\/ */
/* /\*       minimum_kmer_count = atof(argv[5 + firstnoncommandposition]); *\/ */
/* /\*       if (kmer_table == 1 && argc >= 6 + firstnoncommandposition) { *\/ */
/* /\* 	strcpy(searchstring, argv[6 + firstnoncommandposition]); *\/ */

/* /\* /\\* CHECKS IF SEARCHSTRING IS A NUCLEOTIDE SEQUENCE, IF NOT, ASSUMES IT IS A FILENAME FOR PWM *\\/ *\/ */
/* /\* 	for (counter = 0; counter < strlen(searchstring); counter++) *\/ */
/* /\* 	  if (searchstring[counter] != 'A' & searchstring[counter] != 'C' & searchstring[counter] != 'G' & searchstring[counter] != 'T') *\/ */
/* /\* 	    load_pwm = 1; *\/ */

/* /\* 	if (load_pwm == 1) { *\/ */
/* /\* 	  iupac_query = 1; *\/ */
/* /\* 	  load_pwm = 0; *\/ */
/* /\* 	  for (counter = 0; counter < strlen(searchstring); counter++) { *\/ */
/* /\* 	    for (counter2 = 0; (counter2 < strlen(nucleotide_iupac)) & (searchstring[counter] != nucleotide_iupac[counter2]); counter2++) ; *\/ */
/* /\* 	    if (counter2 == strlen(nucleotide_iupac)) { *\/ */
/* /\* 	      iupac_query = 0; *\/ */
/* /\* 	      load_pwm = 1; *\/ */
/* /\* 	      break; *\/ */
/* /\* 	    } *\/ */
/* /\* 	  } *\/ */
/* /\* 	} *\/ */

/* /\* 	if (load_pwm == 1) { *\/ */
/* /\* 	  if (strstr(searchstring + strlen(searchstring) - 4, ".adm") == 0) { *\/ */
/* /\* 	    Load_pwm(&qp, searchstring, 1); *\/ */
/* /\* 	    Log_ratio_pwm(&qp); *\/ */
/* /\* 	  } else { *\/ */
/* /\* 	    load_adm = 1; *\/ */
/* /\* 	    printf("\nloads ADM"); *\/ */
/* /\* 	    adjacent_dinucleotide_model_init(&unflanked_adm, "unflanked_adm", Nlength); *\/ */
/* /\* 	    adjacent_dinucleotide_model_init(&flanked_adm, "unflanked_adm", Nlength + 2 * (too_long_kmer - 2)); *\/ */
/* /\* 	    Load_ADM(&unflanked_adm, searchstring); *\/ */
/* /\* 	    Log_ratio_ADM(&unflanked_adm); *\/ */
/* /\* 	    PWM_from_ADM(&unflanked_adm, &qp); *\/ */
/* /\* 	    Add_flanks_to_ADM(&flanked_adm, &unflanked_adm, too_long_kmer - 2, -0.60206); *\/ */
/* /\* 	    /\\* short int *jj; *\/ */
/* /\* 	       Kmerscore_ADM(&flanked_adm, 5762,8, jj); *\/ */
/* /\* 	       Kmerscore_ADM(&flanked_adm, 24763,8, jj); *\/ */
/* /\* 	       Kmerscore_ADM(&flanked_adm, 24762,8, jj); *\/ */
/* /\* 	       exit(0); *\\/ *\/ */
/* /\* 	  } *\/ */
/* /\* 	} else { *\/ */
/* /\* 	  Multinomial_pwm(&qp, searchstring); *\/ */
/* /\* 	  if (match_filter != 0) *\/ */
/* /\* 	    Multinomial_pwm(&mfp, searchstring); *\/ */
/* /\* 	} *\/ */
/* /\* 	loaded_pwm_width = qp.width; *\/ */
/* /\* 	Add_flanks(&qp, too_long_kmer - 2, -0.60206); *\/ */
/* /\*       } *\/ */
/* /\*     } *\/ */
/* /\*   } else *\/ */
/* /\* /\\* SOMETHING WRONG, PRINT INSTRUCTIONS *\\/ *\/ */
/* /\*   { *\/ */
/* /\*     printf *\/ */
/* /\*       ("Usage: \n" *\/ */
/* /\*        "./spacek --tool -options [background file name] [sample file name] [shortest kmer] [longest kmer] [minimum count or fold (with option -s) for printing] [PWM file (only with option -kmer)] \n" *\/ */
/* /\*        "\tCount number of occurrences of subsequences (kmers) of the indicated size. Indicates if count is higher than sequences within 1 Huddinge distance.\n" *\/ */
/* /\*        "\tWrites also a graphical summary to an svg file. \n" *\/ */
/* /\*        "\n" *\/ */
/* /\*        "  Other tools:\n" *\/ */
/* /\*        "  --f [background file name or - for none] [sample file name] [seed kmer sequence (consensus or IUPAC) or number for automatic (number specifies rank of local max start seed)] [minimum incidence] or\n" *\/ */
/* /\*        "  --pwm [background file name or - for none] [sample file name] [pwm file name] [score cut-off e.g. 0.1 or 10%%] [minimum incidence] \n" *\/ */
/* /\*        "\tFind occurrences of input sequence / pwm and print flanking preference, spacing and orientation information. \n" *\/ */
/* /\*        "\tGraphical summary of data is also written to an svg file\n" *\/ */
/* /\*        "\n" *\/ */
/* /\*        "  --logo [PWM file name] [output file name (opt)]\n" *\/ */
/* /\*        "\tGenerate SVG logo file from input PWM, also generates .png if convert is installed in path\n" *\/ */
/* /\*        "  --difflogo [PWM1 file name] [PWM2 file name] [offset (opt)] \n" *\/ */
/* /\*        "\tGenerate difference logo (PWM1-PWM2)\n" *\/ */
/* /\*        "  --dist [PWM1 file name] [PWM2 file name] ['gapped' or kmer1 length (opt)] [length of gapped kmer or kmer2 (opt)] \n" *\/ */
/* /\*        "\tCalculate uncentered correlation between maximum scores from PWM1 and PWM2 for all kmers. Default kmers are ungapped, default lengths are PWM widths (8 for gapped)\n" *\/ */
/* /\*        "  --pwmalign [PWM1 file name] [PWM2 file name] [PWM3 file name (optional)]\n" *\/ */
/* /\*        "\tAlign PWMs and generate a svg with logos aligned to PWM2, with boxes indicating similar base positions.\n" *\/ */
/* /\*        "\tOption -contacts shows also base contacts of PWM1 and PWM3 on PWM2.\n" *\/ */
/* /\*        "\n" *\/ */
/* /\*        "Options:\n" *\/ */
/* /\*        "\n" *\/ */
/* /\*        "Logo and sequence formatting:\n" *\/ */
/* /\*        "  -contacts\tShow hydrogen bond contacts (PWMs where contact information is included must be used)\n" *\/ */
/* /\*        "  -paths\tGenerate PWM logos using paths and not fonts\n" *\/ */
/* /\*        "  -neg\tAllow negative values in PWM (using paths)\n" *\/ */
/* /\*        "  -nocall\t do not include kmer-based prediction calls into logos\n" *\/ */
/* /\*        "  -core=[number]\tHilight this number of positions in the middle of the logo\n" *\/ */
/* /\*        "  -noname\tDo not include filename in the logo\n" *\/ */
/* /\*        "  -barcodelogo\tGenerate PWM barcode logos using colored rectangles\n" *\/ */
/* /\*        "  -heightscaledbars\tScale height of barcode logos based on nucleotide frequency\n" *\/ */
/* /\*        "  -maxheightscaledbars\tScale height of barcode logos based on maximum nucleotide frequency at position\n" *\/ */
/* /\*        "  -colorscaledbars\tScale color of barcode logos based on nucleotide frequency\n" *\/ */
/* /\*        "  -label\tInclude white IUPAC label to barcodelogo (cutoff = more than 50%% of max)\n" *\/ */
/* /\*        "  -rna\tInput sequence is from RNA, use uracil (U) instead of thymidine (T) in sequences and logos\n" *\/ */
/* /\*        "  -CpG\tHighlight CpGs in logos if they occur with lower (opaque) or higher (black outline) frequency in sample file than in control file \n" *\/ */
/* /\*        "\n" *\/ */
/* /\*        "Input sequences:\n" *\/ */
/* /\*        "  -14N\tSet sequence length to 14N (default is 20N)\n" *\/ */
/* /\*        "  -30N\tSet sequence length to 30N\n" *\/ */
/* /\*        "  -40N\tSet sequence length to 40N\n" *\/ */
/* /\*        "  -u\tUse only unique input sequences\n" *\/ */
/* /\*        "\n" *\/ */
/* /\*        "Kmer counting:\n" *\/ */
/* /\*        "  -nogaps\tCount only kmers without gaps\n" *\/ */
/* /\*        "  -allgaps\tCount kmers with gaps in any position (default is only middle 1 or 2 positions)\n" *\/ */
/* /\*        "  -longer≈[kmer_length _difference_cutoff for local max] (default = 0.4)]\n" *\/ */
/* /\*        "  -iupac≈[cutoff for making a base to an iupac] (default = 0.25 of maximum base)]\n" *\/ */
/* /\*        "\n" *\/ */
/* /\*        "Matching and PWM generation:\n" *\/ */
/* /\*        "  -m=[number]\tGenerate PWM using multinomial [number] distribution\n" *\/ */
/* /\*        "  -iterate\titerate input seed based on generated one hit pwm allowing up to 2 bp longer seed per round\n" *\/ */
/* /\*        "  -iterate-samesize\titerate seed, limit seed size to specified length\n" *\/ */
/* /\*        "  -iterate-fast\t iterate seed using free seed length\n" *\/ */
/* /\*        "  -lim=[position followed by strand (one or both must be given)]\t Show only hits at indicated strand and/or position (e.g. -lim=F, -lim=5 or -lim=4F)\n" *\/ */
/* /\*        "  -exc=[position][strand],[position][strand]‚..,[position][strand]\tExclude indicated positions and strands (e.g. -exc=4,3F,6R)\n" *\/ */
/* /\*        "  -both\tCount both instances of palindromic hits\n" *\/ */
/* /\*        "  -bothifnotequal\tCount both instances of palindromic hits if hit scores are not equal\n" *\/ */
/* /\*        "  -forwardonly\tCount only forward instance of palindromic hits (default counts hit from strand with better score, if score equal alternates between strands)\n" *\/ */
/* /\*        "  -reverseonly\tCount only reverse instance of palindromic hits\n" *\/ */
/* /\*        "  -forwardifequal\tCount only forward instance of palindromic hits if scores equal (default alternates between strands)\n" *\/ */
/* /\*        "  -reverseifequal\tCount only reverse instance of palindromic hits if scores equal\n" *\/ */
/* /\*        "  -e\tUse even background instead of background from file\n" *\/ */
/* /\*        "  -mono\tUse mononucleotide background instead of multinomial background\n" *\/ */
/* /\*        "  -fk=[position] count 4-mers at given position relative to match (e.g. -fk=-4 counts 4-mers that precede search sequence)\n" *\/ */
/* /\*        "  -mf=[multinomial] match filter: consider only sequences with one hit to indicated multinomial as true one hits\n" *\/ */
/* /\*        "  -kl=[kmer length] kmer length used to calculate lambda and automatic start seed (default 8)\n" *\/ */
/* /\*        "\n" *\/ */
/* /\*        "Text output:\n" *\/ */
/* /\*        "  -q\tPrint frequencies instead of counts\n" *\/ */
/* /\*        "  -n\tPrint nucleotide counts for full sequences\n" *\/ */
/* /\*        "  -s\tPrint values for different gap lengths on same line\n" *\/ */
/* /\*        "  -c\tPrint raw counts from both input files\n" *\/ */
/* /\*        "  -p\tPrint p-values (Winflat program needs to be installed in path)\n" *\/ */
/* /\*        "  -i\tPrint incidence of all input sequences\n" *\/ */
/* /\*        "  -kmer\tPrint kmer count table with scores for each kmer against the PWM indicated\n" *\/ */
/* /\*        "  -editdist[≈min_local_max_percent_cutoff (default = 10)] (for --f option)\tPrint table with max count kmer pairs, all local maxima, and all cloud counts for each edit distance\n" *\/ */
/* /\*        "  -dimer\tPrint counts for dimeric sequences only\n" *\/ */
/* /\*        "  -ic\tOutput information content for all spacings\n" *\/ */
/* /\*        "  -dinuc\tPrint dinucleotide data to output and svg (counting dinucleotides uses multinomial 2)\n" *\/ */
/* /\*        "  -x\tExtended output for debugging\n" *\/ */
/* /\*        "\n" *\/ */
/* /\*        "Svg output:\n" *\/ */
/* /\*        "  -o=[output file name] base name of svg output file(s)\n" *\/ */
/* /\*        "  -match=[kmer1],[kmer2]\tInclude spacing and orientation heatmap for these kmers in kmer summary svg\n" *\/ */
/* /\*        "  -hrows≈[number of rows for kmer svg summary heatmaps (default = 20)]\n" *\/ */
/* /\*        "  -xyplot\tGenerate scatterplot of kmer counts for the shortest kmer, with CpG-containing kmers indicated in red\n" *\/ */
/* /\*        "  -eoplot\tGenerate scatterplot of kmer counts observed and expected from the one hit PWM\n" *\/ */
/* /\*        "\n" *\/ */
/* /\*        " "); *\/ */

/* /\*     exit(0); *\/ */
/* /\*   } *\/ */

/* /\*   /\\* struct normalized_pwm rc; *\\/ *\/ */

/* /\*   /\\* INITIALIZES ALIGN SCORES *\\/ *\/ */
/* /\*   struct alignscore *all_hits_align_scores; *\/ */
/* /\*   all_hits_align_scores = malloc(sizeof(struct alignscore) * 3 + 5); *\/ */
/* /\*   alignscore_structure_init(&all_hits_align_scores[0]); *\/ */
/* /\*   alignscore_structure_init(&all_hits_align_scores[1]); *\/ */
/* /\*   struct alignscore current_align_score; *\/ */
/* /\*   alignscore_structure_init(&current_align_score); *\/ */

/* /\*   /\\* INITIALIZES DINUCLEOTIDE PROPERTIES STRUCTURES IF DINUCLEOTIDE PROPERTIES ARE ANALYZED *\\/ *\/ */

/* /\*   struct dinucleotide_properties_matrix *one_hit_di; *\/ */
/* /\*   struct base_dependency_matrix *dep; *\/ */
/* /\*   struct base_dependency_matrix expected_dinucleotides; *\/ */

/* /\*   if (dinucleotide_properties == 1) { *\/ */
/* /\*     dinucleotide_properties_init(&di); *\/ */
/* /\*     one_hit_di = malloc(sizeof(struct dinucleotide_properties_matrix) * (number_of_files + 1) + 5); *\/ */
/* /\*     dep = malloc(sizeof(struct base_dependency_matrix) * (number_of_files + 1) + 5); *\/ */
/* /\*     dinucleotide_properties_matrix_init(&one_hit_di[0], &di, "dinucleotide1", Nlength * 2, 0, query_sequence_length); *\/ */
/* /\*     dinucleotide_properties_matrix_init(&one_hit_di[1], &di, "dinucleotide2", Nlength * 2, 0, query_sequence_length); *\/ */
/* /\*     base_dependency_matrix_init(&dep[0], "Background_dinucleotide_dependencies", Nlength * 2); *\/ */
/* /\*     base_dependency_matrix_init(&dep[1], "Dinucleotide_dependencies", Nlength * 2); *\/ */
/* /\*     base_dependency_matrix_init(&expected_dinucleotides, "Expected_dinucleotides", Nlength * 2); *\/ */
/* /\*   } *\/ */

/* /\*   /\\* ALLOCATES MEMORY AND CLEARS INDEX IF ONLY UNIQUE SEQUENCES USED *\\/ *\/ */
/* /\*   if (remove_non_unique == 1) { *\/ */
/* /\*     sorted_list[0] = malloc(max_number_of_sequences * sizeof(struct sequence_incidence_table) + 10); *\/ */
/* /\*     sorted_list[1] = malloc(max_number_of_sequences * sizeof(struct sequence_incidence_table) + 10); *\/ */
/* /\*   } *\/ */

/* /\*   /\\* GENERATES mask_ULLS FOR EXTRACTION OF EACH KMER STRING *\\/ *\/ */
/* /\*   /\\* GENERATES KMER mask_ULLS *\\/ *\/ */
/* /\*   for (current_kmer_length = 2; current_kmer_length <= Nlength; current_kmer_length++) { *\/ */
/* /\*     for (start_position = 0; start_position < Nlength - current_kmer_length; start_position++) { *\/ */
/* /\*       for (mask_ULL[current_kmer_length][start_position] = mask_ULL[1][start_position], position = start_position + 1; *\/ */
/* /\* 	   position < current_kmer_length + start_position; position++) { *\/ */
/* /\* 	mask_ULL[current_kmer_length][start_position] += mask_ULL[1][position]; *\/ */
/* /\* 	/\\* printf("\n%i\t%i\t%i", current_kmer_length, start_position, mask_ULL[current_kmer_length][start_position]); *\\/ } *\/ */
/* /\*     } *\/ */
/* /\*   } *\/ */

/* /\*   /\\* GENERATES HIGH AND LOW mask_ULLS FOR DELETIONS *\\/ *\/ */
/* /\*   for (lowmask_ULL[0] = mask_ULL[1][0], position = 1; position < Nlength - 2; position++) *\/ */
/* /\*     lowmask_ULL[position] = lowmask_ULL[position - 1] + mask_ULL[1][position]; *\/ */
/* /\*   for (highmask_ULL[Nlength - 2] = mask_ULL[1][Nlength - 2], position = Nlength - 3; position > 0; position--) *\/ */
/* /\*     highmask_ULL[position] = highmask_ULL[position + 1] + mask_ULL[1][position]; *\/ */

/* /\*   if (flank_with_pwm == 1 || (multinomial != 0 && kmer_table == 0)) { *\/ */
/* /\*     too_long_kmer = shortest_kmer + 1; *\/ */
/* /\*     if (shortest_kmer < 13) *\/ */
/* /\*       kmer_count = 1; *\/ */
/* /\*   } *\/ */

/* /\*   /\\* INITIALIZES RESULT TABLE *\\/ *\/ */
/* /\*   long int *****results; *\/ */
/* /\*   long long int variable_length; *\/ */
/* /\*   results = malloc(sizeof(long int *) * 3 + 5); *\/ */

/* /\*   /\\* long int *results[number_of_files][too_long_kmer+1][too_long_kmer+2][Nlength-shortest_kmer+3]; *\\/ *\/ */
/* /\*   long int *current_value; *\/ */
/* /\*   short int kmer_length_size; *\/ */
/* /\*   short int gap_position_size; *\/ */
/* /\*   /\\* printf("\n\nsizeof int %i\t%i", sizeof(current_value), sizeof(current_value[0])); *\\/ *\/ */
/* /\*   /\\* clears results *\\/ *\/ */
/* /\*   for (counter = 0; counter < number_of_files; counter++) { *\/ */
/* /\*     results[counter] = malloc(sizeof(long int *) * too_long_kmer + 5); *\/ */
/* /\*     for (current_kmer_length = shortest_kmer; current_kmer_length < too_long_kmer; current_kmer_length++) { *\/ */
/* /\*       kmer_length_size = ((current_kmer_length - 2) * (count_also_spaced_kmers != 0)) + 2; *\/ */
/* /\*       results[counter][current_kmer_length] = malloc(sizeof(long int *) * kmer_length_size + 5); *\/ */
/* /\*       for (current_gap_position = 0; current_gap_position < kmer_length_size; current_gap_position++) { *\/ */
/* /\* 	/\\* DO NOT ALLOCATE MEMORY FOR ALL GAPS IF ONLY CENTER GAPS COUNTED *\\/ *\/ */
/* /\* 	if (count_also_spaced_kmers == 0 && current_gap_position > 1) *\/ */
/* /\* 	  continue; *\/ */
/* /\* 	if (count_also_spaced_kmers == 1 *\/ */
/* /\* 	    && (current_gap_position != current_kmer_length / 2 && current_gap_position != current_kmer_length / 2 + current_kmer_length % 2 *\/ */
/* /\* 		&& current_gap_position > 1)) *\/ */
/* /\* 	  continue; *\/ */

/* /\* 	gap_position_size = ((Nlength - current_kmer_length + 1) * (count_also_spaced_kmers != 0)) + 1; *\/ */
/* /\* 	results[counter][current_kmer_length][current_gap_position] = malloc(sizeof(long int *) * gap_position_size + 5); *\/ */
/* /\* 	for (current_gap_length = 0; current_gap_length < gap_position_size; current_gap_length++) { *\/ */
/* /\* 	  /\\* DO NOT ALLOCATE MEMORY FOR ALL GAPS IF ONLY CENTER GAPS COUNTED *\\/ *\/ */
/* /\* 	  if (current_gap_position > 1 && current_gap_length == 0) *\/ */
/* /\* 	    continue; *\/ */
/* /\* 	  if (current_gap_position == 0 && current_gap_length > 1) *\/ */
/* /\* 	    continue; *\/ */
/* /\* 	  if (count_also_spaced_kmers != 2 && current_gap_position == 1 && current_gap_length == 1 && current_kmer_length > 3) *\/ */
/* /\* 	    break; *\/ */

/* /\* 	  variable_length = pow(4, current_kmer_length); *\/ */
/* /\* 	  results[counter][current_kmer_length][current_gap_position][current_gap_length] = malloc(variable_length * sizeof(long int) + 5); *\/ */
/* /\* 	  // memset (results[counter][current_kmer_length][current_gap_position][current_gap_length], 0, pow(4, current_kmer_length) * sizeof(long int *) + 5); *\/ */
/* /\* 	  current_value = results[counter][current_kmer_length][current_gap_position][current_gap_length]; *\/ */
/* /\* 	  for (current_kmer = 0; current_kmer < variable_length; current_kmer++) { *\/ */
/* /\* 	    current_value[current_kmer] = 0; *\/ */
/* /\* 	  } *\/ */
/* /\* 	} *\/ */
/* /\*       } *\/ */
/* /\*     } *\/ */
/* /\*   } *\/ */

/* /\*   short int round; *\/ */
/* /\*   short int max_rounds = 20; *\/ */
/* /\*   char **seed_list; *\/ */
/* /\*   seed_list = malloc(sizeof(char *) * 100 + 5); *\/ */

/* /\*   /\\* SEED ITERATION MAIN LOOP *\\/ *\/ */
/* /\*   if (iterate_seed == 1) *\/ */
/* /\*     remember_iterate_seed = 1; *\/ */
/* /\*   for (round = 1; round < max_rounds + 1; round++) { *\/ */

/* /\*     if (flank == 1) { *\/ */
/* /\*       for (counter = 0; counter < strlen(searchstring); counter++) *\/ */
/* /\* 	if (searchstring[counter] != '0' & searchstring[counter] != '1' & searchstring[counter] != '2' & searchstring[counter] != *\/ */
/* /\* 	    '3' & searchstring[counter] != '4' & searchstring[counter] != '5' & searchstring[counter] != '6' & searchstring[counter] != *\/ */
/* /\* 	    '7' & searchstring[counter] != '8' & searchstring[counter] != '9') { *\/ */
/* /\* 	  auto_seed = 0; *\/ */
/* /\* 	} *\/ */
/* /\*     } else *\/ */
/* /\*       auto_seed = 0; *\/ */

/* /\*     if (auto_seed == 1) { *\/ */
/* /\*       printf("\nAUTOSEED"); *\/ */
/* /\*       seed_from_local_max_number = atoi(searchstring); *\/ */
/* /\*       iterate_seed = 1; *\/ */
/* /\*       iupac_query = 1; *\/ */
/* /\*       remember_iterate_seed = 1; *\/ */
/* /\*       load_pwm = 0; *\/ */
/* /\*       strcpy(searchstring, "CGSATATA"); *\/ */
/* /\*     } *\/ */

/* /\*     /\\* CALCULATES MULTINOMIAL PWM IF MULTINOMIAL INDICATED OR LOADS PWM IF FLANK WITH PWM = 1 *\\/ *\/ */
/* /\*     if (flank_with_pwm == 1 || (multinomial != 0 && kmer_table == 0)) { *\/ */
/* /\*       /\\* shortest_kmer = 8; *\\/ *\/ */
/* /\*       //too_long_kmer = shortest_kmer+1; *\/ */
/* /\*       //if (query_sequence_length < 12) kmer_count = 1; *\/ */
/* /\*       no_insertions = 1; *\/ */
/* /\*       flank_with_pwm = 1; *\/ */
/* /\*       /\\* printf("\nLoads Pwm: %s", pwm_name); *\\/ *\/ */
/* /\*       if (multinomial == 0 && iupac_query == 0) { *\/ */
/* /\* 	Load_pwm(&qp, pwm_name, 1); *\/ */
/* /\*       } else { *\/ */
/* /\* 	if (multinomial < 2 && dinucleotide_properties == 1) *\/ */
/* /\* 	  multinomial_2_dinuc = 1; *\/ */
/* /\* 	if (iupac_query == 1) { *\/ */
/* /\* 	  Iupac_to_pwm(&qp, searchstring);	/\\* CALCULATES PWM FOR IUPAC *\\/ *\/ */
/* /\* 	  if (match_filter != 0) *\/ */
/* /\* 	    Iupac_to_pwm(&mfp, searchstring);	/\\* CALCULATES PWM FOR MATCH FILTER *\\/ *\/ */
/* /\* 	  if (multinomial_2_dinuc == 1) *\/ */
/* /\* 	    Iupac_to_pwm(&mmp, searchstring);	/\\* CALCULATES PWM FOR DINUC MULTINOMIAL 2 FILTER *\\/ *\/ */
/* /\* 	} else { *\/ */
/* /\* 	  Multinomial_pwm(&qp, searchstring);	/\\* CALCULATES PWM FOR MULTINOMIAL *\\/ *\/ */
/* /\* 	  if (match_filter != 0) *\/ */
/* /\* 	    Multinomial_pwm(&mfp, searchstring);	/\\* CALCULATES PWM FOR MATCH FILTER *\\/ *\/ */
/* /\* 	  if (multinomial_2_dinuc == 1) *\/ */
/* /\* 	    Multinomial_pwm(&mmp, searchstring);	/\\* CALCULATES PWM FOR DINUC MULTINOMIAL 2 FILTER *\\/ *\/ */
/* /\* 	} *\/ */
/* /\*       } *\/ */
/* /\*       query_sequence_length = qp.width; *\/ */
/* /\*       /\\* ADDS PWMs TO LOGO LIST *\\/ *\/ */
/* /\*       Copy_normalized_pwm(&np[current_logo], &qp); *\/ */
/* /\*       strcpy(np[current_logo].name, pwm_name); *\/ */

/* /\*       if (multinomial == 0 && iupac_query == 0) *\/ */
/* /\* 	Log_fold_affinity_pwm(&qp); *\/ */
/* /\*       else { *\/ */
/* /\* 	cut_off = 0 - multinomial - 0.0001; *\/ */
/* /\* 	multi_2_cut_off = 0 - 2 - 0.0001; *\/ */
/* /\* 	filter_cut_off = 0 - match_filter - 0.0001; *\/ */
/* /\* 	sprintf(tempstring, " - multinomial %i", multinomial); *\/ */
/* /\* 	strcat(np[current_logo].name, tempstring); *\/ */
/* /\* 	/\\* Reverse_complement_normalized_pwm(&rc, &qp); *\/ */
/* /\* 	   for (counter = 0; counter < 4; counter++)  *\/ */
/* /\* 	   { *\/ */
/* /\* 	   printf("\n"); *\/ */
/* /\* 	   for (counter2 = 0; counter2 < rc.width; counter2++) printf ("\t%.1f", rc.fraction[counter][counter2]);  *\/ */
/* /\* 	   } *\\/ *\/ */
/* /\*       } *\/ */
/* /\*       current_logo++; *\/ */
/* /\*       Invert_pwm(&bqp, &qp, query_sequence_length / 3); *\/ */
/* /\*     } *\/ */

/* /\*     /\\* CALCULATES QUERY SEQUENCE VALUE IF FLANK IS USED *\\/ *\/ */
/* /\*     if (flank == 1) { *\/ */
/* /\*       if ((multinomial ^ flank_with_pwm) == 0 && iupac_query == 0) { *\/ */
/* /\* 	query_sequence_length = strlen(searchstring); *\/ */
/* /\* 	for (query_sequence_value = 0, position = 0, position_value = pow(4, query_sequence_length - 1); position < query_sequence_length; *\/ */
/* /\* 	     position++, position_value /= 4) { *\/ */
/* /\* 	  for (nucleotide_value = 0; nucleotide_value < 4 && searchstring[position] != forward[nucleotide_value]; nucleotide_value++) ; *\/ */
/* /\* 	  if (nucleotide_value == 4) { *\/ */
/* /\* 	    printf("\nSEED ERROR in round %i\n", round); *\/ */
/* /\* 	    exit(1); *\/ */
/* /\* 	  } *\/ */
/* /\* 	  query_sequence_value += position_value * nucleotide_value; *\/ */
/* /\* 	} *\/ */
/* /\*       } *\/ */
/* /\*       /\\* TESTS IF QUERY SEQUENCE IS A PALINDROME *\\/ *\/ */
/* /\*       searchstring_is_palindrome = Is_this_string_iupac_palindrome(searchstring, -1, -1); *\/ */
/* /\*       /\\* (Is_this_sequence_dimer(query_sequence_value, query_sequence_length) == 4); *\\/ *\/ */
/* /\*       if (searchstring_is_palindrome == 1) *\/ */
/* /\* 	printf("\nPalindrome"); *\/ */
/* /\*       if (extendedoutput == 1) *\/ */
/* /\* 	printf("\nsequence_value_ULL = %ju", (uintmax_t) query_sequence_value); *\/ */
/* /\*     } *\/ */

/* /\*     if (dinucleotide_properties == 1) *\/ */
/* /\*       one_hit_di[1].query_sequence_length = query_sequence_length; *\/ */

/* /\*     /\\* FILE MAIN LOOP *\\/ *\/ */
/* /\*     for (; file_number < 2; fclose(open_file), file_number++) { *\/ */
/* /\*       open_file = fopen(file_name[file_number], "r"); *\/ */
/* /\*       if (open_file == '\0') { *\/ */
/* /\* 	printf("File: %s not found\n\n", file_name[file_number]); *\/ */
/* /\* 	if (print_counts == 1) *\/ */
/* /\* 	  number_of_files = 1; *\/ */
/* /\* 	else *\/ */
/* /\* 	  exit(1); *\/ */
/* /\*       } *\/ */

/* /\*       /\\* CLEARS DIFFERENT SEQUENCES LIST *\\/ *\/ */
/* /\*       if (remove_non_unique == 1) { *\/ */
/* /\* 	if (extendedoutput == 1) *\/ */
/* /\* 	  printf("Clears sequence list\n"); *\/ */
/* /\* 	for (counter = 0; counter < 101; counter++) { *\/ */
/* /\* 	  /\\* printf("\t%i", counter); *\\/ *\/ */
/* /\* 	  unordered_list[counter].sequence_value_ULL = 1E9; *\/ */
/* /\* 	  unordered_list[counter].incidence = 0; *\/ */
/* /\* 	} *\/ */

/* /\* 	for (counter = 0; counter < max_number_of_sequences; counter++) { *\/ */
/* /\* 	  sorted_list[0][counter].sequence_value_ULL = 0; *\/ */
/* /\* 	  sorted_list[0][counter].incidence = 0; *\/ */
/* /\* 	  sorted_list[1][counter].sequence_value_ULL = 0; *\/ */
/* /\* 	  sorted_list[1][counter].incidence = 0; *\/ */
/* /\* 	} *\/ */

/* /\* 	number_of_unordered_sequences = 0; *\/ */
/* /\* 	for (counter = 0; counter < 4097; counter++) *\/ */
/* /\* 	  sorted_index[counter] = '\0'; *\/ */
/* /\* 	sorted_index[counter] = sorted_list[0]; *\/ */
/* /\*       } *\/ */
/* /\*       if (extendedoutput == 1) { *\/ */
/* /\* 	printf("Sequence list cleared\n"); *\/ */
/* /\* 	fflush(stdout); *\/ */
/* /\*       } *\/ */

/* /\*       /\\* SEQUENCE LINE LOADING LOOP *\\/ *\/ */
/* /\*       for (eof_reached = 0, line = 1;; line++) { *\/ */

/* /\*  start_of_line_loop: *\/ */

/* /\* 	/\\* TAKES ONE SEQUENCE FROM FILE *\\/ *\/ */
/* /\* 	for (charcounter = 0;;) { *\/ */
/* /\* 	  text1 = getc(open_file); *\/ */
/* /\* 	  /\\* printf("%c2", text1); *\\/ *\/ */
/* /\* 	  if (text1 == EOF) { *\/ */
/* /\* 	    eof_reached = 1; *\/ */
/* /\* 	    break; *\/ */
/* /\* 	  } *\/ */
/* /\* 	  if (text1 == '\n') *\/ */
/* /\* 	    break; *\/ */
/* /\* 	  if (text1 == 'A' || text1 == 'C' || text1 == 'G' || text1 == 'T' || text1 == 'N') {	/\\* ONLY ACCEPTS NUCLEOTIDES IN CAPS, N ONLY ACCEPTED SO THAT ERROR WILL BE DIFFERENT *\\/ *\/ */
/* /\* 	    current_sequence[charcounter] = text1; *\/ */
/* /\* 	    charcounter++; *\/ */
/* /\* 	  } *\/ */
/* /\* 	} *\/ */
/* /\* 	current_sequence[charcounter] = '\0'; *\/ */
/* /\* 	/\\* printf ("\nSequence: %s", current_sequence); *\\/ *\/ */
/* /\* 	if (eof_reached == 0 && strlen(current_sequence) != Nlength - 1) { *\/ */
/* /\* 	  printf("\nWrong sequence length on line %li", line); *\/ */
/* /\* 	  goto start_of_line_loop; *\/ */
/* /\* 	} *\/ */

/* /\* 	/\\* CHECKS IF AT END OF FILE *\\/ *\/ */
/* /\* 	if (eof_reached == 1) { *\/ */
/* /\* 	  printf("\nEOF encountered in file %i on line %li\n", file_number, line); *\/ */
/* /\* 	  break; *\/ */
/* /\* 	} *\/ */

/* /\* 	/\\* STRAND LOOP *\\/ *\/ */
/* /\* 	for (non_unique_flag = 0, strand = 0; strand < 2 && non_unique_flag == 0; strand++) { *\/ */

/* /\* 	  /\\* printf ("\tstrand %i", strand); *\\/ *\/ */

/* /\* 	  /\\* CALCULATES INTEGER VALUE CORRESPONDING TO SEQUENCE N-mer *\\/ *\/ */
/* /\* 	  if (strand == 0) { *\/ */
/* /\* 	    /\\* FORWARD STRAND *\\/ *\/ */
/* /\* 	    for (current_sequence_value_ULL = 0, position = 0, position_value = left_position_value; position < Nlength - 1; *\/ */
/* /\* 		 position++, position_value /= 4) { *\/ */
/* /\* 	      for (nucleotide_value = 0; nucleotide_value < 4 && current_sequence[position] != dnaforward[nucleotide_value]; nucleotide_value++) ; *\/ */
/* /\* 	      if (nucleotide_value == 4) {	/\\* printf("\nSEQUENCE ERROR AT POSITION %i, %i \n", line, position); *\\/ *\/ */
/* /\* 		goto start_of_line_loop; *\/ */
/* /\* 	      } else if (file_number == 0 && (multinomial && complex_background) == 0) *\/ */
/* /\* 		nc[strand].incidence[nucleotide_value][position]++;	/\\* COUNTS NUCLEOTIDE *\\/ *\/ */
/* /\* 	      current_sequence_value_ULL += position_value * nucleotide_value; *\/ */
/* /\* 	    } *\/ */
/* /\* 	    forward_sequence_value_ULL = current_sequence_value_ULL; *\/ */
/* /\* 	    /\\* printf("\nsequence_value_ULL = %llu", current_sequence_value_ULL); *\\/ *\/ */
/* /\* 	  } else { *\/ */
/* /\* 	    /\\* REVERSE STRAND *\\/ *\/ */
/* /\* 	    for (current_sequence_value_ULL = 0, position = Nlength - 2, position_value = left_position_value; position > -1; *\/ */
/* /\* 		 position--, position_value /= 4) { *\/ */
/* /\* 	      for (nucleotide_value = 0; nucleotide_value < 4 && current_sequence[position] != dnareverse[nucleotide_value]; nucleotide_value++) ; *\/ */
/* /\* 	      if (nucleotide_value == 4) *\/ */
/* /\* 		printf("\nSEQUENCE ERROR AT POSITION %li, %i \n", line, position); *\/ */
/* /\* 	      else if (file_number == 0 && (multinomial && complex_background) == 0) *\/ */
/* /\* 		nc[strand].incidence[nucleotide_value][Nlength - position - 2]++;	/\\* COUNTS NUCLEOTIDE *\\/ *\/ */
/* /\* 	      current_sequence_value_ULL += position_value * nucleotide_value; *\/ */
/* /\* 	    } *\/ */
/* /\* 	  } *\/ */

/* /\* 	  /\\* printf("\nsequence, value, sequence : %s\t%ju\t", current_sequence, current_sequence_value_ULL); *\/ */
/* /\* 	     for(position = Nlength-2; position > -1 ; position--) printf("%c", forward[(current_sequence_value_ULL & mask_ULL[1][position]) >> (position * 2)]); *\\/ *\/ */

/* /\* 	  /\\* IF REMOVE NON UNIQUE IS ONE, INDEXES FULL-LENGTH SEQUENCES TO FIND IF SAME SEQUENCE REOCCURS *\\/ *\/ */
/* /\* 	  if (strand == 0 && remove_non_unique == 1) { *\/ */
/* /\* 	    for (shortcounter = 0; shortcounter < number_of_unordered_sequences; shortcounter++) { *\/ */
/* /\* 	      /\\* printf ("\nComparing %i to %i", current_sequence_value_ULL, unordered_list[shortcounter].sequence_value_ULL); *\\/ *\/ */
/* /\* 	      if (current_sequence_value_ULL == unordered_list[shortcounter].sequence_value_ULL) { *\/ */
/* /\* 		non_unique_flag = 1; *\/ */
/* /\* 		break; *\/ */
/* /\* 	      } *\/ */
/* /\* 	    } *\/ */
/* /\* 	    if (non_unique_flag == 1) { *\/ */
/* /\* 	      unordered_list[shortcounter].incidence++; *\/ */
/* /\* 	      /\\* printf("\nUNORDERED LIST HIT: Sequence %i found %i times", unordered_list[shortcounter].sequence_value_ULL, unordered_list[shortcounter].incidence); *\\/ *\/ */
/* /\* 	      continue; *\/ */
/* /\* 	    } *\/ */
/* /\* 	    unordered_list[number_of_unordered_sequences].sequence_value_ULL = current_sequence_value_ULL; *\/ */
/* /\* 	    unordered_list[number_of_unordered_sequences].incidence = 1; *\/ */
/* /\* 	    index_position = sorted_index[current_sequence_value_ULL / 65536]; *\/ */
/* /\* 	    if (index_position != '\0') { *\/ */
/* /\* 	      /\\* printf("\nLOOKING IF %i IS NON-UNIQUE IN SORTED LIST, STARTING FROM %i", current_sequence_value_ULL, (*sorted_index[current_sequence_value_ULL / 65536]).sequence_value_ULL); *\/ */
/* /\* 	         printf("\nRANGE %i", sorted_index[4097] - index_position); *\\/ *\/ */
/* /\* 	      for (; index_position < sorted_index[4097]; index_position++) { *\/ */
/* /\* 		/\\* printf ("\nzzComparing %i to %i", current_sequence_value_ULL, (*index_position).sequence_value_ULL); *\\/ *\/ */
/* /\* 		if (current_sequence_value_ULL == (*index_position).sequence_value_ULL) { *\/ */
/* /\* 		  non_unique_flag = 1; *\/ */
/* /\* 		  break; *\/ */
/* /\* 		} *\/ */
/* /\* 		if (current_sequence_value_ULL < (*index_position).sequence_value_ULL) *\/ */
/* /\* 		  break; *\/ */
/* /\* 	      } *\/ */
/* /\* 	      if (non_unique_flag == 1) { *\/ */
/* /\* 		(*index_position).incidence++; *\/ */
/* /\* 		/\\* printf("\nSORTED LIST HIT: Sequence %i found %i times", (*index_position).sequence_value_ULL, (*index_position).incidence); *\\/ *\/ */
/* /\* 		continue; *\/ */
/* /\* 	      } *\/ */
/* /\* 	    } *\/ */
/* /\* 	    if (number_of_unordered_sequences > 99) { *\/ */
/* /\* 	      /\\* printf("100 found-"); *\\/ *\/ */
/* /\* 	      number_of_unordered_sequences = *\/ */
/* /\* 		  Sort_and_index(sorted_list, unordered_list, current_sorted_list, number_of_unordered_sequences, sorted_index); *\/ */
/* /\* 	      if (current_sorted_list == 1) *\/ */
/* /\* 		current_sorted_list = 0; *\/ */
/* /\* 	      else *\/ */
/* /\* 		current_sorted_list = 1; *\/ */
/* /\* 	      /\\* printf("\nCurrent_number_of_different_sequences: %i", sorted_index[4097] - sorted_list[current_sorted_list]); *\\/ *\/ */
/* /\* 	    } else *\/ */
/* /\* 	      number_of_unordered_sequences++; *\/ */
/* /\* 	  } *\/ */

/* /\* 	  if (extendedoutput == 1) { *\/ */
/* /\* 	    printf("\nLINE %li, STRAND %i, integer_value %ju, sequence: ", line, strand, (uintmax_t) current_sequence_value_ULL); *\/ */
/* /\* 	    /\\* SEQUENCE PRINT *\\/ *\/ */
/* /\* 	    for (position = Nlength - 2; position > -1; position--) *\/ */
/* /\* 	      printf("%c", forward[(current_sequence_value_ULL & mask_ULL[1][position]) >> (position * 2)]); *\/ */
/* /\* 	  } *\/ */

/* /\* 	  /\\* only counts nucleotides from first file if flank is used  *\/ */
/* /\* 	     if (flank == 1 && file_number == 0) continue; *\\/ *\/ */

/* /\* 	  if (flank == 1) { *\/ */
/* /\* 	    /\\* FLANK TOOL *\\/ *\/ */

/* /\* 	    /\\* FINDS IF THIS SEQUENCE HAS MATCH IN EITHER STRAND *\\/ *\/ */
/* /\* 	    if (strand == 0) { *\/ */
/* /\* 	      /\\* COUNTS MATCHES TO FILTER IF FILTER SET *\\/ *\/ */
/* /\* 	      if (match_filter != 0) { *\/ */
/* /\* 		Findpwmmatch(&mfp, filter_cut_off, current_sequence_value_ULL, &match[0]); *\/ */
/* /\* 		Findpwmmatch(&mfp, filter_cut_off, Reverse_complement_sequence_value(current_sequence_value_ULL, Nlength - 1), &match[1]); *\/ */
/* /\* 		Remove_palindromic_matches(&match, query_sequence_length); *\/ */
/* /\* 		matches_to_filter = match[0].position[0] + match[1].position[0]; *\/ */
/* /\* 	      } *\/ */

/* /\* 	      /\\* COUNTS MATCHES TO MULTINOMIAL 2 IF DINUC SET AND MULTINOMIAL IS LESS THAN 2 *\\/ *\/ */
/* /\* 	      if (multinomial_2_dinuc != 0) { *\/ */
/* /\* 		Findpwmmatch(&mmp, multi_2_cut_off, current_sequence_value_ULL, &dinucmatch[0]); *\/ */
/* /\* 		Findpwmmatch(&mmp, multi_2_cut_off, Reverse_complement_sequence_value(current_sequence_value_ULL, Nlength - 1), &dinucmatch[1]); *\/ */
/* /\* 		Remove_palindromic_matches(&dinucmatch, query_sequence_length); *\/ */
/* /\* 		if (dinucmatch[0].position[0] == 1 && dinucmatch[1].position[0] == 0) *\/ */
/* /\* 		  Multinomial_add_to_dinucleotide_matrix(&one_hit_di[file_number], &di, &dep[file_number], &qp, dinucmatch[0].position[1], *\/ */
/* /\* 							 dinucmatch[0].score[1], multi_2_cut_off, forward_sequence_value_ULL); *\/ */
/* /\* 		if (dinucmatch[0].position[0] == 0 && dinucmatch[1].position[0] == 1) *\/ */
/* /\* 		  Multinomial_add_to_dinucleotide_matrix(&one_hit_di[file_number], &di, &dep[file_number], &qp, dinucmatch[1].position[1], *\/ */
/* /\* 							 dinucmatch[1].score[1], multi_2_cut_off, *\/ */
/* /\* 							 Reverse_complement_sequence_value(current_sequence_value_ULL, Nlength - 1)); *\/ */
/* /\* 	      } *\/ */

/* /\* 	      /\\* FINDS number_of_matches of QUERY SEQUENCE OR PWM in CURRENT SEQUENCE *\\/ *\/ */
/* /\* 	      if (flank_with_pwm == 0) { *\/ */
/* /\* 		number_of_matches[0] = Findexactmatch(query_sequence_value, query_sequence_length, current_sequence_value_ULL, &match[0]); *\/ */
/* /\* 		number_of_matches[1] = *\/ */
/* /\* 		    Findexactmatch(query_sequence_value, query_sequence_length, *\/ */
/* /\* 				   Reverse_complement_sequence_value(current_sequence_value_ULL, Nlength - 1), &match[1]); *\/ */
/* /\* 	      } else { *\/ */
/* /\* 		number_of_matches[0] = Findpwmmatch(&qp, cut_off, current_sequence_value_ULL, &match[0]); *\/ */
/* /\* 		number_of_matches[1] = *\/ */
/* /\* 		    Findpwmmatch(&qp, cut_off, Reverse_complement_sequence_value(current_sequence_value_ULL, Nlength - 1), &match[1]); *\/ */
/* /\* 	      } *\/ */
/* /\* 	      if (number_of_matches[0] + number_of_matches[1] > 0) { *\/ */
/* /\* 		current_sequence_contains_match = 1; *\/ */
/* /\* 	      } else { *\/ */
/* /\* 		current_sequence_contains_match = 0; *\/ */
/* /\* 		number_of_sequences_with_no_hits[file_number]++; *\/ */
/* /\* 	      } *\/ */
/* /\* 	    } *\/ */

/* /\* 	    /\\* EXCLUDES ENTIRE STRAND IF INDICATED *\\/ *\/ */
/* /\* 	    if ((limit_hits_to_strand == 'F' || limit_hits_to_strand == 'f') && limit_hits_to_position == 0) *\/ */
/* /\* 	      number_of_matches[1] = 0; *\/ */
/* /\* 	    if ((limit_hits_to_strand == 'R' || limit_hits_to_strand == 'r') && limit_hits_to_position == 0) *\/ */
/* /\* 	      number_of_matches[0] = 0; *\/ */

/* /\* 	    /\\* REMOVES MATCHES FROM REVERSE STRAND IF QUERY SEQUENCE IS A PALINDROME *\/ */
/* /\* 	       if (searchstring_is_palindrome == 1) number_of_matches[1] = 0; *\\/ *\/ */

/* /\* 	    /\\* printf("\n%i number_of_matches in strand %i at positions: ", number_of_matches[strand], strand); *\/ */
/* /\* 	       for (counter = 0; counter < number_of_matches[strand]; counter++)  *\/ */
/* /\* 	       { *\/ */
/* /\* 	       if (strand == 0) printf ("%i\t", match.position[strand][counter]); *\/ */
/* /\* 	       else printf ("%i\t", Nlength - match.position[strand][counter] - query_sequence_length + 1); *\/ */
/* /\* 	       } *\\/ *\/ */

/* /\* 	    /\\* CHECKS IF HIT POSITION & STRAND LIMITED BY USER AND FINDS IF ANY OF THE HITS number_of_matches MATCH USER DEFINED POSITION(S) *\\/ *\/ */
/* /\* 	    if (limit_hits_to_position != 0) { *\/ */
/* /\* 	      correct_position = 0; *\/ */
/* /\* 	      if (limit_hits_to_strand != 'R' && limit_hits_to_strand != 'r') *\/ */
/* /\* 		for (counter = 1; counter <= number_of_matches[0]; counter++) *\/ */
/* /\* 		  if (match[0].position[counter] == limit_hits_to_position) *\/ */
/* /\* 		    correct_position = 1; *\/ */
/* /\* 	      if (limit_hits_to_strand != 'F' && limit_hits_to_strand != 'f') *\/ */
/* /\* 		for (counter = 1; counter <= number_of_matches[1]; counter++) *\/ */
/* /\* 		  if (Nlength - query_sequence_length + 1 - match[1].position[counter] == limit_hits_to_position) *\/ */
/* /\* 		    correct_position = 1; *\/ */
/* /\* 	    } *\/ */

/* /\* 	    /\\* CHECKS IF HIT POSITION & STRAND EXCLUDED BY USER AND FINDS IF ANY OF THE HITS IN EXCLUDED POSITIONS *\\/ *\/ */
/* /\* 	    if (exclude_position[0] != 0) { *\/ */
/* /\* 	      correct_position = 1; *\/ */
/* /\* 	      for (counter2 = 0; exclude_position[counter2] != 0; counter2++) { *\/ */
/* /\* 		if (exclude_strand[counter2] != 'R' && exclude_strand[counter2] != 'r') *\/ */
/* /\* 		  for (counter = 1; counter <= number_of_matches[0]; counter++) *\/ */
/* /\* 		    if (match[0].position[counter] == exclude_position[counter2]) *\/ */
/* /\* 		      correct_position = 0; *\/ */
/* /\* 		if (exclude_strand[counter2] != 'F' && exclude_strand[counter2] != 'f') *\/ */
/* /\* 		  for (counter = 1; counter <= number_of_matches[1]; counter++) *\/ */
/* /\* 		    if (Nlength - query_sequence_length + 1 - match[1].position[counter] == exclude_position[counter2]) *\/ */
/* /\* 		      correct_position = 0; *\/ */
/* /\* 	      } *\/ */

/* /\* 	    } *\/ */

/* /\* 	    /\\* REMOVES PALINDROMIC MATCHES *\\/ *\/ */
/* /\* 	    if (strand == 0 && number_of_matches[0] > 0 && number_of_matches[1] > 0) { *\/ */
/* /\* 	      Remove_palindromic_matches(&match, query_sequence_length); *\/ */
/* /\* 	      palindromic_hits[file_number] += number_of_matches[0] + number_of_matches[1] - match[0].position[0] - match[1].position[0]; *\/ */
/* /\* 	      number_of_matches[0] = match[0].position[0]; *\/ */
/* /\* 	      number_of_matches[1] = match[1].position[0]; *\/ */
/* /\* 	    } *\/ */
/* /\* 	    if (match_filter == 0) *\/ */
/* /\* 	      matches_to_filter = number_of_matches[0] + number_of_matches[1];	/\\* SEPARATE FILTER IS NOT SET *\\/ *\/ */

/* /\* 	    if (strand == 1 && correct_position == 1) { *\/ */
/* /\* 	      if (number_of_matches[0] + number_of_matches[1] > 0) *\/ */
/* /\* 		number_of_sequences_with_hits[file_number]++; *\/ */
/* /\* 	      /\\* FLANK (ONLY SEQUENCES WITH ONE OCCURRENCE) *\\/ *\/ */
/* /\* 	      /\\* ADDS NUCLEOTIDES FROM number_of_matches TO ONE HIT RESULT PWM *\\/ *\/ */

/* /\* 	      /\\* ONE HIT ONLY (MATCHES TO FILTER IS 1) *\\/ *\/ */
/* /\* 	      if (matches_to_filter == 1) { *\/ */
/* /\* 		if (number_of_matches[1] == 0 && number_of_matches[0] == 1) { *\/ */
/* /\* 		  /\\* ADDS FLANK KMER *\\/ *\/ */
/* /\* 		  if (flank_kmer_pos != -100 && (Nlength - match[0].position[1]) - flank_kmer_pos < Nlength *\/ */
/* /\* 		      && (Nlength - match[0].position[1]) - flank_kmer_pos > 3) *\/ */
/* /\* 		    flank_kmer_count[file_number][(forward_sequence_value_ULL >> (((Nlength - match[0].position[1]) - flank_kmer_pos - 4) * 2)) & *\/ */
/* /\* 						  255]++; *\/ */

/* /\* 		  if (multinomial == 0) { *\/ */
/* /\* 		    Add_to_pwm(&one_hit_pwm[file_number], match[0].position[1], forward_sequence_value_ULL, &background_pwm[0]); *\/ */
/* /\* 		    if (dinucleotide_properties == 1 && multinomial_2_dinuc == 0) *\/ */
/* /\* 		      Add_to_dinucleotide_matrix(&one_hit_di[file_number], &di, &dep[file_number], match[0].position[1], forward_sequence_value_ULL); *\/ */
/* /\* 		  } else { *\/ */
/* /\* 		    Multinomial_add_to_pwm(&one_hit_pwm[file_number], &qp, match[0].position[1], match[0].score[1], cut_off, *\/ */
/* /\* 					   forward_sequence_value_ULL, &background_pwm[0]); *\/ */
/* /\* 		    if (dinucleotide_properties == 1 && multinomial_2_dinuc == 0) *\/ */
/* /\* 		      Multinomial_add_to_dinucleotide_matrix(&one_hit_di[file_number], &di, &dep[file_number], &qp, match[0].position[1], *\/ */
/* /\* 							     match[0].score[1], cut_off, forward_sequence_value_ULL); *\/ */
/* /\* 		  } *\/ */
/* /\* 		} *\/ */

/* /\* 		if (number_of_matches[0] == 0 && number_of_matches[1] == 1) { *\/ */
/* /\* 		  /\\* ADDS FLANK KMER *\\/ *\/ */
/* /\* 		  if (flank_kmer_pos != -100 && (Nlength - match[1].position[1]) - flank_kmer_pos < Nlength *\/ */
/* /\* 		      && (Nlength - match[1].position[1]) - flank_kmer_pos > 3) *\/ */
/* /\* 		    flank_kmer_count[file_number][(current_sequence_value_ULL >> (((Nlength - match[1].position[1]) - flank_kmer_pos - 4) * 2)) & *\/ */
/* /\* 						  255]++; *\/ */

/* /\* 		  if (multinomial == 0) { *\/ */
/* /\* 		    Add_to_pwm(&one_hit_pwm[file_number], match[1].position[1], current_sequence_value_ULL, &background_pwm[1]); *\/ */
/* /\* 		    if (dinucleotide_properties == 1 && multinomial_2_dinuc == 0) *\/ */
/* /\* 		      Add_to_dinucleotide_matrix(&one_hit_di[file_number], &di, &dep[file_number], match[1].position[1], current_sequence_value_ULL); *\/ */
/* /\* 		  } else { *\/ */
/* /\* 		    Multinomial_add_to_pwm(&one_hit_pwm[file_number], &qp, match[1].position[1], match[1].score[1], cut_off, *\/ */
/* /\* 					   current_sequence_value_ULL, &background_pwm[1]); *\/ */
/* /\* 		    if (dinucleotide_properties == 1 && multinomial_2_dinuc == 0) *\/ */
/* /\* 		      Multinomial_add_to_dinucleotide_matrix(&one_hit_di[file_number], &di, &dep[file_number], &qp, match[1].position[1], *\/ */
/* /\* 							     match[1].score[1], cut_off, current_sequence_value_ULL); *\/ */
/* /\* 		  } *\/ */
/* /\* 		} *\/ */
/* /\* 	      } *\/ */
/* /\* 	      /\\* END OF ONE HIT ONLY *\\/ *\/ */

/* /\* 	      /\\* FLANK (SEQUENCES WITH TWO OCCURRENCES) *\\/ *\/ */
/* /\* 	      if ((number_of_matches[0] == 2 || number_of_matches[1] == 2) && number_of_matches[0] + number_of_matches[1] == 2 && file_number == 1) { *\/ */
/* /\* 		orientation = head_to_tail; *\/ */
/* /\* 		if (number_of_matches[0] == 2) { *\/ */
/* /\* 		  spacing = match[0].position[1] - match[0].position[2]; *\/ */
/* /\* 		  if (spacing > 49) *\/ */
/* /\* 		    Exit_with_error("spacing too long", error_values); *\/ */
/* /\* 		  first_match_position = match[0].position[2]; *\/ */
/* /\* 		  first_sequence_value_ULL = forward_sequence_value_ULL; *\/ */
/* /\* 		  /\\* printf("\nSame orientation in forward strand "); *\\/ *\/ */
/* /\* 		  Add_to_pwm(&two_hits_pwm[file_number][orientation][spacing], first_match_position, forward_sequence_value_ULL, &background_pwm[0]); *\/ */
/* /\* 		} else { *\/ */
/* /\* 		  spacing = match[1].position[1] - match[1].position[2]; *\/ */
/* /\* 		  first_match_position = match[1].position[2]; *\/ */
/* /\* 		  first_sequence_value_ULL = current_sequence_value_ULL; *\/ */
/* /\* 		  /\\* printf("\nSame orientation in reverse strand "); *\\/ *\/ */
/* /\* 		  Add_to_pwm(&two_hits_pwm[file_number][orientation][spacing], first_match_position, first_sequence_value_ULL, &background_pwm[1]); *\/ */
/* /\* 		} *\/ */
/* /\* 		two_hits_connecting_matrix.incidence[orientation][spacing]++; *\/ */
/* /\* 		two_hits_connecting_matrix.two_hit_matches += 2; *\/ */
/* /\* 		/\\* printf ("%s, spacing %i, %i:th occurrence, position %i", orientation_string[orientation], spacing, two_hits_pwm[file_number][orientation][spacing].max_counts+1, first_match_position); *\\/ *\/ */

/* /\* 	      } *\/ */

/* /\* 	      if (number_of_matches[0] == 1 && number_of_matches[1] == 1 && file_number == 1) { *\/ */
/* /\* 		first_match_position = Nlength - match[1].position[1] - query_sequence_length + 1; *\/ */
/* /\* 		/\\* printf("\nOpposite orientation, values forward %i, reverse %i", match[0].position[1], first_match_position); *\\/ *\/ */
/* /\* 		spacing = match[0].position[1] - first_match_position; *\/ */
/* /\* 		if (spacing > 0) { *\/ */
/* /\* 		  orientation = tail_to_tail; *\/ */
/* /\* 		  first_sequence_value_ULL = forward_sequence_value_ULL; *\/ */
/* /\* 		} else { *\/ */
/* /\* 		  orientation = head_to_head; *\/ */
/* /\* 		  spacing = first_match_position - match[0].position[1]; *\/ */
/* /\* 		  first_match_position = match[0].position[1]; *\/ */
/* /\* 		  first_sequence_value_ULL = forward_sequence_value_ULL; *\/ */
/* /\* 		} *\/ */
/* /\* 		/\\* printf ("\nOpposite orientation %s, spacing %i, %i:th occurrence, position %i", orientation_string[orientation], spacing, two_hits_pwm[file_number][orientation][spacing].max_counts+1, first_match_position); *\\/ *\/ */
/* /\* 		Add_to_pwm(&two_hits_pwm[file_number][orientation][spacing], first_match_position, first_sequence_value_ULL, &background_pwm[0]); *\/ */
/* /\* 		two_hits_connecting_matrix.incidence[orientation][spacing]++; *\/ */
/* /\* 		two_hits_connecting_matrix.two_hit_matches += 2; *\/ */
/* /\* 	      } *\/ */

/* /\* 	      /\\* SPACING AND ORIENTATION *\\/ *\/ */
/* /\* 	      /\\* FILLS CONNECTING MATRIX *\\/ *\/ */
/* /\* 	      /\\* HEAD TO TAIL (SAME) SPACINGS AND PWM FILL *\\/ *\/ */
/* /\* 	      /\\* FORWARD STRAND *\\/ *\/ */
/* /\* 	      if (align_matches == 1) { *\/ */
/* /\* 		max_align_score[0][0] = 0; *\/ */
/* /\* 		max_align_score[1][0] = 0; *\/ */
/* /\* 		max_align_score[0][1] = 0; *\/ */
/* /\* 		max_align_score[1][1] = 0; *\/ */
/* /\* 	      } *\/ */
/* /\* 	      for (counter = 1; counter <= number_of_matches[0]; counter++) { *\/ */
/* /\* 		if (align_matches == 1) *\/ */
/* /\* 		  Add_to_alignscore(forward_sequence_value_ULL, &all_hits_align_scores[file_number], &current_align_score, match[0].position[counter], *\/ */
/* /\* 				    max_align_score[file_number]); *\/ */
/* /\* 		if (multinomial == 0) *\/ */
/* /\* 		  Add_to_pwm(&all_hits_pwm[file_number], match[0].position[counter], forward_sequence_value_ULL, &background_pwm[0]); *\/ */
/* /\* 		else *\/ */
/* /\* 		  Multinomial_add_to_pwm(&all_hits_pwm[file_number], &qp, match[0].position[counter], match[0].score[counter], cut_off, *\/ */
/* /\* 					 forward_sequence_value_ULL, &background_pwm[0]); *\/ */
/* /\* 	      } *\/ */
/* /\* 	      if (file_number == 1) *\/ */
/* /\* 		for (counter = 1; counter <= number_of_matches[0]; counter++) { *\/ */
/* /\* 		  hit_position.incidence[0][match[0].position[counter]]++; *\/ */
/* /\* 		  for (counter2 = counter; counter2 <= number_of_matches[0]; counter2++) *\/ */
/* /\* 		    cm.incidence[head_to_tail][(match[0].position[counter] - match[0].position[counter2])]++; *\/ */
/* /\* 		} *\/ */
/* /\* 	      /\\* REVERSE STRAND *\\/ *\/ */
/* /\* 	      for (counter = 1; counter <= number_of_matches[1]; counter++) { *\/ */
/* /\* 		if (align_matches == 1) *\/ */
/* /\* 		  Add_to_alignscore(current_sequence_value_ULL, &all_hits_align_scores[file_number], &current_align_score, match[1].position[counter], *\/ */
/* /\* 				    max_align_score[file_number]); *\/ */
/* /\* 		if (multinomial == 0) *\/ */
/* /\* 		  Add_to_pwm(&all_hits_pwm[file_number], match[1].position[counter], current_sequence_value_ULL, &background_pwm[1]); *\/ */
/* /\* 		else *\/ */
/* /\* 		  Multinomial_add_to_pwm(&all_hits_pwm[file_number], &qp, match[1].position[counter], match[1].score[counter], cut_off, *\/ */
/* /\* 					 current_sequence_value_ULL, &background_pwm[1]); *\/ */
/* /\* 	      } *\/ */
/* /\* 	      if (file_number == 1) *\/ */
/* /\* 		for (counter = 1; counter <= number_of_matches[1]; counter++) { *\/ */
/* /\* 		  hit_position.incidence[1][match[1].position[counter] - 1]++; *\/ */
/* /\* 		  for (counter2 = counter; counter2 <= number_of_matches[1]; counter2++) *\/ */
/* /\* 		    cm.incidence[head_to_tail][(match[1].position[counter] - match[1].position[counter2])]++; *\/ */
/* /\* 		} *\/ */
/* /\* 	      /\\* HEAD TO HEAD and TAIL TO TAIL (OPPOSITE) *\\/ *\/ */
/* /\* 	      if (file_number == 1) *\/ */
/* /\* 		for (counter = 1; counter <= number_of_matches[0]; counter++) *\/ */
/* /\* 		  for (counter2 = 1; counter2 <= number_of_matches[1]; counter2++) { *\/ */
/* /\* 		    reverse_strand_position = Nlength - match[1].position[counter2] - query_sequence_length + 1; *\/ */
/* /\* 		    if (match[0].position[counter] >= reverse_strand_position) *\/ */
/* /\* 		      cm.incidence[tail_to_tail][(match[0].position[counter] - reverse_strand_position)]++; *\/ */
/* /\* 		    if (match[0].position[counter] <= reverse_strand_position) *\/ */
/* /\* 		      cm.incidence[head_to_head][(reverse_strand_position - match[0].position[counter])]++; *\/ */
/* /\* 		  } *\/ */
/* /\* 	      if (align_matches == 1 && max_align_score[file_number][0] > 0) *\/ */
/* /\* 		align_score_histogram[file_number][0][max_align_score[file_number][0]]++; *\/ */
/* /\* 	      if (align_matches == 1 && max_align_score[file_number][1] > 0) *\/ */
/* /\* 		align_score_histogram[file_number][1][max_align_score[file_number][1]]++; *\/ */
/* /\* 	    } *\/ */

/* /\* 	    /\\* END OF FLANK TOOL *\\/ *\/ */
/* /\* 	  } *\/ */

/* /\* 	  /\\* KMER COUNT TOOL *\\/ *\/ */
/* /\* 	  if (kmer_count == 1) { *\/ */
/* /\* 	    /\\* COUNTS KMERS *\\/ *\/ */

/* /\* 	    /\\* KMERS WITH NO GAPS *\\/ *\/ */

/* /\* 	    for (current_kmer_length = shortest_kmer; current_kmer_length < too_long_kmer; current_kmer_length++) { *\/ */
/* /\* 	      /\\* printf("\n\nSequence length: %i", current_kmer_length); *\\/ *\/ */
/* /\* 	      end_position = Nlength - current_kmer_length; *\/ */
/* /\* 	      for (position = 0; position < end_position; position++) { *\/ */
/* /\* 		/\\* current_value = results[current_kmer_length][0][0]; *\\/ *\/ */

/* /\* 		results[file_number][current_kmer_length][0][0][(current_sequence_value_ULL & mask_ULL[current_kmer_length][position]) >> *\/ */
/* /\* 								(position * 2)]++; *\/ */
/* /\* 		/\\* printf("\n%i", (current_sequence_value_ULL & mask_ULL[current_kmer_length][position]) >> (position * 2)); *\\/ *\/ */
/* /\* 		/\\* fflush(stdout); *\\/ *\/ */
/* /\* 	      } *\/ */
/* /\* 	      if (strand == 1 && current_sequence_contains_match == 0) *\/ */
/* /\* 		for (position = 0; position < end_position; position++) { *\/ */
/* /\* 		  /\\* COUNTS KMERS FROM NON-MATCHING SEQUENCE *\\/ *\/ */
/* /\* 		  results[file_number][current_kmer_length][1][0][(forward_sequence_value_ULL & mask_ULL[current_kmer_length][position]) >> *\/ */
/* /\* 								  (position * 2)]++; *\/ */
/* /\* 		  results[file_number][current_kmer_length][1][0][(current_sequence_value_ULL & mask_ULL[current_kmer_length][position]) >> *\/ */
/* /\* 								  (position * 2)]++; *\/ */
/* /\* 		} *\/ */
/* /\* 	    } *\/ */
/* /\* 	  } *\/ */

/* /\* 	  if (kmer_count == 1 && no_insertions == 0 && count_also_spaced_kmers != 0) { *\/ */
/* /\* 	    /\\* KMERS WITH GAPS *\\/ *\/ */
/* /\* 	    /\\* GENERATES DELETIONS INTO CURRENT SEQUENCE *\\/ *\/ */
/* /\* 	    for (deletion_size = 1; deletion_size < Nlength - shortest_kmer; deletion_size++) { *\/ */
/* /\* 	      Nmer_position = 1; *\/ */
/* /\* 	      too_long_nmer_position = Nlength - 1 - deletion_size; *\/ */
/* /\* 	      for (; Nmer_position < too_long_nmer_position; Nmer_position++) { *\/ */
/* /\* 		deleted_sequence_value_ULL = *\/ */
/* /\* 		    (current_sequence_value_ULL & lowmask_ULL[Nmer_position - 1]) ^ *\/ */
/* /\* 		    ((current_sequence_value_ULL & highmask_ULL[Nmer_position + deletion_size]) >> (deletion_size * 2)); *\/ */
/* /\* 		/\\* printf ("\n\n%i\t%i\t%i", deletion_size, Nmer_position, deleted_sequence_value_ULL); *\\/ *\/ */
/* /\* 		/\\* SEQUENCE PRINT *\\/ *\/ */
/* /\* 		/\\* for(position = Nlength-2; position > -1 ; position--) printf("%c", forward[(deleted_sequence_value_ULL & mask_ULL[1][position]) >> (position * 2)]); *\\/ *\/ */

/* /\* 		/\\* FINDS KMERS WITH GAPS FROM DELETED SEQUENCE *\\/ *\/ */
/* /\* 		for (current_kmer_length = shortest_kmer; current_kmer_length < too_long_kmer; current_kmer_length++) { *\/ */
/* /\* 		  if (count_also_spaced_kmers == 1) { *\/ */
/* /\* 		    position = Nmer_position - current_kmer_length + current_kmer_length / 2; *\/ */
/* /\* 		    end_position = position + 1 + current_kmer_length % 2; *\/ */
/* /\* 		  } else { */
/* 		    position = Nmer_position - current_kmer_length + 1; */
/* 		    end_position = Nmer_position; */
/* 		  } */
/* 		  if (position < 0) */
/* 		    position = 0; */
/* 		  if (end_position > Nlength - current_kmer_length - deletion_size) */
/* 		    end_position = Nlength - current_kmer_length - deletion_size; */
/* 		  for (; position < end_position; position++) { */
/* 		    results[file_number][current_kmer_length][current_kmer_length - Nmer_position + */
/* 							      position][deletion_size][((deleted_sequence_value_ULL & */
/* 											 mask_ULL[current_kmer_length][position])) >> (position * */
/* 																       2)]++; */
/* 		  } */

/* 		} */

/* 	      } */
/* 	    } */

/* 	    /\\* END OF KMER COUNT TOOL *\\/ */
/* 	  } */

/* 	  /\\* END OF STRAND LOOP *\\/ */
/* 	} */

/* 	/\\* END OF LINE LOOP *\\/ */
/*       } */

/*       if (remove_non_unique == 1) { */
/* 	number_of_unordered_sequences = Sort_and_index(sorted_list, unordered_list, current_sorted_list, number_of_unordered_sequences, sorted_index); */
/* 	if (current_sorted_list == 1) */
/* 	  current_sorted_list = 0; */
/* 	else */
/* 	  current_sorted_list = 1; */
/* 	/\\* printf("\nCurrent_number_of_different_sequences , unordered: %i , %i", sorted_index[4097] - sorted_list[current_sorted_list], number_of_unordered_sequences); *\\/ */
/* 	number_of_sequences_analyzed[file_number] = (double)((sorted_index[4097] - sorted_list[current_sorted_list]) + number_of_unordered_sequences); */
/* 	/\\* (sorted_index[4097] - sorted_list[current_sorted_list]) + number_of_unordered_sequences; *\\/ */
/*       } else */
/* 	number_of_sequences_analyzed[file_number] = (double)(line - 1); */
/*       printf("\nNumber of sequences in file %i: %.0f", file_number, number_of_sequences_analyzed[file_number]); */

/*       /\\* PRINTS INPUT SEQUENCES AND THEIR INCIDENCE *\\/ */
/*       if (print_input_sequences == 1 && remove_non_unique == 1) { */
/* 	printf("\n\nUnique_sequences: %.2f\nSequence", number_of_sequences_analyzed[file_number]); */
/* 	if (remove_non_unique == 1) */
/* 	  printf("\tIncidence"); */
/* 	for (counter = 0; counter < number_of_sequences_analyzed[file_number]; counter++) { */
/* 	  printf("\n"); */
/* 	  for (position = Nlength - 2; position > -1; position--) */
/* 	    printf("%c", forward[(sorted_list[1 - current_sorted_list][counter].sequence_value_ULL & mask_ULL[1][position]) >> (position * 2)]); */
/* 	  if (remove_non_unique == 1) */
/* 	    printf("\t%li", sorted_list[1 - current_sorted_list][counter].incidence); */
/* 	} */

/*       } */

/*       /\\* printf ("\n1-%i", np[0].width); *\\/ */

/*       /\\* CALCULATES BACKGROUND NUCLEOTIDE DISTRIBUTION FROM FIRST FILE *\\/ */
/*       if (file_number == 0 && even_background == 0 && (multinomial && complex_background) == 0) { */
/* 	Count_to_normalized_pwm(&background_pwm[0], &nc[0]); */
/* 	Count_to_normalized_pwm(&background_pwm[1], &nc[1]); */
/*       } */

/*       /\\* printf ("\n2-%i", np[0].width); *\\/ */
/*       /\\* END OF MAIN (FILE) LOOP *\\/ */
/*     } */

/*  exit_with_style: */

/*     if (flank == 1) { */

/*       if (round == 1) { */
/* 	/\\* ESTIMATES LAMBDA FROM KMER COUNT DATA (TAKES TOTAL COUNT OF LOWER HALF OF default 8-MERS can be set by -kl=[kmer length]) *\\/ */
/* 	number_of_kmers = pow(4, shortest_kmer); */
/* 	current_kmer_length = shortest_kmer; */

/* 	signal_kmer_count_p = results[1][current_kmer_length][0][0]; */
/* 	background_kmer_count_p = results[0][current_kmer_length][0][0]; */

/* 	/\\* COPIES KMERS FOR QSORT IF KMERS NEEDED LATER *\\/ */
/* 	if (kmer_count == 1 || expected_observed_plot == 1 || xyplot == 1) { */
/* 	  signal_kmer_count_p = malloc(number_of_kmers * sizeof(long int) + 5); */
/* 	  background_kmer_count_p = malloc(number_of_kmers * sizeof(long int) + 5); */
/* 	  for (counter = 0; counter < number_of_kmers; counter++) { */
/* 	    signal_kmer_count_p[counter] = results[1][current_kmer_length][0][0][counter]; */
/* 	    background_kmer_count_p[counter] = results[0][current_kmer_length][0][0][counter]; */
/* 	  } */
/* 	} */

/* 	/\\* SORT KMERS ACCORDING TO INCIDENCE *\\/ */
/* 	if (lowpaircount == 1) { */
/* 	  struct sumtable sums; */
/* 	  sumtable_init(&sums, current_kmer_length); */
/* 	  Lowpaircounts(&sums, results[1][current_kmer_length][0][0], 0, current_kmer_length); */
/* 	  sumtable_init(&sums, current_kmer_length); */
/* 	} */
/* 	qsort(background_kmer_count_p, number_of_kmers, sizeof(long int), Numeric_sort_long_ints); */
/* 	qsort(signal_kmer_count_p, number_of_kmers, sizeof(long int), Numeric_sort_long_ints); */
/* 	for (background_kmer_count = 0, signal_kmer_count = 0, counter = number_of_kmers * 0.25; counter < number_of_kmers * 0.75; counter++) { */
/* 	  signal_kmer_count += signal_kmer_count_p[counter]; */
/* 	  background_kmer_count += background_kmer_count_p[counter]; */
/* 	} */
/* 	lambda = (double)(signal_kmer_count * number_of_sequences_analyzed[0]) / (double)(background_kmer_count * number_of_sequences_analyzed[1]); */

/* 	//else lambda =  (number_of_sequences_with_hits[0] * number_of_sequences_analyzed[1]) / (number_of_sequences_with_hits[1] * number_of_sequences_analyzed[0]); */
/*       } */

/*       two_hits_connecting_matrix.one_hit_matches = one_hit_pwm[1].max_counts; */

/*       Count_to_normalized_pwm(&np[current_logo], &all_hits_pwm[0]); */
/*       strcpy(np[current_logo].name, "All hits background "); */
/*       sprintf(tempstring, " : %li", all_hits_pwm[0].max_counts); */
/*       strcat(np[current_logo].name, tempstring); */
/*       current_logo++; */

/*       Count_to_normalized_pwm(&np[current_logo], &all_hits_pwm[1]); */
/*       strcpy(np[current_logo].name, "All hits uncorrected "); */
/*       sprintf(tempstring, " : %li", all_hits_pwm[1].max_counts); */
/*       strcat(np[current_logo].name, tempstring); */
/*       current_logo++; */

/*       /\\* GENERATES EXPONENTIAL PWM *\\/ */
/*       double pseudocount; */
/*       double pseudocount2; */
/*       double swap2; */
/*       double sizefactor; */
/*       sizefactor = (number_of_sequences_analyzed[0] / number_of_sequences_analyzed[1]); */
/*       if (multinomial != 0 && even_background == 0 && complex_background == 1) { */
/* 	pseudocount = all_hits_pwm[1].max_counts * 0.01; */
/* 	pseudocount2 = one_hit_pwm[1].max_counts * 0.01; */
/* 	for (counter = 0; counter < all_hits_pwm[0].width; counter++) */
/* 	  for (counter2 = 0; counter2 < 4; counter2++) { */
/* 	    swap2 = */
/* 		10000 * (((double)all_hits_pwm[1].incidence[counter2][counter] + pseudocount) * sizefactor / */
/* 			 ((double)all_hits_pwm[0].incidence[counter2][counter] + pseudocount * sizefactor) - lambda); */
/* 	    all_hits_exponential[1].incidence[counter2][counter] = (int)swap2; */
/* 	    swap2 = */
/* 		10000 * (((double)one_hit_pwm[1].incidence[counter2][counter] + pseudocount2) * sizefactor / */
/* 			 (((double)one_hit_pwm[0].incidence[counter2][counter]) + pseudocount2 * sizefactor) - lambda); */
/* 	    one_hit_exponential[1].incidence[counter2][counter] = (int)swap2; */
/* 	  } */
/*       } */

/*       /\\* SUBTRACTS BACKGROUND MULTINOMIAL DISTRIBUTION FROM SIGNAL *\\/ */
/*       if (multinomial != 0 && even_background == 0 && complex_background == 1) { */
/* 	for (counter = 0; counter < all_hits_pwm[0].width; counter++) */
/* 	  for (counter2 = 0; counter2 < 4; counter2++) { */
/* 	    swap = (all_hits_pwm[1].incidence[counter2][counter]) * sizefactor - lambda * all_hits_pwm[0].incidence[counter2][counter]; */
/* 	    /\\* printf("\nPosition %i %i: background counts %.3f, signal counts %.3f, corrected value %.3f", counter2, counter,  all_hits_pwm[0].incidence[counter2][counter],  all_hits_pwm[1].incidence[counter2][counter], swap); *\\/ */
/* 	    all_hits_pwm[1].incidence[counter2][counter] = swap; */
/* 	    swap = (one_hit_pwm[1].incidence[counter2][counter]) * sizefactor - lambda * one_hit_pwm[0].incidence[counter2][counter]; */
/* 	    one_hit_pwm[1].incidence[counter2][counter] = swap; */
/* 	  } */
/*       } */
/*       if (same_seed_size == 1) */
/* 	max_seed_size = strlen(searchstring); */
/*       else */
/* 	max_seed_size = strlen(searchstring) + 2; */
/*       if (end_trim_seed == 1) */
/* 	max_seed_size = Nlength - 6; */
/*     } else */
/*       break;			// BREAKS FROM SEED ITERATION LOOP IF FLANK TOOL NOT USED */

/*     /\\* SEED ITERATION *\\/ */
/*     /\\* CHECKS IF SEED HAS CONVERGED *\\/ */

/*     if (auto_seed == 0) { */
/*       /\\* SEED IS ALREADY ESTABLISHED (FROM INPUT OR FROM AUTOSEED) *\\/ */
/*       if (round == 1) { */
/* 	sprintf(tempstring, "\n** user specified seed %s", searchstring); */
/* 	strcat(seed_story, tempstring); */
/*       } */
/*       if (iterate_seed == 0) { */
/* 	sprintf(tempstring, " with no seed iteration"); */
/* 	strcat(seed_story, tempstring); */
/* 	break; */
/*       } */

/*       sprintf(tempstring, "\n** ROUND %i **", round); */
/*       strcat(seed_story, tempstring); */

/*       seed_list[round] = Seed_from_count_PWM(&one_hit_pwm[1], 0.35, max_seed_size); */

/*       if (strcmp(searchstring, seed_list[round]) == 0) { */
/* 	sprintf(tempstring, "\n\t* CONVERGED on %s\n", searchstring); */
/* 	strcat(seed_story, tempstring); */
/* 	iterate_seed = 0; */
/*       } else */
/* 	for (counter = 1; counter < round; counter++) */
/* 	  if (strcmp(seed_list[round], seed_list[counter]) == 0) { */
/* 	    sprintf(tempstring, "\n\t* OSCILLATES: round %li and %i are both %s\n", counter, round, searchstring); */
/* 	    strcat(seed_story, tempstring); */
/* 	    iterate_seed = 0; */
/* 	  } */

/*       if (iterate_seed == 0) */
/* 	break; */

/*       sprintf(tempstring, "\n\t* SEED REFINED input:%s output:%s\n", searchstring, seed_list[round]); */
/*       strcat(seed_story, tempstring); */
/*     } else { */
/*       /\\* AUTOSEED FROM KMER COUNTS *\\/ */
/*       /\\* RANKS KMERS BY LOCAL MAX INCIDENCE *\\/ */
/*       struct kmer_incidence_table *top_kmers; */
/*       top_kmers = malloc(sizeof(struct kmer_incidence_table) * number_of_kmers + 5); */
/*       for (current_kmer = 0; current_kmer < number_of_kmers; current_kmer++) { */
/* 	(top_kmers[current_kmer]).kmer = current_kmer; */
/* 	if (Localmax */
/* 	    (results, 1, shortest_kmer, too_long_kmer, current_kmer_length, 0, 0, current_kmer, count_also_spaced_kmers, */
/* 	     kmer_length_difference_cutoff, tempstring)) */
/* 	  (top_kmers[current_kmer]).incidence = results[1][current_kmer_length][0][0][current_kmer]; */
/* 	else */
/* 	  (top_kmers[current_kmer]).incidence = 0; */
/*       } */
/*       qsort(top_kmers, number_of_kmers, sizeof(struct kmer_incidence_table), Sort_according_to_incidence); */

/*       short int new_seeds_skipped_over = seed_from_local_max_number - 1; */
/*       /\\* ELIMINATES REVERSE COMPLEMENTS TO FIND INDEX TO THE nTH LOCAL MAX SEED *\\/ */
/*       for (current_kmer = 0; current_kmer < number_of_kmers && new_seeds_skipped_over > 0;) { */
/* 	current_kmer++; */
/* 	for (counter = 0; counter < current_kmer; counter++) */
/* 	  if (top_kmers[current_kmer].kmer == Reverse_complement_sequence_value_li(top_kmers[counter].kmer, shortest_kmer)) { */
/* 	    // printf("\npalindrome found current,old: "); Kmerprint(top_kmers[current_kmer].kmer, shortest_kmer); printf(","); Kmerprint(top_kmers[counter].kmer, shortest_kmer);  */
/* 	    break; */
/* 	  } */
/* 	if (counter == current_kmer) */
/* 	  new_seeds_skipped_over--; */
/*       } */
/*       seed_list[round] = Stringfromkmervalue(top_kmers[current_kmer].kmer, shortest_kmer); */
/*       strcpy(searchstring, seed_list[round]);	/\\* UPDATES SEED *\\/ */
/*       sprintf(tempstring, "\n** ROUND 1 **\n\t* seed discovery from %i-mers using AUTOSEED: takes local max ranked %i %s as starting seed\n", */
/* 	      shortest_kmer, seed_from_local_max_number, seed_list[round]); */
/*       strcat(seed_story, tempstring); */
/*     } */

/*     /\\* CLEARS PWMs and other data *\\/ */
/*     current_logo = 0; */
/*     kmer_count = 0; */
/*     count_connecting_matrix_clear(&cm, "all hits connecting matrix", Nlength, 0); */
/*     count_connecting_matrix_clear(&two_hits_connecting_matrix, "two hit connecting matrix", Nlength, 0); */
/*     hit_position_matrix_clear(&hit_position, "hit positions", max_Nlength, 0); */
/*     for (file_number = 0; file_number < number_of_files; file_number++) { */
/*       number_of_sequences_with_no_hits[file_number] = 0; */
/*       number_of_sequences_with_hits[file_number] = 0; */
/*       count_pwm_clear(&all_hits_pwm[file_number], "empty", Nlength * 2, 0); */
/*       count_pwm_clear(&nc[file_number], "empty", Nlength * 2, 0); */
/*       count_pwm_clear(&one_hit_pwm[file_number], "empty", Nlength * 2, 0); */
/*       count_pwm_clear(&one_hit_exponential[file_number], "empty", Nlength * 2, 0); */
/*       count_pwm_clear(&all_hits_exponential[file_number], "empty", Nlength * 2, 0); */
/*       for (counter2 = 0; counter2 < 4; counter2++) */
/* 	for (counter3 = 0; counter3 < Nlength * 2; counter3++) */
/* 	  count_pwm_clear(&two_hits_pwm[file_number][counter2][counter3], "empty", Nlength * 2, 0); */
/*     } */
/*     if (dinucleotide_properties == 1) { */
/*       dinucleotide_properties_matrix_clear(&one_hit_di[0], &di, "dinucleotide1", Nlength * 2, 0, query_sequence_length); */
/*       dinucleotide_properties_matrix_clear(&one_hit_di[1], &di, "dinucleotide2", Nlength * 2, 0, query_sequence_length); */
/*       base_dependency_matrix_clear(&dep[0], "Background_dinucleotide_dependencies", Nlength * 2); */
/*       base_dependency_matrix_clear(&dep[1], "Dinucleotide_dependencies", Nlength * 2); */
/*       base_dependency_matrix_clear(&expected_dinucleotides, "Expected_dinucleotides", Nlength * 2); */
/*     } */

/*     file_number = 0; */
/*     iupac_query = 1; */

/*     strcpy(searchstring, seed_list[round]);	/\\* UPDATES SEED *\\/ */
/*     for (counter = 0; counter < strlen(searchstring); counter++) */
/*       if (searchstring[counter] == '?') { */
/* 	sprintf(tempstring, "\n** SEED ERROR in round %i\n", round); */
/* 	strcat(seed_story, tempstring); */
/* 	printf("%s", seed_story); */
/* 	exit(1); */
/*       } */

/*     /\\* END OF SEED ITERATION LOOP *\\/ */
/*   } */
/*   if (round > max_rounds && iterate_seed == 1) { */
/*     sprintf(tempstring, "\n** FAILED TO CONVERGE IN %i ROUNDS\n", max_rounds); */
/*     strcat(seed_story, tempstring); */
/*     printf("%s", seed_story); */
/*     exit(1); */
/*   } */

/*   if (flank == 1) { */

/*     file_number = 1; */
/*     Count_to_normalized_pwm(&p, &all_hits_pwm[file_number]); */

/*     /\\* PRINTS PWMs *\\/ */
/*     /\\* ALL HITS *\\/ */
/*     short int left_flank_length = 10; */
/*     short int right_flank_length = 12; */
/*     first = Nlength - left_flank_length - 1; */
/*     last = Nlength + right_flank_length + query_sequence_length - 1; */
/*     long int position; */
/*     short int nucleotide; */
/*     fflush(stdout); */
/*     printf("\nAll Hits \tPosition"); */
/*     printf("\nAll Hits               "); */
/*     for (position = first; position < last; position++) { */
/*       printf("\t%li", (position - first - left_flank_length + 1)); */
/*     } */
/*     for (counter2 = 0; counter2 < 4; counter2++) { */
/*       printf("\nAll Hits"); */
/*       printf("\t%c", forward[counter2]); */
/*       for (counter = first; counter < last; counter++) { */
/* 	if (print_frequencies == 0) */
/* 	  printf("\t%.0f", all_hits_pwm[file_number].incidence[counter2][counter]); */
/* 	else */
/* 	  printf("\t%0.3f", p.fraction[counter2][counter]); */
/*       } */
/*     } */

/*     fflush(stdout); */
/*     /\\* ONLY ONE HIT *\\/ */
/*     Count_to_normalized_pwm(&p, &one_hit_pwm[file_number]); */
/*     printf("\n\nOne Hit   \tPosition"); */
/*     printf("\nOne Hit              "); */
/*     for (counter = first; counter < last; counter++) */
/*       printf("\t%li", counter - first - left_flank_length + 1); */
/*     for (counter2 = 0; counter2 < 4; counter2++) { */
/*       printf("\nOne Hit     "); */
/*       printf("\t%c", forward[counter2]); */
/*       for (counter = first; counter < last; counter++) { */
/* 	if (print_frequencies == 0) */
/* 	  printf("\t%.0f", one_hit_pwm[file_number].incidence[counter2][counter]); */
/* 	else */
/* 	  printf("\t%0.3f", p.fraction[counter2][counter]); */
/*       } */
/*     } */

/*     /\\* ONLY ONE HIT, EXPONENTIAL *\\/ */
/*     Count_to_normalized_pwm(&p, &one_hit_exponential[file_number]); */
/*     printf("\n\nOneHit_Exp\tPosition"); */
/*     printf("\nOneHit_Exp          "); */
/*     for (counter = first; counter < last; counter++) */
/*       printf("\t%li", counter - first - left_flank_length + 1); */
/*     for (counter2 = 0; counter2 < 4; counter2++) { */
/*       printf("\nOneHit_Exp "); */
/*       printf("\t%c", forward[counter2]); */
/*       for (counter = first; counter < last; counter++) { */
/* 	if (print_frequencies == 0) */
/* 	  printf("\t%.0f", one_hit_exponential[file_number].incidence[counter2][counter]); */
/* 	else */
/* 	  printf("\t%0.3f", p.fraction[counter2][counter]); */
/*       } */
/*     } */

/*     /\\* ADDS PWMs TO LOGO LIST *\\/ */
/*     Count_to_normalized_pwm(&np[current_logo], &all_hits_pwm[file_number]); */
/*     strcpy(np[current_logo].name, "All hits"); */
/*     sprintf(tempstring, " : %li", all_hits_pwm[file_number].max_counts); */
/*     strcat(np[current_logo].name, tempstring); */
/*     current_logo++; */

/*     Count_to_normalized_pwm(&np[current_logo], &one_hit_pwm[file_number]); */
/*     strcpy(np[current_logo].name, "One hit  "); */
/*     sprintf(tempstring, " : %li", one_hit_pwm[file_number].max_counts); */
/*     strcat(np[current_logo].name, tempstring); */
/*     current_logo++; */

/*     Count_to_normalized_pwm(&np[current_logo], &one_hit_exponential[file_number]); */
/*     strcpy(np[current_logo].name, "One hit exponential"); */
/*     sprintf(tempstring, " : %li", one_hit_pwm[file_number].max_counts);	/\\* ACTUAL MAX COUNTS THAT WERE USED, EXPONENTIAL IS NORMALIZED TO 10000 *\\/ */
/*     strcat(np[current_logo].name, tempstring); */
/*     current_logo++; */

/*     /\\* TWO HITS *\\/ */
/*     for (total_possible_spacings = 0, total_number_of_two_hits = 0, orientation = 0; orientation < 3; orientation++) */
/*       for (spacing = query_sequence_length; spacing < Nlength - query_sequence_length; spacing++) { */
/* 	total_number_of_two_hits += two_hits_pwm[file_number][orientation][spacing].max_counts; */
/* 	total_possible_spacings += (Nlength - query_sequence_length - spacing) * 2; */
/*       } */
/*     printf("\nTotal number of two hits: %li", total_number_of_two_hits); */
/*     printf("\nTotal possible spacings: %li", total_possible_spacings); */
/*     for (orientation = 0; orientation < 3; orientation++) */
/*       for (spacing = 0; spacing < 30; spacing++) { */
/* 	Count_to_normalized_pwm(&p, &two_hits_pwm[file_number][orientation][spacing]); */
/* 	if (two_hits_pwm[file_number][orientation][spacing].max_counts != 0) */
/* 	  current_fold_expected = */
/* 	      (two_hits_pwm[file_number][orientation][spacing].max_counts * total_possible_spacings) / (all_hits_pwm[file_number].max_counts * */
/* 													(Nlength - query_sequence_length - */
/* 													 spacing) * */
/* 													(pow */
/* 													 (all_hits_pwm[file_number].max_counts / */
/* 													  number_of_sequences_analyzed[1], */
/* 													  2) * */
/* 													 exp(-all_hits_pwm[file_number].max_counts / */
/* 													     number_of_sequences_analyzed[1]) / 2)); */
/* 	else */
/* 	  current_fold_expected = 0; */
/* 	two_hits_fold_connecting_matrix.incidence[orientation][spacing] = (long int)(current_fold_expected * 100); */
/* 	/\\* printf("\n%i,%i,%i,%f, %i",orientation, spacing, two_hits_pwm[file_number][orientation][spacing].max_counts, current_fold_expected, two_hits_fold_connecting_matrix.incidence[orientation][spacing]); *\\/ */
/* 	if (two_hits_pwm[file_number][orientation][spacing].max_counts > minimum_kmer_count) { */
/* 	  if (two_hits_pwm[file_number][orientation][spacing].max_counts > max_counts) { */
/* 	    max_counts = two_hits_pwm[file_number][orientation][spacing].max_counts; */
/* 	    max_orientation = orientation; */
/* 	    max_spacing = spacing; */
/* 	  } */
/* 	  Count_to_normalized_pwm(&np[current_logo], &two_hits_pwm[file_number][orientation][spacing]); */
/* 	  strcpy(np[current_logo].name, orientation_string[orientation]); */
/* 	  sprintf(tempstring, "-%i", spacing); */
/* 	  strcat(np[current_logo].name, tempstring); */
/* 	  sprintf(tempstring, " : %li", two_hits_pwm[file_number][orientation][spacing].max_counts); */
/* 	  strcat(np[current_logo].name, tempstring); */
/* 	  sprintf(tempstring, " (%.1f)", current_fold_expected); */
/* 	  strcat(np[current_logo].name, tempstring); */
/* 	  current_logo++; */
/* 	  printf("\n%s-%i", orientation_string[orientation], spacing); */
/* 	  printf("\n%s-%i", orientation_string[orientation], spacing); */
/* 	  printf("\tTwo hits with spacing %i, orientation %s, max_counts %li", spacing, orientation_string[orientation], */
/* 		 two_hits_pwm[file_number][orientation][spacing].max_counts); */
/* 	  printf("\n%s-%i", orientation_string[orientation], spacing); */
/* 	  for (counter = 20; counter <= 43; counter++) */
/* 	    printf("\t%li", counter - 26); */
/* 	  for (counter2 = 0; counter2 < 4; counter2++) { */
/* 	    printf("\n%s-%i", orientation_string[orientation], spacing); */
/* 	    printf("\t%c", forward[counter2]); */
/* 	    for (counter = 20; counter < 43; counter++) { */
/* 	      if (print_frequencies == 0) */
/* 		printf("\t%.0f", two_hits_pwm[file_number][orientation][spacing].incidence[counter2][counter]); */
/* 	      else */
/* 		printf("\t%0.3f", p.fraction[counter2][counter]); */
/* 	    } */
/* 	  } */
/* 	} */
/*       } */

/*     /\\* PRINTS MULTIPLE HIT CONNECTING MATRIX *\\/ */
/*     Count_to_normalized_connecting_matrix(&cp, &cm); */
/*     printf("\n\n! Two or more hits connecting matrix\n! Orientation            \tSpacing\n!                  "); */
/*     for (counter = 0; counter < Nlength - query_sequence_length; counter++) */
/*       printf("\t%li", counter); */
/*     printf */
/* 	("\n! ------------------------------------------------------------------------------------------------------------------------\n! Head to tail > > :"); */
/*     for (counter = 0; counter < Nlength - query_sequence_length; counter++) { */
/*       if (print_frequencies == 0) */
/* 	printf("\t%li", cm.incidence[head_to_tail][counter]); */
/*       else */
/* 	printf("\t%0.3f", cp.fraction[head_to_tail][counter]); */
/*     } */
/*     printf("\n! Head to head > < :"); */
/*     for (counter = 0; counter < Nlength - query_sequence_length; counter++) { */
/*       if (print_frequencies == 0) */
/* 	printf("\t%li", cm.incidence[head_to_head][counter]); */
/*       else */
/* 	printf("\t%0.3f", cp.fraction[head_to_head][counter]); */
/*     } */
/*     printf("\n! Tail to tail < > :"); */
/*     for (counter = 0; counter < Nlength - query_sequence_length; counter++) { */
/*       if (print_frequencies == 0) */
/* 	printf("\t%li", cm.incidence[tail_to_tail][counter]); */
/*       else */
/* 	printf("\t%0.3f", cp.fraction[tail_to_tail][counter]); */
/*     } */

/*     /\\* PRINTS TWO HIT CONNECTING MATRIX *\\/ */
/*     Count_to_normalized_connecting_matrix(&cp, &two_hits_connecting_matrix); */
/*     printf("\n\n2! Exactly two hits connecting matrix\n2! Orientation            \tSpacing\n2!                  "); */
/*     for (counter = 0; counter < Nlength - query_sequence_length; counter++) */
/*       printf("\t%li", counter); */
/*     printf */
/* 	("\n2! ------------------------------------------------------------------------------------------------------------------------\n2! Head to tail > > :"); */
/*     for (counter = 0; counter < Nlength - query_sequence_length; counter++) { */
/*       if (print_frequencies == 0) */
/* 	printf("\t%li", two_hits_connecting_matrix.incidence[head_to_tail][counter]); */
/*       else */
/* 	printf("\t%0.3f", cp.fraction[head_to_tail][counter]); */
/*     } */
/*     printf("\n2! Head to head > < :"); */
/*     for (counter = 0; counter < Nlength - query_sequence_length; counter++) { */
/*       if (print_frequencies == 0) */
/* 	printf("\t%li", two_hits_connecting_matrix.incidence[head_to_head][counter]); */
/*       else */
/* 	printf("\t%0.3f", cp.fraction[head_to_head][counter]); */
/*     } */
/*     printf("\n2! Tail to tail < > :"); */
/*     for (counter = 0; counter < Nlength - query_sequence_length; counter++) { */
/*       if (print_frequencies == 0) */
/* 	printf("\t%li", two_hits_connecting_matrix.incidence[tail_to_tail][counter]); */
/*       else */
/* 	printf("\t%0.3f", cp.fraction[tail_to_tail][counter]); */
/*     } */

/*   } else */
/*     /\\* PRINTS COUNTS OR FOLD CHANGES *\\/ */
/*     for (file_number = 0; file_number < number_of_files; file_number++) { */
/*       if (print_counts == 0 && file_number > 0) */
/* 	break; */

/*       if (output_all_gap_lengths_in_one_line == 0) { */
/* 	/\\* ALL LINES SEPARATELY *\\/ */
/* 	if (center_spaced_kmer_table == 0) */
/* 	  printf */
/* 	      ("\nBackground\tSignal\tSequence_length\tGap_position\tGap length\tKmer\tIUPAC\tredundancy\tBackground count\tSignal count\tLength_normalized_count\tFold_change\tLocal_max\tRepeat\n"); */
/* 	else */
/* 	  printf("\nSequence\tBackground_count\tSignal_count\tNormalized_incidence\tScore\tRepeat\t:SK\n"); */
/* 	for (current_kmer_length = shortest_kmer; current_kmer_length < too_long_kmer; current_kmer_length++) { */
/* 	  current_gap_position = 0; */
/* 	  kmer_length_size = current_kmer_length; */
/* 	  gap_position_size = (Nlength - current_kmer_length - 1) * (count_also_spaced_kmers != 0) + 1; */
/* 	  for (; current_gap_position < kmer_length_size; current_gap_position++) { */
/* 	    for (current_gap_length = 0; current_gap_length < gap_position_size; current_gap_length++) { */
/* 	      if (current_gap_length == 0 && current_gap_position == 0) ; */
/* 	      else { */
/* 		if (count_also_spaced_kmers == 1 && current_gap_position != current_kmer_length / 2 */
/* 		    && current_gap_position != current_kmer_length / 2 + current_kmer_length % 2) */
/* 		  continue; */
/* 		if (current_gap_length == 0 && current_gap_position > 0) */
/* 		  continue; */
/* 		if (current_gap_position == 0 && current_gap_length > 0) */
/* 		  break; */
/* 	      } */

/* 	      for (current_kmer = 0; current_kmer < pow(4, current_kmer_length); current_kmer++) { */
/* 		if (only_palindromes == 1 && Is_this_sequence_dimer(current_kmer, current_kmer_length) == 0) */
/* 		  continue; */
/* 		top_normalized_count = */
/* 		    results[file_number + 1 - */
/* 			    print_counts][current_kmer_length][current_gap_position][current_gap_length][current_kmer] >> ((too_long_kmer - */
/* 															    current_kmer_length - */
/* 															    1) << 1); */
/* 		if (top_normalized_count >= minimum_kmer_count) { */
/* 		  if (center_spaced_kmer_table == 0) { */
/* 		    printf("%s\t", file_name[file_number]); */
/* 		    if (print_counts == 0) */
/* 		      printf("%s\t", file_name[file_number + 1]); */
/* 		    printf("%i\t%i\t%i\t", current_kmer_length, current_gap_position, current_gap_length); */
/* 		  } */
/* 		  for (position = current_kmer_length - 1; position > -1; position--) { */
/* 		    if (current_kmer_length - position - 1 == current_gap_position) */
/* 		      for (counter = 0; counter < current_gap_length; counter++) */
/* 			printf("n"); */
/* 		    printf("%c", forward[(current_kmer & mask_ULL[1][position]) >> (position * 2)]); */
/* 		  } */
/* 		  kmer_incidence = results[file_number][current_kmer_length][current_gap_position][current_gap_length][current_kmer]; */

/* 		  /\\* PRINTS IUPAC *\\/ */
/* 		  printf("\t"); */

/* 		  for (represents_n_kmers = 1, position = current_kmer_length - 1; position > -1; position--) { */
/* 		    if (current_kmer_length - position - 1 == current_gap_position) */
/* 		      for (counter = 0; counter < current_gap_length; counter++) */
/* 			printf("n"); */
/* 		    for (iupac_bits = 0, max_nucleotide = 0, max_kmer_count = 0, nucleotide = 0; nucleotide < 4; nucleotide++) { */
/* 		      test_kmer = */
/* 			  ((current_kmer & (~(mask_ULL[1][position]))) | (nucleotide << (position * 2))) & lowmask_ULL[current_kmer_length - 1]; */
/* // printf("  compared to: "); for(position2 = current_kmer_length-1; position2 > -1 ; position2--) {if(current_kmer_length - position2 - 1  == current_gap_position) for(counter = 0; counter < current_gap_length; counter++) printf("n");printf("%c", forward[(test_kmer & mask_ULL[1][position2]) >> (position2 * 2)]);} */
/* 		      test_kmer_count = results[file_number + 1][current_kmer_length][current_gap_position][current_gap_length][test_kmer]; */
/* 		      if (test_kmer_count > max_kmer_count) { */
/* 			max_nucleotide = nucleotide; */
/* 			max_kmer_count = test_kmer_count; */
/* 			// printf("\tMAX"); */
/* 		      } */
/* 		    } */
/* 		    iupac_bits |= (1 << max_nucleotide); */
/* 		    for (represents_n_nucleotides = 1, nucleotide = 0; nucleotide < 4; nucleotide++) { */
/* 		      if (nucleotide == max_nucleotide) */
/* 			continue; */
/* 		      test_kmer = (current_kmer & (~(mask_ULL[1][position]))) | (nucleotide << (position * 2)) & lowmask_ULL[current_kmer_length - 1]; */
/* 		      test_kmer_count = results[file_number + 1][current_kmer_length][current_gap_position][current_gap_length][test_kmer]; */
/* 		      /\\* IUPACs BASE IF IT IS MORE THAN CUTOFF FRACTION OF MAX *\\/ */
/* 		      if (test_kmer_count > max_kmer_count * iupac_cutoff) { */
/* 			iupac_bits |= (1 << nucleotide); */
/* 			represents_n_nucleotides++; */
/* 		      } */
/* 		    } */
/* 		    printf("%c", nucleotide_bitiupac[iupac_bits]); */
/* 		    represents_n_kmers *= represents_n_nucleotides; */
/* 		  } */
/* 		  printf("\t%li", represents_n_kmers); */
/* 		  if (print_counts == 1) */
/* 		    printf("\t%li\t%.1f\t%s\n", kmer_incidence, */
/* 			   ((double)kmer_incidence) * pow(4, current_kmer_length) / (Nlength - current_kmer_length - current_gap_length), */
/* 			   repeat_report[Is_this_sequence_dimer(current_kmer, current_kmer_length)]); */
/* 		  else { */
/* 		    kmer2_incidence = results[file_number + 1][current_kmer_length][current_gap_position][current_gap_length][current_kmer]; */
/* 		    printf("\t%li\t%li\t%li\t%.1f\t%s", kmer_incidence, kmer2_incidence, top_normalized_count, */
/* 			   (double)(number_of_sequences_analyzed[file_number] / number_of_sequences_analyzed[file_number + 1]) * */
/* 			   (double)kmer2_incidence / (double)kmer_incidence, */
/* 			   repeat_report[Is_this_sequence_dimer(current_kmer, current_kmer_length)]); */
/* 		    if (print_p_values == 1) */
/* 		      printf("\t%f", */
/* 			     Winflat(kmer_incidence, kmer2_incidence, number_of_sequences_analyzed[file_number], */
/* 				     number_of_sequences_analyzed[file_number + 1])); */
/* 		    if (Localmax */
/* 			(results, file_number + 1, shortest_kmer, too_long_kmer, current_kmer_length, current_gap_position, current_gap_length, */
/* 			 current_kmer, count_also_spaced_kmers, kmer_length_difference_cutoff, tempstring)) */
/* 		      printf("\tlocal_max"); */
/* 		    else */
/* 		      printf("\t"); */
/* 		    if (center_spaced_kmer_table == 0) */
/* 		      printf("\n"); */
/* 		    else */
/* 		      printf("\t:SK\n"); */
/* 		  } */
/* 		} */
/* 	      } */
/* 	    } */
/* 	  } */
/* 	} */
/*       } else { */
/* 	/\\* ALL GAP LENGTHS IN ONE LINE *\\/ */
/* 	printf("\n__%s\t\t\t\tGAP LENGTHS\nFile            \tLength\tSequence\tTotal\tMax\t", file_name[file_number]); */
/* 	for (current_gap_length = 0; current_gap_length < Nlength - shortest_kmer; current_gap_length++) */
/* 	  printf("\t%i", current_gap_length); */
/* 	printf("\tRepeat\n"); */
/* 	for (current_kmer_length = shortest_kmer; current_kmer_length < too_long_kmer; current_kmer_length++) { */
/* 	  for (current_gap_position = 0; current_gap_position < current_kmer_length; current_gap_position++) { */
/* 	    if (current_gap_length == 0 && current_gap_position == 0) ; */
/* 	    else if (count_also_spaced_kmers == 1 && current_gap_position != current_kmer_length / 2 */
/* 		     && current_gap_position != current_kmer_length / 2 + current_kmer_length % 2) */
/* 	      continue; */

/* 	    for (current_kmer = 0; current_kmer < pow(4, current_kmer_length); current_kmer++) { */
/* 	      localmaxes = 0; */
/* 	      max_incidence = 0; */
/* 	      total_incidence = 0; */
/* 	      if (only_palindromes == 1 && Is_this_sequence_dimer(current_kmer, current_kmer_length) == 0) */
/* 		continue; */
/* 	      for (current_gap_length = 0; current_gap_length < Nlength - current_kmer_length; current_gap_length++) { */
/* 		/\\* printf("\n%i\t%i\t%i\t%i", current_kmer_length, current_gap_position, current_gap_length, current_kmer); */
/* 		   fflush(stdout); *\\/ */
/* 		if (current_gap_length == 0) */
/* 		  kmer_incidence = results[file_number][current_kmer_length][0][0][current_kmer]; */
/* 		else */
/* 		  kmer_incidence = results[file_number][current_kmer_length][current_gap_position][current_gap_length][current_kmer]; */
/* 		if (print_counts == 0) { */
/* 		  if (current_gap_length == 0) */
/* 		    kmer2_incidence = results[file_number + 1][current_kmer_length][0][0][current_kmer]; */
/* 		  else */
/* 		    kmer2_incidence = results[file_number + 1][current_kmer_length][current_gap_position][current_gap_length][current_kmer]; */
/* 		  total_incidence += kmer2_incidence; */
/* 		  if (((number_of_sequences_analyzed[file_number] / number_of_sequences_analyzed[file_number + 1]) * (double)(kmer2_incidence + 1) / */
/* 		       (double)(kmer_incidence + 1)) > max_incidence) */
/* 		    max_incidence = */
/* 			((number_of_sequences_analyzed[file_number] / number_of_sequences_analyzed[file_number + 1]) * (double)(kmer2_incidence + 1) / */
/* 			 (double)(kmer_incidence + 1)); */
/* 		} else { */
/* 		  total_incidence += kmer_incidence; */
/* 		  if (kmer_incidence > max_incidence) */
/* 		    max_incidence = kmer_incidence; */
/* 		} */
/* 	      } */
/* 	      if (max_incidence >= minimum_kmer_count) { */

/* 		printf("__%s\t", file_name[file_number]); */
/* 		if (print_counts == 0) */
/* 		  printf("%s\t", file_name[file_number + 1]); */
/* 		printf("%i\t", current_kmer_length); */
/* 		for (position = current_kmer_length - 1; position > -1; position--) { */
/* 		  if (current_kmer_length - position - 1 == current_gap_position) */
/* 		    printf("n"); */
/* 		  printf("%c", forward[(current_kmer & mask_ULL[1][position]) >> (position * 2)]); */
/* 		} */
/* 		palindrome_correction = 1 + (Is_this_sequence_dimer(current_kmer, current_kmer_length) == 4); */
/* 		printf("\t%li\t%.0f\t", total_incidence / palindrome_correction, max_incidence / palindrome_correction); */
/* 		for (current_gap_length = 0; current_gap_length < Nlength - current_kmer_length; current_gap_length++) { */
/* 		  if (current_gap_length == 0) */
/* 		    kmer_incidence = results[file_number][current_kmer_length][0][0][current_kmer]; */
/* 		  else */
/* 		    kmer_incidence = results[file_number][current_kmer_length][current_gap_position][current_gap_length][current_kmer]; */
/* 		  if (print_counts == 1) */
/* 		    printf("\t%li", kmer_incidence / palindrome_correction); */
/* 		  else { */
/* 		    if (current_gap_length == 0) */
/* 		      kmer2_incidence = results[file_number + 1][current_kmer_length][0][0][current_kmer]; */
/* 		    else */
/* 		      kmer2_incidence = results[file_number + 1][current_kmer_length][current_gap_position][current_gap_length][current_kmer]; */
/* 		    if (print_p_values == 1) */
/* 		      printf("\t%f", */
/* 			     Winflat(kmer_incidence, kmer2_incidence, number_of_sequences_analyzed[file_number], */
/* 				     number_of_sequences_analyzed[file_number + 1])); */
/* 		    else */
/* 		      printf("\t%f", */
/* 			     (number_of_sequences_analyzed[file_number] / number_of_sequences_analyzed[file_number + 1]) * (double)kmer2_incidence / */
/* 			     (double)kmer_incidence); */
/* 		  } */
/* 		  if (kmer2_incidence >= minimum_kmer_count && print_counts == 0) */
/* 		    if (Localmax */
/* 			(results, file_number + 1, shortest_kmer, too_long_kmer, current_kmer_length, current_gap_position, current_gap_length, */
/* 			 current_kmer, count_also_spaced_kmers, kmer_length_difference_cutoff, tempstring)) */
/* 		      localmaxes++; */
/* 		} */
/* 		printf("\t%s", repeat_report[Is_this_sequence_dimer(current_kmer, current_kmer_length)]); */
/* 		if (localmaxes != 0) */
/* 		  printf("\tlocal_max"); */
/* 		printf("\n"); */
/* 	      } */
/* 	    } */
/* 	  } */
/* 	} */
/*       } */
/*       if (count_also_spaced_kmers != 0) */
/* 	Kmer_svg(user_specified_output_file, results, file_number + 1, shortest_kmer, too_long_kmer, count_also_spaced_kmers, */
/* 		 kmer_length_difference_cutoff, minimum_kmer_count, number_of_sequences_analyzed[file_number], */
/* 		 number_of_sequences_analyzed[file_number + 1], number_of_heatmap_rows, tempstring, match_orientations, match_length); */
/*     } */

/* /\\* PRINTS NUCLEOTIDE COUNTS *\\/ */
/*   if (print_nucleotides == 1) { */
/*     printf("\n\n# Nucleotide counts"); */
/*     for (strand = 0; strand < 2; strand++) { */
/*       if (strand == 0) */
/* 	printf("\n# Forward strand\n#\t"); */
/*       else */
/* 	printf("\n# Reverse strand\n#\t"); */
/*       for (counter = 1; counter < Nlength; counter++) */
/* 	printf("\t%li", counter); */
/*       for (nucleotide_value = 0; nucleotide_value < 4; nucleotide_value++) { */
/* 	printf("\n#\t%c", forward[nucleotide_value]); */
/* 	for (position = 0; position < Nlength - 1; position++) */
/* 	  printf("\t%.0f", nc[strand].incidence[nucleotide_value][position]); */
/*       } */
/*     } */
/*   } */

/*   if (print_nucleotides == 1) { */
/*     printf("\n\n# Nucleotide frequencies"); */
/*     for (strand = 0; strand < 2; strand++) { */
/*       if (strand == 0) */
/* 	printf("\n# Forward strand\n#\t"); */
/*       else */
/* 	printf("\n# Reverse strand\n#\t"); */
/*       for (counter = 1; counter < Nlength; counter++) */
/* 	printf("\t%li", counter); */
/*       for (nucleotide_value = 0; nucleotide_value < 4; nucleotide_value++) { */
/* 	printf("\n#\t%c", forward[nucleotide_value]); */
/* 	for (position = 0; position < Nlength - 1; position++) */
/* 	  printf("\t%.3f", background_pwm[strand].fraction[nucleotide_value][position]); */
/*       } */
/*     } */
/*   } */

/*   fflush(stdout); */

/*   if (flank == 1) { */
/*     strcpy(background_pwm[0].name, "Forward backgr"); */
/*     strcpy(background_pwm[1].name, "Reverse backgr"); */
/*     np_p[current_logo] = &background_pwm[0]; */
/*     current_logo++; */
/*     np_p[current_logo] = &background_pwm[1]; */
/*     current_logo++; */
/*     if (user_specified_output_file[0] == '\0') { */
/*       strcpy(tempstring, file_name[1]); */
/*       strcat(tempstring, "_logo.svg"); */
/*     } else */
/*       strcpy(tempstring, user_specified_output_file); */

/*     cm.number_of_total_matches = all_hits_pwm[1].max_counts; */
/*     two_hits_connecting_matrix.number_of_total_matches = all_hits_pwm[1].max_counts; */
/*     two_hits_fold_connecting_matrix.number_of_total_matches = all_hits_pwm[1].max_counts; */
/*     cms_to_heatmap[1] = &two_hits_connecting_matrix; */
/*     cms_to_heatmap[2] = &two_hits_fold_connecting_matrix; */

/*     /\\* NORMALIZES HIT POSITIONS AND PRINTS THEM *\\/ */
/*     for (max_counts = 0, counter2 = 0; counter2 < 2; counter2++) */
/*       for (counter = 1; counter < Nlength - query_sequence_length; counter++) */
/* 	max_counts += hit_position.incidence[counter2][counter]; */
/*     for (counter2 = 0; counter2 < 2; counter2++) */
/*       for (counter = 1; counter < Nlength - query_sequence_length; counter++) */
/* 	hit_position.fraction[counter2][counter] = */
/* 	    ((Nlength - query_sequence_length) * 2 - 2) * (double)hit_position.incidence[counter2][counter] / (double)max_counts; */

/*     printf("\n\n+\n+ Hit positions"); */
/*     for (counter = 1; counter <= Nlength - query_sequence_length; counter++) */
/*       printf("\t%li", counter); */
/*     printf("\tInformation content\n+ Forward: "); */
/*     for (counter = 1; counter <= Nlength - query_sequence_length; counter++) { */
/*       if (print_frequencies == 1) */
/* 	printf("\t%.1f", hit_position.fraction[0][counter]); */
/*       else */
/* 	printf("\t%li", hit_position.incidence[0][counter]); */
/*     } */
/*     forward_pos_ic = Information_content(hit_position.incidence[0], 1, Nlength - query_sequence_length + 1); */
/*     if (forward_pos_ic > 1) */
/*       warning++; */
/*     if (forward_pos_ic > 2) */
/*       warning++; */
/*     printf("\t%.2f", forward_pos_ic); */
/*     printf("\n+ Reverse: "); */
/*     for (counter = 1; counter <= Nlength - query_sequence_length; counter++) { */
/*       if (print_frequencies == 1) */
/* 	printf("\t%.1f", hit_position.fraction[1][Nlength - query_sequence_length - counter]); */
/*       else */
/* 	printf("\t%li", hit_position.incidence[1][Nlength - query_sequence_length - counter]); */
/*     } */
/*     reverse_pos_ic = Information_content(hit_position.incidence[1], 0, Nlength - query_sequence_length - 1); */
/*     if (reverse_pos_ic > 1) */
/*       warning++; */
/*     if (reverse_pos_ic > 2) */
/*       warning++; */
/*     printf("\t%.2f", reverse_pos_ic); */

/*     /\\* PRINTS HIT STATISTICS AND LAMBDA *\\/ */
/*     printf("\n\n\tTotal sequences\tHits\tFraction\tLower half 8mer total\tMax 8mer count"); */
/*     if (kmer_count == 1) */
/*       printf("\nBackground\t%.0f\t%.0f\t%.0f\t%li\t%li", number_of_sequences_analyzed[0], number_of_sequences_with_hits[0], */
/* 	     100 * number_of_sequences_with_hits[0] / number_of_sequences_analyzed[0], background_kmer_count, */
/* 	     results[0][current_kmer_length][0][0][number_of_kmers - 1]); */
/*     if (kmer_count == 1) */
/*       printf("\nSignal    \t%.0f\t%.0f\t%.0f\t%li\t%li", number_of_sequences_analyzed[1], number_of_sequences_with_hits[1], */
/* 	     100 * number_of_sequences_with_hits[1] / number_of_sequences_analyzed[1], signal_kmer_count, */
/* 	     results[1][current_kmer_length][0][0][number_of_kmers - 1]); */
/*     /\\* printf ("\nLower half 8mer counts: background %i, signal %i", background_kmer_count,  signal_kmer_count); *\\/ */

/*     printf("\nLambda %.3f", lambda); */
/*     /\\* printf("\nshortest kmer = %i and too long kmer = %i bak %.0f sig %.0f", shortest_kmer, too_long_kmer, number_of_sequences_analyzed[0], number_of_sequences_analyzed[1]); *\\/ */

/* /\\* DINUCLEOTIDE PROPERTIES *\\/ */
/*     if (dinucleotide_properties == 1) { */
/*       Expected_dinucleotides(&dep[1], &expected_dinucleotides, searchstring, 1000000000); */

/*       for (counter = end_trim; counter < dep[1].width - end_trim; counter++) { */
/* 	for (counter2 = end_trim; counter2 < dep[1].width - end_trim; counter2++) { */
/* 	  ic_score = dep[1].total_relative_deviation[counter][counter2]; */

/* 	  if (ic_score > max_ic_score) { */
/* 	    for (total_exp = 0, total_obs = 0, counter3 = 0; counter3 < 16; counter3++) { */
/* 	      total_obs += (double)dep[1].incidence[counter][counter2][counter3]; */
/* 	      total_exp += (double)expected_dinucleotides.incidence[counter][counter2][counter3]; */
/* 	    } */
/* 	    if (total_obs > minimum_kmer_count) { */
/* 	      max_ic_score = ic_score; */
/* 	      delta_ic = ic_score; */
/* 	      max_first = counter; */
/* 	      max_second = counter2; */
/* 	      for (min_fold = 1000, max_fold = -1000, counter3 = 0; counter3 < 16; counter3++) { */
/* 		current_fold = */
/* 		    (total_exp * (double)dep[1].incidence[counter][counter2][counter3]) / (total_obs * */
/* 											   (double)expected_dinucleotides. */
/* 											   incidence[counter][counter2][counter3]); */
/* 		if (current_fold < min_fold) { */
/* 		  min_fold = current_fold; */
/* 		  min_dinucleotide = counter3; */
/* 		} */
/* 		if (current_fold > max_fold) { */
/* 		  max_fold = current_fold; */
/* 		  max_dinucleotide = counter3; */
/* 		} */
/* 	      } */
/* 	    } */
/* 	  } */
/* 	} */
/*       } */

/*       (*(np_p[3])).pairwise_correlation[1].first_base = max_first; */
/*       (*(np_p[3])).pairwise_correlation[1].second_base = max_second; */
/*       (*(np_p[3])).pairwise_correlation[1].max_dinucleotide = max_dinucleotide; */
/*       (*(np_p[3])).pairwise_correlation[1].min_dinucleotide = min_dinucleotide; */
/*       (*(np_p[3])).pairwise_correlation[1].max_fold_change = max_fold; */
/*       (*(np_p[3])).pairwise_correlation[1].min_fold_change = min_fold; */
/*       (*(np_p[3])).pairwise_correlation[1].delta_ic = delta_ic; */
/*       strcpy((*(np_p[3])).seed, searchstring); */
/*       Svg_logo(tempstring, current_logo, np_p, cms_to_heatmap, &one_hit_di[1], &dep[0], &dep[1], &expected_dinucleotides, all_hits_align_scores, */
/* 	       lambda, warning); */
/*     } else { */
/*       Svg_logo(tempstring, current_logo, np_p, cms_to_heatmap, '\0', '\0', '\0', '\0', all_hits_align_scores, lambda, warning); */
/*     } */

/*   } */

/*   /\\* OUTPUTS KMER TABLE *\\/ */
/*   if (kmer_table == 1 && flank != 1) { */
/*     current_gap_position = current_kmer_length / 2; */
/*     current_gap_length = 0; */
/*     last_kmer[0] = pow(4, (too_long_kmer - 1)); */
/*     for (file_number = 0; file_number < number_of_files; file_number++) { */

/*       for (current_gap_length = 0; current_gap_length < loaded_pwm_width - current_kmer_length + 2; current_gap_length++) { */
/* 	for (current_kmer_length = shortest_kmer; current_kmer_length < too_long_kmer; current_kmer_length++) { */
/* 	  last_kmer[current_kmer_length] = pow(4, current_kmer_length); */
/* 	  if (file_number == 1) */
/* 	    printf("\nKMER\tsequence"); */
/* 	  else */
/* 	    printf("\nBACKGROUND KMER\tsequence"); */
/* 	  printf("\t%imer\t%imer_expected", current_kmer_length, current_kmer_length); */
/* 	} */
/* 	if (file_number == 1) */
/* 	  printf("\t::: %i:", current_gap_length); */
/* 	else */
/* 	  printf("\t::B %i:", current_gap_length); */
/* 	for (current_kmer = 0; current_kmer < last_kmer[0]; current_kmer++) { */
/* 	  printf("\n%li\t", current_kmer); */
/* 	  for (position = too_long_kmer - 2; position >= shortest_kmer & current_kmer <= last_kmer[position]; position--) { */
/* 	    printf("%c", forward_lc[(current_kmer & mask_ULL[1][position]) >> (position * 2)]); */
/* 	    if (position == current_kmer_length / 2) */
/* 	      for (counter = 0; counter < current_gap_length; counter++) */
/* 		printf("n"); */
/* 	  } */
/* 	  for (; position > -1; position--) { */
/* 	    printf("%c", forward[(current_kmer & mask_ULL[1][position]) >> (position * 2)]); */
/* 	    if (position == current_kmer_length / 2) */
/* 	      for (counter = 0; counter < current_gap_length; counter++) */
/* 		printf("n"); */
/* 	  } */
/* 	  for (current_kmer_length = shortest_kmer; current_kmer_length < too_long_kmer; current_kmer_length++) { */
/* 	    if (current_kmer < last_kmer[current_kmer_length]) { */
/* 	      /\\* TESTS IF CURRENT SEQUENCE IS A PALINDROME *\\/ */
/* 	      palindrome_correction = 1 + (Is_this_sequence_dimer(current_kmer, current_kmer_length) == 4); */
/* 	      if (current_gap_length == 0) */
/* 		kmer_incidence = results[file_number][current_kmer_length][0][0][current_kmer]; */
/* 	      else */
/* 		kmer_incidence = results[file_number][current_kmer_length][current_gap_position][current_gap_length][current_kmer]; */
/* 	      printf("\t%li", kmer_incidence / palindrome_correction); */
/* 	      printf("\t%.2f", Kmerscore(&qp, current_kmer, current_kmer_length, kmermatch_position)); */
/* 	      printf("\t%.2f", */
/* 		     fastgappedKmerscore(qp.fraction, qp.width, current_kmer & lowmask_ULL[current_kmer_length / 2 - 1], */
/* 					 (current_kmer >> current_kmer_length) & lowmask_ULL[current_kmer_length / 2 - 1], current_kmer_length / 2, */
/* 					 current_gap_length)); */
/* 	      if (load_adm == 1) */
/* 		printf("\t%.2f", Kmerscore_ADM(&flanked_adm, current_kmer, current_kmer_length, kmermatch_position)); */
/* 	    } else */
/* 	      printf("\t-\t-"); */
/* 	  } */
/* 	  if (file_number == 1) */
/* 	    printf("\t::: %i:", current_gap_length); */
/* 	  else */
/* 	    printf("\t::B %i:", current_gap_length); */
/* 	} */
/*       } */
/*     } */
/*   } */

/*   /\\* COUNTS INFORMATION CONTENT OF HIT AND NON-HIT K-MER DISTRIBUTIONS *\\/ */
/*   current_gap_position = 0; */
/*   current_gap_length = 0; */

/*   for (current_kmer_length = shortest_kmer; current_kmer_length < too_long_kmer; current_kmer_length++) { */
/*     background_info = 2 * current_kmer_length; */
/*     signal_info = 2 * current_kmer_length; */
/*     if (flank == 1) { */
/*       background_nonhit_info = 2 * current_kmer_length; */
/*       signal_nonhit_info = 2 * current_kmer_length; */
/*     } */
/*     number_of_kmers = pow(4, current_kmer_length); */
/*     for (counter = 0; counter < number_of_kmers; counter++) { */
/*       fraction = results[1][current_kmer_length][0][0][counter] / (number_of_sequences_analyzed[1] * 2 * (Nlength - current_kmer_length)); */
/*       if (results[1][current_kmer_length][0][0][counter] != 0) */
/* 	signal_info += fraction * log2(fraction); */
/*       fraction = results[0][current_kmer_length][0][0][counter] / (number_of_sequences_analyzed[0] * 2 * (Nlength - current_kmer_length)); */
/*       if (results[0][current_kmer_length][0][0][counter] != 0) */
/* 	background_info += fraction * log2(fraction); */
/*       if (flank == 1) { */
/* 	fraction = results[1][current_kmer_length][1][0][counter] / (number_of_sequences_with_no_hits[1] * 2 * (Nlength - current_kmer_length)); */
/* 	if (results[1][current_kmer_length][1][0][counter] != 0) */
/* 	  signal_nonhit_info += fraction * log2(fraction); */
/* 	fraction = results[0][current_kmer_length][1][0][counter] / (number_of_sequences_with_no_hits[0] * 2 * (Nlength - current_kmer_length)); */
/* 	if (results[0][current_kmer_length][1][0][counter] != 0) */
/* 	  background_nonhit_info += fraction * log2(fraction); */
/*       } */
/*     } */
/*     printf */
/* 	("\n\nInformation content for kmer length %i\tAll sequences\tSequences not hit\nSignal          \t%.2f\t%.2f\tbits\nBackground      \t%.2f\t%.2f\tbits", */
/* 	 current_kmer_length, signal_info, signal_nonhit_info, background_info, background_nonhit_info); */
/*   } */

/* /\\* CALCULATES AND PRINTS ALL INFORMATION CONTENTS *\\/ */
/*   if (information_content_output == 1 && flank != 1) { */
/*     double signal_information_content[too_long_kmer][too_long_kmer][too_long_kmer]; */
/*     double background_information_content[too_long_kmer][too_long_kmer][too_long_kmer]; */
/*     long int signal_count_sum[too_long_kmer][too_long_kmer][too_long_kmer]; */
/*     long int background_count_sum[too_long_kmer][too_long_kmer][too_long_kmer]; */
/*     printf("\n\nINFORMATION CONTENT\tKmer\tGap_position\tGap_length\tSignal_bits\tBackground_bits\tSignal_counts\tBackground_counts\tIC"); */

/*     for (current_kmer_length = shortest_kmer; current_kmer_length < too_long_kmer; current_kmer_length++) { */
/*       for (current_gap_position = 0; current_gap_position < current_kmer_length; current_gap_position++) { */
/* 	for (current_gap_length = (current_gap_position > 0); current_gap_length < Nlength - current_kmer_length; current_gap_length++) { */
/* 	  signal_information_content[current_kmer_length][current_gap_position][current_gap_length] = 2 * current_kmer_length; */
/* 	  background_information_content[current_kmer_length][current_gap_position][current_gap_length] = 2 * current_kmer_length; */
/* 	  signal_count_sum[current_kmer_length][current_gap_position][current_gap_length] = 0; */
/* 	  background_count_sum[current_kmer_length][current_gap_position][current_gap_length] = 0; */
/* 	  if (current_gap_position == 0 & current_gap_length != 0) */
/* 	    continue; */

/* 	  for (current_kmer = 0; current_kmer < pow(4, current_kmer_length); current_kmer++) { */
/* 	    signal_count_sum[current_kmer_length][current_gap_position][current_gap_length] += */
/* 		results[1][current_kmer_length][current_gap_position][current_gap_length][current_kmer]; */
/* 	    background_count_sum[current_kmer_length][current_gap_position][current_gap_length] += */
/* 		results[0][current_kmer_length][current_gap_position][current_gap_length][current_kmer]; */
/* 	    fraction = */
/* 		results[1][current_kmer_length][current_gap_position][current_gap_length][current_kmer] / (number_of_sequences_analyzed[1] * 2 * */
/* 													   (Nlength - current_kmer_length - */
/* 													    current_gap_length)); */
/* 	    /\\* printf("\nkmer length %i,position %i,gap length %i, current_kmer %i, fraction %f", current_kmer_length, current_gap_position, current_gap_length, current_kmer, fraction); *\\/ */
/* 	    if (results[1][current_kmer_length][current_gap_position][current_gap_length][current_kmer] != 0) */
/* 	      signal_information_content[current_kmer_length][current_gap_position][current_gap_length] += fraction * log2(fraction); */
/* 	    fraction = */
/* 		results[0][current_kmer_length][current_gap_position][current_gap_length][current_kmer] / (number_of_sequences_analyzed[0] * 2 * */
/* 													   (Nlength - current_kmer_length - */
/* 													    current_gap_length)); */
/* 	    if (results[0][current_kmer_length][current_gap_position][current_gap_length][current_kmer] != 0) */
/* 	      background_information_content[current_kmer_length][current_gap_position][current_gap_length] += fraction * log2(fraction); */
/* 	  } */
/* 	  printf("\n\t%i\t%i\t%i\t%.3f\t%.3f\t%li\t%li\tIC", current_kmer_length, current_gap_position, current_gap_length, */
/* 		 signal_information_content[current_kmer_length][current_gap_position][current_gap_length], */
/* 		 background_information_content[current_kmer_length][current_gap_position][current_gap_length], */
/* 		 signal_count_sum[current_kmer_length][current_gap_position][current_gap_length], */
/* 		 background_count_sum[current_kmer_length][current_gap_position][current_gap_length]); */

/* 	} */
/*       } */
/*     } */
/*   } */

/* /\\* PRINTS DINUCLEOTIDE MATRIX *\\/ */
/*   if (dinucleotide_properties == 1) { */
/*     for (counter = 0; counter < (one_hit_di[1]).number_of_dinucleotide_properties; counter++) { */
/*       printf("\n%s", (one_hit_di[1]).dinucleotide_property_string[counter]); */
/*       for (counter2 = 0; counter2 < (one_hit_di[1]).width; counter2++) */
/* 	printf("\t%.2f", (one_hit_di[1]).score[counter][counter2]); */
/*     } */

/* /\\* */
/* for (counter = 0; counter < np[2].width; counter++) */
/* { */
/* printf("\n%i", counter); */
/* for (counter2 = 0; counter2 < np[2].width; counter2++) printf("\t%.2f", dep.information_content[counter][counter2] - expected_dinucleotides.information_content[counter][counter2] ); */
/* } */
/* *\\/ */
/*     char *in_out_string = "oi"; */
/*     short int first_is_in = 0; */
/*     short int second_is_in = 0; */
/*     short int first_pos_relative_to_seed = 0; */
/*     short int second_pos_relative_to_seed = 0; */
/*     long int max_obs; */
/*     long int current_total_dinuc_count; */
/*     long int current_total_dinuc_background_count; */
/*     double difference; */

/* /\\* PRINTS DINUCLEOTIDE DATA *\\/ */
/*     printf */
/* 	("\ndinucleotide\tpos1\tpos2\tmononuc1_fraction\tmononuc2_fraction\tdinuc_fraction\tdinuc_background_fraction\tcount\tbackground_count\texpected(mono)\texpected(h1)\tobsIC\texpIC\tdeltaIC\tscore\tfold_dev\ttotal_counts\tmax_count\ttotal_rel_dev\texp_min_deviation\teo_correlation\tpermutated_correlation\tic\tin/out\tCpG"); */
/*     for (counter = 0; counter < dep[1].width; counter++) { */
/*       for (counter2 = 0; counter2 < dep[1].width; counter2++) { */
/* /\\* CALCULATES SUMMARY STATISTICS *\\/ */
/* 	for ( /\\*total_relative_deviation[counter][counter2] = 0, uncentered_correlation[counter][counter2] = 0, sum_x_squared = 0, sum_y_squared = 0, sum_xy = 0, *\\/ max_obs = 0, total_exp = 0, total_obs = 0, counter3 = 0; counter3 < 16; counter3++) { */
/* 	  total_obs += (double)dep[1].incidence[counter][counter2][counter3]; */
/* 	  total_exp += (double)expected_dinucleotides.incidence[counter][counter2][counter3]; */
/* 	  if (dep[1].incidence[counter][counter2][counter3] > max_obs) */
/* 	    max_obs = dep[1].incidence[counter][counter2][counter3]; */
/* 	  /\\*   printf("\n%f\t%f", total_obs, total_exp); *\\/ */
/* /\\* CALCULATES TOTAL RELATIVE DEVIATION */
/* total_relative_deviation[counter][counter2] += abs(total_obs * ((double) expected_dinucleotides.incidence[counter][counter2][counter3] / total_exp) - dep[1].incidence[counter][counter2][counter3])/total_obs; *\\/ */
/* 	} */
/* /\\*for (counter3 = 0; counter3 < 16; counter3++) total_relative_deviation[counter][counter2] += abs((double) expected_dinucleotides.incidence[counter][counter2][counter3] - (double) dep[1].incidence[counter][counter2][counter3]) / total_obs;*\\/ */

/* /\\* PRINTS LARGE TABLE *\\/ */
/* 	for (current_total_dinuc_count = 0, current_total_dinuc_background_count = 0, counter3 = 0; counter3 < 16; counter3++) { */
/* 	  current_total_dinuc_count += dep[1].incidence[counter][counter2][counter3]; */
/* 	  current_total_dinuc_background_count += dep[0].incidence[counter][counter2][counter3]; */
/* 	} */
/* 	for (counter3 = 0; counter3 < 16; counter3++) { */
/* 	  first_pos_relative_to_seed = counter2 - Nlength + 2; */
/* 	  second_pos_relative_to_seed = counter - Nlength + 2; */

/* /\\* DETERMINES IF POSITIONS ARE INSIDE SEED *\\/ */
/* 	  first_is_in = 0; */
/* 	  second_is_in = 0; */
/* 	  if (first_pos_relative_to_seed > 0 && first_pos_relative_to_seed <= query_sequence_length) */
/* 	    first_is_in = 1; */
/* 	  if (second_pos_relative_to_seed > 0 && second_pos_relative_to_seed <= query_sequence_length) */
/* 	    second_is_in = 1; */
/* 	  if (counter2 < counter) { */
/* 	    printf */
/* 		("\n%c%c\t%i\t%i\t_%i:%i_\t%.2f\t%.2f\t%.2f\t%.2f\t%li\t%li\t%.1f\t%.1f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.0f\t%li\t%.3f\t%.3f\t%.3f\t%.3f\tic\t%c%c\t", */
/* 		 forward[(counter3 & 12) >> 2], forward[counter3 & 3], first_pos_relative_to_seed, second_pos_relative_to_seed, */
/* 		 first_pos_relative_to_seed, second_pos_relative_to_seed, np[2].fraction[(counter3 & 12) >> 2][counter2], */
/* 		 np[2].fraction[(counter3 & 3)][counter], ((double)dep[1].incidence[counter][counter2][counter3] / current_total_dinuc_count), */
/* 		 ((double)dep[0].incidence[counter][counter2][counter3] / current_total_dinuc_background_count), */
/* 		 dep[1].incidence[counter][counter2][counter3], dep[0].incidence[counter][counter2][counter3], */
/* 		 total_obs * ((double)expected_dinucleotides.incidence[counter][counter2][counter3] / total_exp), */
/* 		 Expected_mismatched_dinucleotides(counter3, dep[1].incidence[counter][counter2]), dep[1].information_content[counter][counter2], */
/* 		 expected_dinucleotides.information_content[counter][counter2], */
/* 		 dep[1].information_content[counter][counter2] - expected_dinucleotides.information_content[counter][counter2], */
/* 		 dep[1].information_content[counter][counter2] * (dep[1].information_content[counter][counter2] - */
/* 								  expected_dinucleotides.information_content[counter][counter2]), */
/* 		 (total_exp * (double)dep[1].incidence[counter][counter2][counter3]) / (total_obs * */
/* 											(double)expected_dinucleotides. */
/* 											incidence[counter][counter2][counter3]), total_obs, max_obs, */
/* 		 dep[1].total_relative_deviation[counter][counter2], dep[1].count_statistic_expected_total_relative_deviation[counter][counter2], */
/* 		 dep[1].eo_correlation[counter][counter2], dep[1].permutated_correlation[counter][counter2], in_out_string[first_is_in], */
/* 		 in_out_string[second_is_in]); */

/* 	    difference = */
/* 		(double)dep[1].incidence[counter][counter2][counter3] / current_total_dinuc_count - */
/* 		(double)dep[0].incidence[counter][counter2][counter3] / current_total_dinuc_background_count; */
/* 	    if (difference < -0.2) */
/* 	      printf("low"); */
/* 	    if (difference > 0.2) */
/* 	      printf("high"); */

/* 	  } */
/* 	  if ((counter3 == 6) && counter - counter2 == 1) */
/* 	    printf("\tCpG"); */

/* 	} */

/*       } */
/*     } */

/*     printf */
/* 	("\n\nSum of relative deviations between expected and observed dinucleotide counts\npos1\tpos2\tdeviation\texp_min_deviation\tuncentered correlation\texcess information\tscore\trdv"); */
/*     for (counter2 = Nlength - 1; counter2 - Nlength <= query_sequence_length; counter2++) { */
/*       for (counter = counter2 + 1; counter - Nlength + 2 <= query_sequence_length; counter++) { */
/* 	first_pos_relative_to_seed = counter2 - Nlength + 2; */
/* 	second_pos_relative_to_seed = counter - Nlength + 2; */
/* 	printf("\n%i\t%i\t%.3f\t%.3f\t%.3f\t%.3f\t%.2f\t%.2f\trdv", first_pos_relative_to_seed, second_pos_relative_to_seed, */
/* 	       dep[1].total_relative_deviation[counter][counter2], dep[1].count_statistic_expected_total_relative_deviation[counter][counter2], */
/* 	       dep[1].eo_correlation[counter][counter2], dep[1].permutated_correlation[counter][counter2], */
/* 	       dep[1].information_content[counter][counter2] - expected_dinucleotides.information_content[counter][counter2], */
/* 	       dep[1].information_content[counter][counter2] * (dep[1].information_content[counter][counter2] - */
/* 								expected_dinucleotides.information_content[counter][counter2])); */
/* 	/\\* printf("\t%f",dep[1].total_relative_deviation[counter][counter2]); *\\/ */
/*       } */
/*     } */

/*     struct adjacent_dinucleotide_model adm; */
/*     adjacent_dinucleotide_model_init(&adm, "adm", query_sequence_length); */
/*     //Generate_ADM(&adm, &dep[1], query_sequence_length); */
/*     Generate_background_corrected_ADM(&adm, &dep[0], &dep[1], query_sequence_length, lambda); */
/*     Print_ADM(&adm, "bc_adm"); */
/*     strcpy(tempstring, file_name[1]); */
/*     strcat(tempstring, "_ADM"); */
/*     strcpy(adm.name, tempstring); */
/*     strcat(tempstring, "_logo.svg"); */
/*     Svg_logo_ADM(tempstring, 0, 0, &adm, 0.05, 0.1, 1000, 0.1); */
/*     strcpy(tempstring, adm.name); */
/*     strcat(tempstring, "_background_corrected_riverlake_logo.svg"); */
/*     Svg_riverlake_logo(tempstring, 0, 0, &adm, 0.05, 0.1, 1000, 0.1); */
/*     strcpy(tempstring, adm.name); */
/*     strcat(tempstring, "_riverlake_logo.svg"); */
/*     Generate_ADM(&adm, &dep[1], query_sequence_length); */
/*     Print_ADM(&adm, "ADM"); */
/*     Svg_riverlake_logo(tempstring, 0, 0, &adm, 0.05, 0.1, 1000, 0.1); */
/*   } */

/* /\\* */
/* for (counter = 0; counter < (one_hit_di[1]).number_of_dinucleotide_properties; counter++)  */
/* { */
/* printf ("\n%s", (one_hit_di[1]).dinucleotide_property_string[counter]); */
/* for (counter2 = 0; counter2 < (one_hit_di[1]).width; counter2++) printf ("\t%i", (one_hit_di[1]).count[counter][counter2]); */
/* } */
/* *\\/ */
/*   double bak_correl; */
/*   double sig_correl; */
/*   /\\* PRINTS FLANK KMER DATA *\\/ */
/*   if (flank_kmer_pos != -100) { */
/*     bak_correl = Kmer_mono_correlation(flank_kmer_count[0], flank_kmer_expected_count[0], 4); */
/*     sig_correl = Kmer_mono_correlation(flank_kmer_count[1], flank_kmer_expected_count[1], 4); */
/*     printf("\n\nFLANK KMERS AT POSITION %i", flank_kmer_pos); */
/*     printf("\nKmer\tback_obs\tback_exp\tsignal_obs\tsignal_exp"); */
/*     for (counter = 0; counter < 256; counter++) */
/*       printf("\n%c%c%c%c\t%li\t%.1f\t%li\t%.1f\tfK", forward[counter >> 6], forward[(counter >> 4) & 3], forward[(counter >> 2) & 3], */
/* 	     forward[counter & 3], flank_kmer_count[0][counter], flank_kmer_expected_count[0][counter], flank_kmer_count[1][counter], */
/* 	     flank_kmer_expected_count[1][counter]); */
/*     printf("\nCorrelation with mononucleotide-based expected kmers:\t%.3f\t%.3f", bak_correl, sig_correl); */
/*   } */

/*   short int gap; */
/*   /\\* PRINTS ALIGN SCORES *\\/ */
/*   if (align_matches == 1) { */
/*     /\\* for(counter = 0; counter < Nlength; counter++) for (counter2 = 0; counter2 < Nlength; counter2++) */
/*        { */
/*        printf("\n%i,%i\t", counter, counter2); */
/*        Seqprint(mask_ULL[counter][counter2], Nlength-1);     */
/*        } *\\/ */
/*     printf("\n\nALIGN SCORES ::a\n"); */
/*     for (counter = 0; counter < Nlength * 2 - 2; counter++) */
/*       printf("%li\t", counter - Nlength + 1); */
/*     printf("::a\n"); */
/*     for (strand = 0; strand < 2; strand++) { */
/*       for (counter = 0; counter < Nlength * 2 - 2; counter++) */
/* 	printf("\t %.1f", */
/* 	       ((double)all_hits_align_scores[1].score[counter][strand]) / ((double)all_hits_align_scores[1].count[counter][strand] + 0.01)); */
/*       printf("\t::a\n"); */
/*     } */
/*     printf("\n\nDirect Repeat abundance"); */
/*     for (counter = 0; counter < Nlength; counter++) { */
/*       printf("\n%li", counter); */
/*       for (gap = 0; gap < Nlength; gap++) */
/* 	printf("\t%li", all_hits_align_scores[1].direct_repeat[counter][gap]); */
/*     } */
/*     printf("\n\nInverted Repeat abundance"); */
/*     for (counter = 0; counter < Nlength; counter++) { */
/*       printf("\n%li", counter); */
/*       for (gap = 0; gap < Nlength; gap++) */
/* 	printf("\t%li", all_hits_align_scores[1].inverted_repeat[counter][gap]); */
/*     } */

/*     printf("\n\nDirect repeat length histogram\tSignal\tBackground\tRatio"); */
/*     for (counter = 0; counter < Nlength; counter++) { */
/*       printf("\n%li\t%li\t%li\t%.1f", counter, align_score_histogram[1][0][counter], align_score_histogram[0][0][counter], */
/* 	     (double)align_score_histogram[1][0][counter] / ((double)align_score_histogram[0][0][counter] + 0.001)); */
/*     } */

/*     printf("\n\nInverted repeat length histogram\tSignal\tBackground\tRatio"); */
/*     for (counter = 0; counter < Nlength; counter++) { */
/*       printf("\n%li\t%li\t%li\t%.1f", counter, align_score_histogram[1][1][counter], align_score_histogram[0][1][counter], */
/* 	     (double)align_score_histogram[1][1][counter] / ((double)align_score_histogram[0][1][counter] + 0.001)); */
/*     } */

/*   } */

/* /\\* MAKES KMER COUNT XY PLOT *\\/ */
/*   if ((xyplot == 1 && kmer_count == 1) || expected_observed_plot == 1) { */
/*     if (user_specified_output_file[0] == '\0') { */
/*       strcpy(plotfilename, file_name[0]); */
/*       strcat(plotfilename, "_"); */
/*       strcat(plotfilename, file_name[1]); */
/*     } else */
/*       strcpy(plotfilename, user_specified_output_file); */
/*     strcat(plotfilename, "_"); */
/*     sprintf(valuestring, "%i", shortest_kmer); */
/*     strcat(plotfilename, valuestring); */
/*     if (expected_observed_plot == 0) { */
/*       strcat(plotfilename, "mer_xyplot.svg"); */
/*       XYplot(plotfilename, results, shortest_kmer, count_also_spaced_kmers, kmer_length_difference_cutoff, file_name[0], file_name[1], tempstring); */
/*     } else { */

/* /\\* MAKES EXPECTED-OBSERVED XYPLOT *\\/ */
/*       strcat(plotfilename, "mer_eo_plot.svg"); */
/* /\\* GENERATES EXPECTED COUNTS *\\/ */
/*       number_of_kmers = pow(4, shortest_kmer); */

/*       long long int total_expected_kmer_count; */
/*       long long int total_observed_kmer_count; */
/*       long long int current_count; */
/*       double size_factor; */

/*       struct normalized_pwm *eopwm; */
/*       struct normalized_pwm pwm_generated_from_input_sequences; */
/*       normalized_pwm_init(&pwm_generated_from_input_sequences, "new", Nlength * 2, 0); */
/*       Count_to_normalized_pwm(&pwm_generated_from_input_sequences, &one_hit_pwm[1]); */
/*       Log_ratio_pwm(&pwm_generated_from_input_sequences); */
/*       if (kmer_table == 1) */
/* 	eopwm = &qp; */
/*       else */
/* 	eopwm = &pwm_generated_from_input_sequences; */

/*       long int max_observed_count = 0; */
/*       long int max_expected_count = 0; */

/*       for (total_expected_kmer_count = 0, total_observed_kmer_count = 0, current_kmer = 0; current_kmer < number_of_kmers; current_kmer++) */
/* 	max_observed_count = MAX(max_observed_count, results[1][shortest_kmer][0][0][current_kmer]); */

/*       for (total_expected_kmer_count = 0, total_observed_kmer_count = 0, current_kmer = 0; current_kmer < number_of_kmers; current_kmer++) { */
/* 	if (load_adm == 1) */
/* 	  current_count = (long long int)pow(2, Kmerscore_ADM(&flanked_adm, current_kmer, shortest_kmer, kmermatch_position)); */
/* 	else */
/* 	  current_count = (long long int)pow(2, Kmerscore(eopwm, current_kmer, shortest_kmer, kmermatch_position)); */
/* 	results[0][shortest_kmer][0][0][current_kmer] = current_count; */
/* 	max_expected_count = MAX(max_expected_count, current_count); */
/* 	total_expected_kmer_count += current_count; */
/* 	total_observed_kmer_count += results[1][shortest_kmer][0][0][current_kmer]; */
/*       } */

/*       size_factor = (double)max_observed_count / max_expected_count; */
/*       printf("\nTotal observed %lli Expected %lli Size factor %.6f", total_observed_kmer_count, total_expected_kmer_count, size_factor); */
/*       for (current_kmer = 0; current_kmer < number_of_kmers; current_kmer++) */
/* 	results[0][shortest_kmer][0][0][current_kmer] *= size_factor; */
/*       XYplot(plotfilename, results, shortest_kmer, count_also_spaced_kmers, kmer_length_difference_cutoff, "Expected", "Observed", tempstring); */
/*     } */
/*   } */

/*   if (flank == 1) { */
/*     printf("\n\nPalindromic hits\nBackground:\t%li\nSignal    :\t%li", palindromic_hits[0], palindromic_hits[1]); */
/*     printf("\n"); */
/*     /\\* PRINTS SEEDS OF MONOMER AND DIMERS THAT ARE MORE THAN 5% OF MATCHES *\\/ */
/*     printf("\n\nType\t%%_of_all\t%%_of_1or2\tOrSeed\tAndSeed\tPWMSeed"); */
/*     /\\* ONE HIT *\\/ */
/*     if (remember_iterate_seed == 0) { */
/*       strcat(seed_story, "\n\n** PWM SEED CALCULATION"); */
/*       strcpy(tempstring, Seed_from_count_PWM(&one_hit_pwm[1], 0.35, max_seed_size)); */
/*     } else */
/*       strcpy(tempstring, seed_list[round]); */
/*     printf("\nOneHit\t%.1f%%\t%.1f%%\t%s\t%s\t%s\tonehit_seed", */
/* 	   100 * (double)two_hits_connecting_matrix.one_hit_matches / two_hits_connecting_matrix.number_of_total_matches, */
/* 	   100 * (double)two_hits_connecting_matrix.one_hit_matches / (two_hits_connecting_matrix.one_hit_matches + */
/* 								       two_hits_connecting_matrix.two_hit_matches), Seed_from_input_seed(searchstring, */
/* 																	 0, 0, '|'), */
/* 	   Seed_from_input_seed(searchstring, 0, 0, '&'), tempstring); */
/*     /\\* TWO HITS *\\/ */
/*     for (orientation = 0; orientation < 3; orientation++) */
/*       for (spacing = 1; spacing < Nlength - query_sequence_length; spacing++) { */
/* 	if ((double)two_hits_connecting_matrix.incidence[orientation][spacing] / */
/* 	    (two_hits_connecting_matrix.one_hit_matches + two_hits_connecting_matrix.two_hit_matches) > 0.05) { */
/* 	  sprintf(tempstring, "\n** %s-%i END TRIMMED PWM SEED", orientation_string[orientation], spacing); */
/* 	  strcat(seed_story, tempstring); */
/* 	  printf("\n%s-%i\t%.1f%%\t%.1f%%\t%s\t%s\t%s\ttwohit_seed", orientation_string[orientation], spacing, */
/* 		 100 * (double)two_hits_connecting_matrix.incidence[orientation][spacing] / two_hits_connecting_matrix.number_of_total_matches, */
/* 		 100 * (double)two_hits_connecting_matrix.incidence[orientation][spacing] / (two_hits_connecting_matrix.one_hit_matches + */
/* 											     two_hits_connecting_matrix.two_hit_matches), */
/* 		 Seed_from_input_seed(searchstring, orientation, spacing, '|'), Seed_from_input_seed(searchstring, orientation, spacing, '&'), */
/* 		 Seed_from_count_PWM(&two_hits_pwm[1][orientation][spacing], 0.35, strlen(searchstring) * 2 + spacing + 2)); */
/* 	  if (spacing < query_sequence_length) */
/* 	    printf("_overlap"); */
/* 	} */
/*       } */
/*   } */

/*   printf("%s", seed_story); */

/*   t1 = time('\0'); */
/*   printf("\n\nversion %s : command %s", svgsafe(VERSION), svgsafe(COMMAND)); */
/*   printf("\n\nTime: %ld seconds\n", (long)(t1 - t0)); */
/* } */
