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
short int adjacent_dinucleotide_model_init(struct adjacent_dinucleotide_model *d, char *name, short int width)
{
  short int dinucleotide;
  short int first;
  short int second;
  (*d).width = width + 1;
  (*d).name = malloc(1000);
  strcpy((*d).name, name);
  (*d).fraction = malloc(sizeof(double *) * 16 + 5);
  (*d).mononuc_fraction = malloc(sizeof(double *) * 4 + 5);
  for (first = 0; first < 16; first++) {
    (*d).fraction[first] = malloc(sizeof(double) * (width + 1) + 5);
    for (second = 0; second <= width; second++)
      (*d).fraction[first][second] = 0;
  }
  for (first = 0; first < 4; first++) {
    (*d).mononuc_fraction[first] = malloc(sizeof(double) * (width + 1) + 5);
    for (second = 0; second <= width; second++)
      (*d).mononuc_fraction[first][second] = 0;
  }
  return (0);
}


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
  //  short int counter2;
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
  //  signed long int counter2;

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
  return 0;
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
  //  short int CG_difference;
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

/* SUBROUTINE THAT GENERATES A RIVERLAKE LOGO FILE FOR ADMS */
short int Svg_riverlake_logo(char *filename, long int offset, long int yoffset, struct adjacent_dinucleotide_model *a,
			     double mononucleotide_frequency_cutoff, double log_fold_change_cutoff, double absolute_deviation_cutoff,
			     double gray_dinucleotide_cutoff)
{
  absolute_deviation_cutoff = -1;
  FILE *outfile;
  outfile = fopen(filename, "w");

  short int rivers;
  double width;
  double tot_deviation;
  char *rivercolor;
  short int max_riverwidth = 20;
  short int max_radius = 20;
  short int nucleotide_height = max_radius * 1.2;
  short int nucleotide_width = max_radius * 1.7;
  short int top_position;
  short int counter;
  short int first;
  short int second;
  double starty;
  double endy;
  short int nucleotide_value;
  short int pwm_position = 0;
  short int warning = 0;
  short int shift_right = 0;
  double **order;
  order = malloc(sizeof(double *) * 4 + 5);
  for (counter = 0; counter < 3; counter++)
    order[counter] = malloc(sizeof(double) * 6 + 5);
  double **position_memory;
  position_memory = malloc(sizeof(double *) * 4 + 5);
  for (counter = 0; counter < 3; counter++) {
    position_memory[counter] = malloc(sizeof(double) * 6 + 5);
  }
  for (counter = 0; counter < 4; counter++)
    position_memory[1][counter] = 0;

  double swap;
  double font_position;
  char *forward;
  if (rna == 0)
    forward = dnaforward;
  else
    forward = rnaforward;

  char **colors;
  colors = malloc(sizeof(char *) * 7 + 5);
  for (counter = 0; counter < 7; counter++)
    colors[counter] = malloc(sizeof(char) * 20 + 5);
  strcpy(colors[0], "green");
  strcpy(colors[1], "blue");
  strcpy(colors[2], "orange");
  strcpy(colors[3], "red");
  strcpy(colors[4], "black");
  strcpy(colors[5], "black");

  short int total_width = nucleotide_width * ((*a).width - 1) + 2*(offset + max_radius);
  short int total_height = (nucleotide_height * 3 + 2*(yoffset + max_radius)) * 1.25;
  char *font = "Courier";
  fprintf(outfile,
	  "<?xml version=\"1.0\" standalone=\"no\"?> <!DOCTYPE svg PUBLIC \"-//W3C//DTD SVG 1.1//EN\" \"http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd\">");
  fprintf(outfile, "<!--%s : command %s -->\n", svgsafe(VERSION), svgsafe(COMMAND));

  if (noname == 1) {
    fprintf(outfile, "<svg ");
    //    fprintf(outfile, "%i", (*a).width);
    fprintf(outfile, "width=\"%i\" height=\"%i\" ", total_width, total_height);   // JARKKO 26.4.2017
    fprintf(outfile,
	    "x=\"0\" y=\"%i\" version=\"1.1\" xmlns=\"http://www.w3.org/2000/svg\" xmlns:xlink=\"http://www.w3.org/1999/xlink\">",
	    300);
  }

  else
    fprintf(outfile,
	    "<svg width=\"2000\" height=\"200\" x=\"0\" y=\"%i\" version=\"1.1\" xmlns=\"http://www.w3.org/2000/svg\" xmlns:xlink=\"http://www.w3.org/1999/xlink\">",
	    300);

  /* SCALES SVG RIVERLAKE LOGO */
  fprintf(outfile, "<title>%s</title><g transform=\"scale(%f,%f)\">", (*a).name, 1.0, 1.25);

  //fprintf(outfile, "<defs><filter id=\"fractalnoise\" in=\"SourceGraphic\"> <feTurbulence type=\"fractalNoise\" baseFrequency=\"0.4\" numOctaves=\"4\"/></filter></defs>");
  //filterUnits=\"objectBoundingBox\" x=\"0%%\" y=\"0%%\" width=\"100%%\" height=\"100%%\"

  Add_nucleotide_paths(outfile);	/* Adds nucleotide paths */

  /* GENERATES LOGO */
  top_position = 20;

  for (pwm_position = 0; pwm_position < (*a).width; pwm_position++) {

    /* USES ALPHABETIC ORDER (NO SORT), AND SHIFTS BASE GRAPHICAL POSITION MEMORY */
    for (counter = 0; counter < 4; counter++) {
      order[0][counter] = counter;
      order[1][counter] = (*a).mononuc_fraction[counter][pwm_position];
    }

    /* DETERMINES IF ADJACENT BASES AFFECT EACH OTHER MORE THAN CUTOFF */
    //for(shift_right = 0, first = 0; first < 4; first++) if((*a).mononuc_fraction[first][pwm_position] > mononucleotide_frequency_cutoff) for(second = 0; second < 4; second++) if((*a).mononuc_fraction[second][pwm_position+1] > mononucleotide_frequency_cutoff && (((log10((*a).fraction[first*4+second][pwm_position] / ((*a).mononuc_fraction[second][pwm_position+1] + pseudocount)))) > log_fold_change_cutoff || (*a).fraction[first*4+second][pwm_position] - (*a).mononuc_fraction[second][pwm_position+1] > absolute_deviation_cutoff)) shift_right = 1;

    shift_right = 1;

    if (shift_right == 1) {
      for (starty = max_radius, first = 0; first < 4; first++) {

	//(*a).mononuc_fraction[first][pwm_position] * 50;
	if ((*a).mononuc_fraction[first][pwm_position] > mononucleotide_frequency_cutoff) {
	  printf("\nBase %c at position %i over cutoff", forward[first], pwm_position);
	  for (endy = max_radius, second = 0; second < 4; second++) {
	    /* DRAWS LINE TO CONNECT PREFERENTIAL PAIR AND SETS SHIFT RIGHT FLAG IF DEVIATION FROM PWM IS DETECTED */
	    if ((*a).mononuc_fraction[second][pwm_position + 1] > mononucleotide_frequency_cutoff) {
	      printf("\nDinucleotide %c%c at position %i over both cutoffs; cond %.3f vs mono %.3f: log fold %.3f", forward[first], forward[second],
		     pwm_position, (*a).fraction[first * 4 + second][pwm_position], ((*a).mononuc_fraction[second][pwm_position + 1] + pseudocount),
		     (log10((*a).fraction[first * 4 + second][pwm_position] / ((*a).mononuc_fraction[second][pwm_position + 1] + pseudocount))));
	      // HUOM. ALLA OLEVA ON AINA TOSI, SILLÄ FUNKTION ALUSSA ASETETAAN absolute_deviation_cutoff = -1
	      if (((log10((*a).fraction[first * 4 + second][pwm_position] / ((*a).mononuc_fraction[second][pwm_position + 1] + pseudocount)))) >
		  log_fold_change_cutoff
		  || (*a).fraction[first * 4 + second][pwm_position] - (*a).mononuc_fraction[second][pwm_position + 1] > absolute_deviation_cutoff) {
		printf("\t ***over LOGFOLD or ABSDEV cutoff");
		shift_right = 1;
		tot_deviation =
		    ((*a).fraction[first * 4 + second][pwm_position] -
		     (*a).mononuc_fraction[second][pwm_position + 1]) * (*a).mononuc_fraction[first][pwm_position] * 2; // KERROIN 2 TÄYSIN TURHA

		/* DRAWS WIDER RIVER */
		width = (*a).mononuc_fraction[second][pwm_position + 1] * (*a).mononuc_fraction[first][pwm_position] * max_riverwidth;   // KULTA
		if (tot_deviation >= 0)
		  width = (*a).fraction[first * 4 + second][pwm_position] * (*a).mononuc_fraction[first][pwm_position] * max_riverwidth; // HARMAA
		fprintf(outfile, "<g><title>observed %.2f expected %.2f</title><line x1=\"%.2f\" y1=\"%.2f\" x2=\"%.2f\" y2=\"%.2f\" ",
			(*a).fraction[first * 4 + second][pwm_position] * (*a).mononuc_fraction[first][pwm_position],
			(*a).mononuc_fraction[second][pwm_position + 1] * (*a).mononuc_fraction[first][pwm_position],
			(float)pwm_position * nucleotide_width + offset + max_radius, starty,
			(float)(pwm_position + 1) * nucleotide_width + offset + max_radius, endy);
		/* SETS WIDER RIVER COLOR AND DASH */
		rivercolor = "#A0B4CE";                                                                                                  // BLUEISH GRAY
		if (tot_deviation < 0) {
		  rivercolor = "gold";
		  fprintf(outfile, "stroke-dasharray=\"2,2\" stroke-opacity=\"0.25\" ");
		}
		// filter=\"url(#fractalnoise)\" stroke-dasharray=\"0.25,0.75\"
		fprintf(outfile, "stroke = \"%s\" stroke-width = \"%f\"/></g>\n", rivercolor, width);

		/* DRAWS DEEPER RIVER */
		width = (tot_deviation / 2) * max_riverwidth;                                                                            // SININEN
		if (tot_deviation < 0)
		  width = (*a).fraction[first * 4 + second][pwm_position] * (*a).mononuc_fraction[first][pwm_position] * max_riverwidth; // HARMAA
		fprintf(outfile, "<g><title>observed %.2f expected %.2f</title><polyline points =\"%.2f,%.2f %.2f,%.2f\" ",
			(*a).fraction[first * 4 + second][pwm_position] * (*a).mononuc_fraction[first][pwm_position],
			(*a).mononuc_fraction[second][pwm_position + 1] * (*a).mononuc_fraction[first][pwm_position],
			(float)pwm_position * nucleotide_width + offset + max_radius, starty,
			(float)(pwm_position + 1) * nucleotide_width + offset + max_radius, endy);
		/* SETS DEEPER RIVER COLOR */
		printf("\n**Total deviation at dinucleotide %c%c is %.2f", forward[first], forward[second], tot_deviation);
		rivercolor = "blue";
		if (tot_deviation < 0)
		  rivercolor = "#A0B4CE";
		// if (tot_deviation > 0.50) rivercolor = "#00A4BE";
		// if (tot_deviation < -0.25) rivercolor = "lemonchiffon";
		// if (tot_deviation < -0.50) rivercolor = "gold";
		fprintf(outfile, "fill = \"none\" stroke = \"%s\" stroke-width = \"%f\"/></g>\n", rivercolor, width);

	      } // end if 
	      /*else 
	         if ((*a).fraction[first*4+second][pwm_position] * (*a).mononuc_fraction[first][pwm_position] > gray_dinucleotide_cutoff) 
	         {
	         fprintf(outfile, "<polyline points =\"");
	         fprintf(outfile, "%li,%.0f %li,%.0f\" ", pwm_position * 20 + offset +22, starty, pwm_position * 20 + offset +38, endy); 
	         fprintf(outfile, "fill = \"none\" stroke = \"lightgrey\" stroke-width = \"%f\"/>\n", (*a).fraction[first*4+second][pwm_position] * (*a).mononuc_fraction[first][pwm_position] * 15);   
	         } */

	    }
	    endy += nucleotide_height;
	  } // end for second

	  //printf("\nDEBUG--pwm position %i, base %c value %.6f", pwm_position, dnaforward[first], order[1][first]);
	  if (order[1][first] > 0) {
	    /* PRINTS LAKES */
	    fprintf(outfile, "<g><title>%.2f</title><circle cx=\"%.2f\" cy=\"%.2f\" r=\"%.2f\" fill=\"%s\" stroke-width=\"%.2f\"/></g>\n",
		    order[1][first], (float)pwm_position * nucleotide_width + offset + max_radius,
		    (float)first * nucleotide_height + yoffset + max_radius, (float)order[1][first] * max_radius, "lightsteelblue", (float)0);
	    // font_position += (order[1][nucleotide_value] * 100);
	    /* PRINTS OUT SCALED NUCLEOTIDE LABELS */
	    fprintf(outfile, "<use xlink:href=\"#%c\" ", forward[(int)order[0][first]]);
	    fprintf(outfile, " transform=\"translate(%f,%f) scale(%f,%f)\" visibility=\"visible\" />\n",
		    (float)pwm_position * nucleotide_width + offset + max_radius - 10 * order[1][first],
		    (float)first * nucleotide_height + yoffset + max_radius + 10 * order[1][first], order[1][first] * 2, order[1][first] * 2);
	  }

	}
	starty += nucleotide_height;
      } // end for first
      //offset += 20;
    } // end if shift_right

  }  // for pwm_position

  /* PRINTS OUT NAME OF LOGO AND OTHER DATA */
  fprintf(outfile, "</g>");
  if (noname == 0) {
    fprintf(outfile, "<text  x=\"");
    fprintf(outfile, "%li", pwm_position * nucleotide_width + 20 + offset);
    fprintf(outfile, "\" y=\"");
    fprintf(outfile, "%i", top_position + 30);
    fprintf(outfile, "\" fill = \"black\" stroke = \"%s\" font-size=\"30\" font-family = \"", colors[4 - 2 * warning + 3 * (warning == 2)]);
    fprintf(outfile, "%s", font);
    fprintf(outfile, "\" >");
    fprintf(outfile, "%s", (*a).name);
    fprintf(outfile, "</text>\n");
  }
  
  fprintf(outfile, "</svg>");
  fclose(outfile);
  return (0);
}

/* SUBROUTINE THAT GENERATES A RIVERLAKE LOGO FILE FOR ADMS */
short int Svg_riverlake_logo_jarkko(char *filename, long int offset, long int yoffset, struct adjacent_dinucleotide_model *a,
				    double mononucleotide_frequency_cutoff, double log_fold_change_cutoff, double absolute_deviation_cutoff,
				    double gray_dinucleotide_cutoff,
				    int use_constant_radius,
				    int use_transition_probability_width)
{
  absolute_deviation_cutoff = -1;
  FILE *outfile;
  outfile = fopen(filename, "w");

  short int rivers;
  double width;
  double tot_deviation;
  char *rivercolor;
  short int max_riverwidth = 20;
  short int max_radius = 20;
  short int constant_radius = 5;  // Jarkko
 
  short int nucleotide_height = max_radius * 1.2;
  short int nucleotide_width = max_radius * 1.7;
  short int top_position;
  short int counter;
  short int first;
  short int second;
  double starty;
  double endy;
  short int nucleotide_value;
  short int pwm_position = 0;
  short int warning = 0;
  //  short int shift_right = 0;
  double **order;
  order = malloc(sizeof(double *) * 4 + 5);
  for (counter = 0; counter < 3; counter++)
    order[counter] = malloc(sizeof(double) * 6 + 5);
  double **position_memory;
  position_memory = malloc(sizeof(double *) * 4 + 5);
  for (counter = 0; counter < 3; counter++) {
    position_memory[counter] = malloc(sizeof(double) * 6 + 5);
  }
  for (counter = 0; counter < 4; counter++)
    position_memory[1][counter] = 0;

  double swap;
  double font_position;
  char *forward;
  if (rna == 0)
    forward = dnaforward;
  else
    forward = rnaforward;

  char **colors;
  colors = malloc(sizeof(char *) * 7 + 5);
  for (counter = 0; counter < 7; counter++)
    colors[counter] = malloc(sizeof(char) * 20 + 5);
  strcpy(colors[0], "green");
  strcpy(colors[1], "blue");
  strcpy(colors[2], "orange");
  strcpy(colors[3], "red");
  strcpy(colors[4], "black");
  strcpy(colors[5], "black");

  short int total_width = nucleotide_width * ((*a).width - 1) + 2*(offset + max_radius);
  short int total_height = (nucleotide_height * 3 + 2*(yoffset + max_radius)) * 1.25;
  char *font = "Courier";
  fprintf(outfile,
	  "<?xml version=\"1.0\" standalone=\"no\"?> <!DOCTYPE svg PUBLIC \"-//W3C//DTD SVG 1.1//EN\" \"http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd\">");
  fprintf(outfile, "<!--%s : command %s -->\n", svgsafe(VERSION), svgsafe(COMMAND));

  if (noname == 1) {
    fprintf(outfile, "<svg ");
    //    fprintf(outfile, "%i", (*a).width);
    fprintf(outfile, "width=\"%i\" height=\"%i\" ", total_width, total_height);   // JARKKO 26.4.2017
    fprintf(outfile,
	    "x=\"0\" y=\"%i\" version=\"1.1\" xmlns=\"http://www.w3.org/2000/svg\" xmlns:xlink=\"http://www.w3.org/1999/xlink\">",
	    300);
  }

  else
    fprintf(outfile,
	    "<svg width=\"2000\" height=\"200\" x=\"0\" y=\"%i\" version=\"1.1\" xmlns=\"http://www.w3.org/2000/svg\" xmlns:xlink=\"http://www.w3.org/1999/xlink\">",
	    300);

  /* SCALES SVG RIVERLAKE LOGO */
  fprintf(outfile, "<title>%s</title><g transform=\"scale(%f,%f)\">", (*a).name, 1.0, 1.25);

  //fprintf(outfile, "<defs><filter id=\"fractalnoise\" in=\"SourceGraphic\"> <feTurbulence type=\"fractalNoise\" baseFrequency=\"0.4\" numOctaves=\"4\"/></filter></defs>");
  //filterUnits=\"objectBoundingBox\" x=\"0%%\" y=\"0%%\" width=\"100%%\" height=\"100%%\"

  Add_nucleotide_paths(outfile);	/* Adds nucleotide paths */

  /* GENERATES LOGO */
  top_position = 20;

  for (pwm_position = 0; pwm_position < (*a).width; pwm_position++) {

    /* USES ALPHABETIC ORDER (NO SORT), AND SHIFTS BASE GRAPHICAL POSITION MEMORY */
    for (counter = 0; counter < 4; counter++) {
      order[0][counter] = counter;
      order[1][counter] = (*a).mononuc_fraction[counter][pwm_position];
    }


    for (starty = max_radius, first = 0; first < 4; first++) {

      //(*a).mononuc_fraction[first][pwm_position] * 50;
      if ((*a).mononuc_fraction[first][pwm_position] > mononucleotide_frequency_cutoff || use_transition_probability_width) {
	printf("\nBase %c at position %i over cutoff", forward[first], pwm_position);
	for (endy = max_radius, second = 0; second < 4; second++) {
	  /* DRAWS LINE TO CONNECT PREFERENTIAL PAIR AND SETS SHIFT RIGHT FLAG IF DEVIATION FROM PWM IS DETECTED */
	  if ((*a).mononuc_fraction[second][pwm_position + 1] > mononucleotide_frequency_cutoff || use_transition_probability_width) {
	    printf("\nDinucleotide %c%c at position %i over both cutoffs; cond %.3f vs mono %.3f: log fold %.3f", forward[first], forward[second],
		   pwm_position, (*a).fraction[first * 4 + second][pwm_position], ((*a).mononuc_fraction[second][pwm_position + 1] + pseudocount),
		   (log10((*a).fraction[first * 4 + second][pwm_position] / ((*a).mononuc_fraction[second][pwm_position + 1] + pseudocount))));
	    // HUOM. ALLA OLEVA ON AINA TOSI, SILLÄ FUNKTION ALUSSA ASETETAAN absolute_deviation_cutoff = -1
	    if (((log10((*a).fraction[first * 4 + second][pwm_position] / ((*a).mononuc_fraction[second][pwm_position + 1] + pseudocount)))) >
		log_fold_change_cutoff
		|| (*a).fraction[first * 4 + second][pwm_position] - (*a).mononuc_fraction[second][pwm_position + 1] > absolute_deviation_cutoff) {
	      printf("\t ***over LOGFOLD or ABSDEV cutoff");

	      tot_deviation =
		((*a).fraction[first * 4 + second][pwm_position] -
		 (*a).mononuc_fraction[second][pwm_position + 1]) * (*a).mononuc_fraction[first][pwm_position] * 2; // KERROIN 2 TÄYSIN TURHA

	      /* DRAWS WIDER RIVER */
	      float multiplier = use_transition_probability_width ? 0.5 : (*a).mononuc_fraction[first][pwm_position];
	      width = (*a).fraction[first * 4 + second][pwm_position] * multiplier * max_riverwidth; // HARMAA
	      fprintf(outfile, "<g><title>observed %.2f expected %.2f</title><line x1=\"%.2f\" y1=\"%.2f\" x2=\"%.2f\" y2=\"%.2f\" ",
		      (*a).fraction[first * 4 + second][pwm_position] * (*a).mononuc_fraction[first][pwm_position],
		      (*a).mononuc_fraction[second][pwm_position + 1] * (*a).mononuc_fraction[first][pwm_position],
		      (float)pwm_position * nucleotide_width + offset + max_radius, starty,
		      (float)(pwm_position + 1) * nucleotide_width + offset + max_radius, endy);
	      /* SETS WIDER RIVER COLOR AND DASH */
	      rivercolor = "#A0B4CE";                                                                                                  // BLUEISH GRAY
	      // filter=\"url(#fractalnoise)\" stroke-dasharray=\"0.25,0.75\"
	      fprintf(outfile, "stroke = \"%s\" stroke-width = \"%f\"/></g>\n", rivercolor, width);


	    } // end if 

	  }
	  endy += nucleotide_height;
	} // end for second
      } // end if

      if ((*a).mononuc_fraction[first][pwm_position] > mononucleotide_frequency_cutoff || use_constant_radius) {
	  //printf("\nDEBUG--pwm position %i, base %c value %.6f", pwm_position, dnaforward[first], order[1][first]);
	if (order[1][first] > 0) {
	  /* PRINTS LAKES */
	  float radius = use_constant_radius ? constant_radius : (float)order[1][first] * max_radius; 
	  float radius2 = use_constant_radius ? 0.1*constant_radius : (float)order[1][first] * max_radius; 
	  fprintf(outfile, "<g><title>%.2f</title><circle cx=\"%.2f\" cy=\"%.2f\" r=\"%.2f\" fill=\"%s\" stroke-width=\"%.2f\"/></g>\n",
		  order[1][first],
		  (float)pwm_position * nucleotide_width + offset + max_radius,
		  (float)first * nucleotide_height + yoffset + max_radius,
		  radius,
		  "lightsteelblue", (float)0);
	  // font_position += (order[1][nucleotide_value] * 100);
	  /* PRINTS OUT SCALED NUCLEOTIDE LABELS */
	  fprintf(outfile, "<use xlink:href=\"#%c\" ", forward[(int)order[0][first]]);
	  if (use_constant_radius) 
	    fprintf(outfile, " transform=\"translate(%f,%f) scale(%f,%f)\" visibility=\"visible\" />\n",
		    (float)pwm_position * nucleotide_width + offset + max_radius - constant_radius/2.0,
		    (float)first * nucleotide_height + yoffset + max_radius + constant_radius/2.0,
		    radius2,
		    radius2);
	  else {
	    float temp = order[1][first];
	    fprintf(outfile, " transform=\"translate(%f,%f) scale(%f,%f)\" visibility=\"visible\" />\n",
		    (float)pwm_position * nucleotide_width + offset + max_radius - 10 * temp,
		    (float)first * nucleotide_height + yoffset + max_radius + 10 * temp, temp * 2, temp * 2);
	  }
	} // end if

      } // end if
      
      starty += nucleotide_height;
    } // end for first

  }  // for pwm_position

  /* PRINTS OUT NAME OF LOGO AND OTHER DATA */
  fprintf(outfile, "</g>");
  if (noname == 0) {
    fprintf(outfile, "<text  x=\"");
    fprintf(outfile, "%li", pwm_position * nucleotide_width + 20 + offset);
    fprintf(outfile, "\" y=\"");
    fprintf(outfile, "%i", top_position + 30);
    fprintf(outfile, "\" fill = \"black\" stroke = \"%s\" font-size=\"30\" font-family = \"", colors[4 - 2 * warning + 3 * (warning == 2)]);
    fprintf(outfile, "%s", font);
    fprintf(outfile, "\" >");
    fprintf(outfile, "%s", (*a).name);
    fprintf(outfile, "</text>\n");
  }
  
  fprintf(outfile, "</svg>");
  fclose(outfile);
  return (0);
} // Svg_riverlake_logo_jarkko


/* SUBROUTINE THAT GENERATES A RIVERLAKE LOGO FILE FOR ADMS */
short int Svg_riverlake_logo_jarkko_diff(char *filename, long int offset, long int yoffset, struct adjacent_dinucleotide_model *a,
					 double mononucleotide_frequency_cutoff, double log_fold_change_cutoff, double absolute_deviation_cutoff,
					 double gray_dinucleotide_cutoff,
					 int use_constant_radius,
					 int use_transition_probability_width)
{
  absolute_deviation_cutoff = -1;
  FILE *outfile;
  outfile = fopen(filename, "w");

  short int rivers;
  double width;
  double tot_deviation;
  char *rivercolor;
  short int max_riverwidth = 20;
  short int max_radius = 20;
  short int constant_radius = 5;  // Jarkko
 
  short int nucleotide_height = max_radius * 1.2;
  short int nucleotide_width = max_radius * 1.7;
  short int top_position;
  short int counter;
  short int first;
  short int second;
  double starty;
  double endy;
  short int nucleotide_value;
  short int pwm_position = 0;
  short int warning = 0;
  //  short int shift_right = 0;
  double **order;
  order = malloc(sizeof(double *) * 4 + 5);
  for (counter = 0; counter < 3; counter++)
    order[counter] = malloc(sizeof(double) * 6 + 5);
  double **position_memory;
  position_memory = malloc(sizeof(double *) * 4 + 5);
  for (counter = 0; counter < 3; counter++) {
    position_memory[counter] = malloc(sizeof(double) * 6 + 5);
  }
  for (counter = 0; counter < 4; counter++)
    position_memory[1][counter] = 0;

  double swap;
  double font_position;
  char *forward;
  if (rna == 0)
    forward = dnaforward;
  else
    forward = rnaforward;

  char **colors;
  colors = malloc(sizeof(char *) * 7 + 5);
  for (counter = 0; counter < 7; counter++)
    colors[counter] = malloc(sizeof(char) * 20 + 5);
  strcpy(colors[0], "green");
  strcpy(colors[1], "blue");
  strcpy(colors[2], "orange");
  strcpy(colors[3], "red");
  strcpy(colors[4], "black");
  strcpy(colors[5], "black");

  short int total_width = nucleotide_width * ((*a).width - 1) + 2*(offset + max_radius);
  short int total_height = (nucleotide_height * 3 + 2*(yoffset + max_radius)) * 1.25;
  char *font = "Courier";
  fprintf(outfile,
	  "<?xml version=\"1.0\" standalone=\"no\"?> <!DOCTYPE svg PUBLIC \"-//W3C//DTD SVG 1.1//EN\" \"http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd\">");
  fprintf(outfile, "<!--%s : command %s -->\n", svgsafe(VERSION), svgsafe(COMMAND));

  if (noname == 1) {
    fprintf(outfile, "<svg ");
    //    fprintf(outfile, "%i", (*a).width);
    fprintf(outfile, "width=\"%i\" height=\"%i\" ", total_width, total_height);   // JARKKO 26.4.2017
    fprintf(outfile,
	    "x=\"0\" y=\"%i\" version=\"1.1\" xmlns=\"http://www.w3.org/2000/svg\" xmlns:xlink=\"http://www.w3.org/1999/xlink\">",
	    300);
  }

  else
    fprintf(outfile,
	    "<svg width=\"2000\" height=\"200\" x=\"0\" y=\"%i\" version=\"1.1\" xmlns=\"http://www.w3.org/2000/svg\" xmlns:xlink=\"http://www.w3.org/1999/xlink\">",
	    300);

  /* SCALES SVG RIVERLAKE LOGO */
  fprintf(outfile, "<title>%s</title><g transform=\"scale(%f,%f)\">", (*a).name, 1.0, 1.25);

  //fprintf(outfile, "<defs><filter id=\"fractalnoise\" in=\"SourceGraphic\"> <feTurbulence type=\"fractalNoise\" baseFrequency=\"0.4\" numOctaves=\"4\"/></filter></defs>");
  //filterUnits=\"objectBoundingBox\" x=\"0%%\" y=\"0%%\" width=\"100%%\" height=\"100%%\"

  Add_nucleotide_paths(outfile);	/* Adds nucleotide paths */

  /* GENERATES LOGO */
  top_position = 20;

  for (pwm_position = 0; pwm_position < (*a).width; pwm_position++) {

    /* USES ALPHABETIC ORDER (NO SORT), AND SHIFTS BASE GRAPHICAL POSITION MEMORY */
    for (counter = 0; counter < 4; counter++) {
      order[0][counter] = counter;
      order[1][counter] = (*a).mononuc_fraction[counter][pwm_position];
    }

    /* DETERMINES IF ADJACENT BASES AFFECT EACH OTHER MORE THAN CUTOFF */
    //for(shift_right = 0, first = 0; first < 4; first++) if((*a).mononuc_fraction[first][pwm_position] > mononucleotide_frequency_cutoff) for(second = 0; second < 4; second++) if((*a).mononuc_fraction[second][pwm_position+1] > mononucleotide_frequency_cutoff && (((log10((*a).fraction[first*4+second][pwm_position] / ((*a).mononuc_fraction[second][pwm_position+1] + pseudocount)))) > log_fold_change_cutoff || (*a).fraction[first*4+second][pwm_position] - (*a).mononuc_fraction[second][pwm_position+1] > absolute_deviation_cutoff)) shift_right = 1;

    for (starty = max_radius, first = 0; first < 4; first++) {

      //(*a).mononuc_fraction[first][pwm_position] * 50;
      if (fabs((*a).mononuc_fraction[first][pwm_position]) > mononucleotide_frequency_cutoff || use_transition_probability_width) {
	printf("\nBase %c at position %i over cutoff", forward[first], pwm_position);
	for (endy = max_radius, second = 0; second < 4; second++) {
	  /* DRAWS LINE TO CONNECT PREFERENTIAL PAIR AND SETS SHIFT RIGHT FLAG IF DEVIATION FROM PWM IS DETECTED */
	  if (fabs((*a).mononuc_fraction[second][pwm_position + 1]) > mononucleotide_frequency_cutoff || use_transition_probability_width) {
	    printf("\nDinucleotide %c%c at position %i over both cutoffs; cond %.3f vs mono %.3f: log fold %.3f", forward[first], forward[second],
		   pwm_position, (*a).fraction[first * 4 + second][pwm_position], ((*a).mononuc_fraction[second][pwm_position + 1] + pseudocount),
		   (log10((*a).fraction[first * 4 + second][pwm_position] / ((*a).mononuc_fraction[second][pwm_position + 1] + pseudocount))));
	    // HUOM. ALLA OLEVA ON AINA TOSI, SILLÄ FUNKTION ALUSSA ASETETAAN absolute_deviation_cutoff = -1
	    /*
	    if (((log10((*a).fraction[first * 4 + second][pwm_position] / ((*a).mononuc_fraction[second][pwm_position + 1] + pseudocount)))) >
		log_fold_change_cutoff
		|| (*a).fraction[first * 4 + second][pwm_position] - (*a).mononuc_fraction[second][pwm_position + 1] > absolute_deviation_cutoff) {
	    */
	    if (1) {
	      printf("\t ***over LOGFOLD or ABSDEV cutoff");

	      /*
	      tot_deviation =
		((*a).fraction[first * 4 + second][pwm_position] -
		 (*a).mononuc_fraction[second][pwm_position + 1]) * (*a).mononuc_fraction[first][pwm_position] * 2; // KERROIN 2 TÄYSIN TURHA
	      */
	      
	      /* DRAWS WIDER RIVER */
	      //float multiplier = use_transition_probability_width ? 0.5 : (*a).mononuc_fraction[first][pwm_position];
	      width = fabs((*a).fraction[first * 4 + second][pwm_position]) * max_riverwidth; // HARMAA
	      fprintf(outfile, "<g><title>observed %.2f expected %.2f</title><line x1=\"%.2f\" y1=\"%.2f\" x2=\"%.2f\" y2=\"%.2f\" ",
		      (*a).fraction[first * 4 + second][pwm_position] * (*a).mononuc_fraction[first][pwm_position],
		      (*a).mononuc_fraction[second][pwm_position + 1] * (*a).mononuc_fraction[first][pwm_position],
		      (float)pwm_position * nucleotide_width + offset + max_radius, starty,
		      (float)(pwm_position + 1) * nucleotide_width + offset + max_radius, endy);
	      /* SETS WIDER RIVER COLOR AND DASH */
	      rivercolor = (*a).fraction[first * 4 + second][pwm_position] >= 0.0 ? "#A0B4CE" : "purple";                          // BLUEISH GRAY
	      // filter=\"url(#fractalnoise)\" stroke-dasharray=\"0.25,0.75\"
	      fprintf(outfile, "stroke = \"%s\" stroke-width = \"%f\"/></g>\n", rivercolor, width);


	    } // end if 

	  }
	  endy += nucleotide_height;
	} // end for second
      } // end if

      if (fabs((*a).mononuc_fraction[first][pwm_position]) > mononucleotide_frequency_cutoff || use_constant_radius) {
	  //printf("\nDEBUG--pwm position %i, base %c value %.6f", pwm_position, dnaforward[first], order[1][first]);
	if (order[1][first] != 0.0) {
	  /* PRINTS LAKES */
	  float radius = use_constant_radius ? constant_radius : (float)order[1][first] * max_radius; 
	  float radius2 = use_constant_radius ? 0.1*constant_radius : (float)order[1][first] * max_radius; 
	  fprintf(outfile, "<g><title>%.2f</title><circle cx=\"%.2f\" cy=\"%.2f\" r=\"%.2f\" fill=\"%s\" stroke-width=\"%.2f\"/></g>\n",
		  order[1][first],
		  (float)pwm_position * nucleotide_width + offset + max_radius,
		  (float)first * nucleotide_height + yoffset + max_radius,
		  fabsf(radius),
		  radius >= 0.0 ? "lightsteelblue" : "purple", (float)0);
	  // font_position += (order[1][nucleotide_value] * 100);
	  /* PRINTS OUT SCALED NUCLEOTIDE LABELS */
	  fprintf(outfile, "<use xlink:href=\"#%c\" ", forward[(int)order[0][first]]);
	  if (use_constant_radius) 
	    fprintf(outfile, " transform=\"translate(%f,%f) scale(%f,%f)\" visibility=\"visible\" />\n",
		    (float)pwm_position * nucleotide_width + offset + max_radius - constant_radius/2.0,
		    (float)first * nucleotide_height + yoffset + max_radius + constant_radius/2.0,
		    radius2,
		    radius2);
	  else {
	    float temp = fabs(order[1][first]);
	    fprintf(outfile, " transform=\"translate(%f,%f) scale(%f,%f)\" visibility=\"visible\" />\n",
		    (float)pwm_position * nucleotide_width + offset + max_radius - 10 * temp,
		    (float)first * nucleotide_height + yoffset + max_radius + 10 * temp, temp * 2, temp * 2);
	  }
	} // end if

      } // end if
      
      starty += nucleotide_height;
    } // end for first

  }  // for pwm_position

  /* PRINTS OUT NAME OF LOGO AND OTHER DATA */
  fprintf(outfile, "</g>");
  if (noname == 0) {
    fprintf(outfile, "<text  x=\"");
    fprintf(outfile, "%li", pwm_position * nucleotide_width + 20 + offset);
    fprintf(outfile, "\" y=\"");
    fprintf(outfile, "%i", top_position + 30);
    fprintf(outfile, "\" fill = \"black\" stroke = \"%s\" font-size=\"30\" font-family = \"", colors[4 - 2 * warning + 3 * (warning == 2)]);
    fprintf(outfile, "%s", font);
    fprintf(outfile, "\" >");
    fprintf(outfile, "%s", (*a).name);
    fprintf(outfile, "</text>\n");
  }
  
  fprintf(outfile, "</svg>");
  fclose(outfile);
  printf("\n");
  return (0);
} // Svg_riverlake_logo_jarkko_diff


/* SUBROUTINE THAT SUBTRACTS NORMALIZED c FROM NORMALIZED PWM n */
short int Subtract_normalized_pwm(struct normalized_pwm *n, struct normalized_pwm *c, short int offset)
{
  short int counter;
  short int position;
  //  double total_nucleotides = 0;
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


/* SUBROUTINE THAT LOADS AN ADM */
short int Load_ADM(struct adjacent_dinucleotide_model *a, char *filename)
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
  for (line = 0; line < 16;) {
    for (counter = 0; counter < 30; counter++) {
      text1 = getc(pwmfile);
      if (text1 == EOF || text1 == '\n' || text1 == '\t') {
	current_string[counter] = '\0';
	if (counter > 0
	    && (current_string[0] == '0' || current_string[0] == '1' || current_string[0] == '2' || current_string[0] == '3'
		|| current_string[0] == '4' || current_string[0] == '5' || current_string[0] == '6' || current_string[0] == '7'
		|| current_string[0] == '8' || current_string[0] == '9' || current_string[0] == ' ' || current_string[0] == '-')) {
	  (*a).fraction[line][pwm_position] = atof(current_string);
	  printf("%f\t", (*a).fraction[line][pwm_position]);
	  fflush(stdout);
	  pwm_position++;
	}
	if (text1 == '\n' || text1 == EOF) {
	  (*a).width = pwm_position;
	  line++;
	  pwm_position = 0;
	  printf("\n");
	}
	break;
      }
      current_string[counter] = text1;
      /* printf ("%c", text1); */
    }
  }
  pwm_position = 0;
  for (line = 0; line < 4;) {
    for (counter = 0; counter < 30; counter++) {
      text1 = getc(pwmfile);
      if (text1 == EOF || text1 == '\n' || text1 == '\t') {
	current_string[counter] = '\0';
	if (counter > 0
	    && (current_string[0] == '0' || current_string[0] == '1' || current_string[0] == '2' || current_string[0] == '3'
		|| current_string[0] == '4' || current_string[0] == '5' || current_string[0] == '6' || current_string[0] == '7'
		|| current_string[0] == '8' || current_string[0] == '9' || current_string[0] == ' ' || current_string[0] == '-')) {
	  (*a).mononuc_fraction[line][pwm_position] = atof(current_string);
	  printf("%f\t", (*a).mononuc_fraction[line][pwm_position]);
	  fflush(stdout);
	  pwm_position++;
	}
	if (text1 == '\n' || text1 == EOF) {
	  (*a).width = pwm_position;
	  line++;
	  pwm_position = 0;
	  printf("\n");
	}
	break;
      }
      current_string[counter] = text1;
      /* printf ("%c", text1); */
    }
  }

  free(current_string);
  //Print_ADM(a);
  if (text1 == EOF && line != 19)
    return (1);
  else
    return (0);
}

/* SUBROUTINE THAT SUBTRACTS NORMALIZED c FROM NORMALIZED ADM n */
short int Subtract_normalized_ADM(struct adjacent_dinucleotide_model *n, struct adjacent_dinucleotide_model *c)
{
  short int counter;
  short int position;
  //  double total_nucleotides = 0;
  printf("Subtracted ADM\n");
  for (counter = 0; counter < 16; counter++) {
    for (position = 0; position < (*n).width-1; position++) {// difference of dinucleotide probabilities
      (*n).fraction[counter][position] = (*n).mononuc_fraction[counter/4][position] * (*n).fraction[counter][position] - (*c).mononuc_fraction[counter/4][position] * (*c).fraction[counter][position];
      printf("%f\t", (*n).fraction[counter][position]);
    }
    printf("\n");
  }
  for (counter = 0; counter < 4; counter++) {
    for (position = 0; position < (*n).width; position++) {
      (*n).mononuc_fraction[counter][position] -= (*c).mononuc_fraction[counter][position];
      printf("%f\t", (*n).mononuc_fraction[counter][position]);
    }
    printf("\n");
  }

  return (0);
}



/* ########################## MAIN PROGRAM ############################# */
/* ########################## MAIN PROGRAM ############################# */
/* ########################## MAIN PROGRAM ############################# */
/* ########################## MAIN PROGRAM ############################# */
/* ########################## MAIN PROGRAM ############################# */
int main(int argc, char *argv[])
{

  Nlength++;

  /* time_t t0, t1; */
  /* t0 = time('\0'); */

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



  long int counter;
  long int counter2;
  long int counter3;

  /* clears pvalue cache */
  for (counter = 0; counter < 1000; counter++)
    for (counter2 = 0; counter2 < 1000; counter2++)
      pvalue_cache[counter][counter2] = 2;

  //  FILE *open_file;
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

  char *linecommand;

  seed_story = malloc(10000 * sizeof(char) + 5);
  strcpy(seed_story, "\n");
  char *current_sequence = malloc(1000);
  strcpy(current_sequence, "INITIALVALUE");
  char *tempstring = malloc(10000);


  /* FORWARD AND REVERSE NUCLEOTIDE COUNT PWM STRUCTURES */
  struct count_pwm *nc;
  nc = malloc(sizeof(struct count_pwm) * 3 + 5);
  count_pwm_init(&nc[0], "EMPTY", max_Nlength, 0);
  count_pwm_init(&nc[1], "EMPTY", max_Nlength, 0);



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


  /* VARIABLES AND STRUCTURES FOR HIT POSITION DATA GENERATION */
  struct hit_position_matrix hit_position;
  hit_position_matrix_init(&hit_position, "hit positions", max_Nlength, 0);

  char *exclude_strand;
  exclude_strand = malloc(sizeof(char *) * max_Nlength + 5);
  short int *exclude_position;
  exclude_position = malloc(sizeof(short int *) * max_Nlength + 5);
  for (counter = 0; counter < max_Nlength; counter++) {
    exclude_position[counter] = 0;
    exclude_strand[counter] = 'N';
  }

  /* FLAGS */

  short int svg_only = 0;

  short int difference_logo = 0;
  short int negative_values_allowed = 0;
  short int offset = 0;


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

	if (strcmp(linecommand, "--logo") == 0)
	  svg_only = 1;
	if (strcmp(linecommand, "--difflogo") == 0) {
	  svg_only = 1;
	  difference_logo = 1;
	  paths = 1;
	}

	if (strcmp(linecommand, "-paths") == 0)
	  paths = 1;

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
  //  short int reverse_strand_position;
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
  struct adjacent_dinucleotide_model unflanked_adm2;  // for subtraction of two ADMs, Jarkko
  struct adjacent_dinucleotide_model flanked_adm;
  struct adjacent_dinucleotide_model *flanked_adm_p;


  /* VARIABLES AND STRUCTURES FOR LOGO GENERATION */
  //  short int current_logo = 0;
  struct normalized_pwm *np;
  np = malloc(sizeof(struct normalized_pwm) * 200 + 5);
  struct normalized_pwm **np_p;
  np_p = malloc(sizeof(struct normalized_pwm *) * 200 + 5);
  for (counter = 0; counter < 200; counter++) {
    normalized_pwm_init(&np[counter], "empty", Nlength * 2, 0);
    np_p[counter] = &np[counter];
  }


/* double **total_relative_deviation;
total_relative_deviation = malloc(sizeof(double *) * Nlength * 2 + 5);
for(counter = 0; counter < Nlength * 2; counter++) total_relative_deviation[counter] = malloc(sizeof(double) * Nlength * 2 + 5); */
  double **uncentered_correlation;
  uncentered_correlation = malloc(sizeof(double *) * Nlength * 2 + 5);
  for (counter = 0; counter < Nlength * 2; counter++)
    uncentered_correlation[counter] = malloc(sizeof(double) * Nlength * 2 + 5);


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
    else {
      printf("\nloads ADM\n");
      adjacent_dinucleotide_model_init(&unflanked_adm, "unflanked_adm", Nlength);
      Load_ADM(&unflanked_adm, searchstring);
      strcpy(unflanked_adm.name, searchstring);

      if (difference_logo == 1) {
	int use_constant_radius = 0;
	int use_transition_probability_width = 0;

	strcpy(tempstring, argv[2 + firstnoncommandposition]);
	printf("\nloads ADM2\n");
	adjacent_dinucleotide_model_init(&unflanked_adm2, "unflanked_adm2", Nlength);
	Load_ADM(&unflanked_adm2, tempstring);
	Subtract_normalized_ADM(&unflanked_adm, &unflanked_adm2);
	//unflanked_adm.negative_values_allowed = 1;
	strcat(unflanked_adm.name, "_minus_");
	strcat(unflanked_adm.name, tempstring);
	strcpy(tempstring, unflanked_adm.name);
	strcat(tempstring, ".svg");
	Svg_riverlake_logo_jarkko_diff(tempstring, 0, 0, &unflanked_adm, 0.03, 0.1, 1000, 0.1, use_constant_radius, use_transition_probability_width);
      }
      else {
	int use_constant_radius = 0;
	int use_transition_probability_width = 0;
	strcpy(tempstring, unflanked_adm.name);
	strcpy(searchstring, unflanked_adm.name);
	strcat(tempstring, ".svg");
	if (argc >= 3 + firstnoncommandposition) {   // JARKKO 28.4.2017
	  strcpy(tempstring, argv[2 + firstnoncommandposition]);
	  strcpy(searchstring, argv[2 + firstnoncommandposition]);
	}
	Svg_riverlake_logo(tempstring, 0, 0, &unflanked_adm, 0.05, 0.1, 1000, 0.1);
	//Svg_riverlake_logo_jarkko(tempstring, 0, 0, &unflanked_adm, 0.05, 0.1, 1000, 0.1, use_constant_radius, use_transition_probability_width);
      }
  }

    char *system_command;
    system_command = malloc(2000);
    strcpy(system_command, "convert ");
    strcat(system_command, tempstring);
    strcat(system_command, " ");
    strcat(system_command, searchstring);
    strcat(system_command, ".png");
//printf("\n%s\n", system_command);
//system(system_command);
    return 0;
  }
  /* SOMETHING WRONG, PRINT INSTRUCTIONS */
  else {
    printf
      ("Usage: \n"
       "./spacek --tool -options PWM1 [ PWM2 | svgoutfile ] \n"
       "\n"
       "Tool:\n"
       "  --logo [PWM file name] [output file name (opt)]\n"
       "\tGenerate SVG logo file from input PWM, also generates .png if convert is installed in path\n"
       "  --difflogo [PWM1 file name] [PWM2 file name] [offset (opt)] \n"
       "\tGenerate difference logo (PWM1-PWM2)\n"
       "Options:\n"
       "\n"
       "Logo and sequence formatting:\n"
       "  -paths\tGenerate PWM logos using paths and not fonts\n"
       "  -neg\tAllow negative values in PWM (using paths)\n"
       "  -core=[number]\tHilight this number of positions in the middle of the logo\n"
       "  -noname\tDo not include filename in the logo\n"
       "  -rna\tInput sequence is from RNA, use uracil (U) instead of thymidine (T) in sequences and logos\n"
       "\n"
       " ");

    exit(0);
  }

  return 0;
}
/* END OF LOGO OPTION */

