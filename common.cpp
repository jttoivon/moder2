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
#include "common.hpp"
#include "iupac.hpp"

#include <libgen.h> // for dirname and basename
#include <fstream>
#include <bitset>
#include <sstream>
#include <cstdarg>
#include <cstring>

#include <boost/tuple/tuple.hpp>

//std::vector<std::string> sequences;
std::string header;

/*
template
boost::tuple<int,int,int,int>
get_ranges<>(const boost::multi_array<double, 2>& a);
*/


double
cut(double x)
{
  if (x < 0)
    return 0;
  else
    return x;
}


const char*
yesno(bool b)
{
  return b ? "yes" : "no";
}

void
error(bool b, const std::string& message)
{
  if (b) {
    fprintf(stderr, "Error: %s\n", message.c_str());
    exit(1);
  }
}

// Number of 1-bits in the range
int
bits_in_range(const std::bitset<64>& occ, int begin, int end)
{
  int count=0;
  for (int j=begin; j < end; ++j)
    count += occ[j];
  return count;
}

template <>
std::string
print_vector<bool>(const std::vector<bool>& v)
{
  if (v.size() == 0)
    return std::string("[]");
  std::ostringstream result;
  result << "[";
  for (int i=0; i < v.size()-1; ++i)
    result << yesno(v[i]) << ", ";
  result << yesno(v[v.size()-1]) << "]";

  return result.str();
}

std::string
print_bitvector(unsigned int x, int number_of_bits)
{
  std::string bit_repr(number_of_bits, '-');
  for (int pos=0; pos < number_of_bits; ++pos) {
    unsigned char mask = 1u << (number_of_bits-pos-1);
    bit_repr[pos] = x & mask ? '1' : '0';
  }
  
  return std::string(bit_repr);
}

std::string
mybasename(const std::string& path)
{
  char *basec = strdup(path.c_str());
  
  std::string result = basename(basec);
  free(basec);
  return result;
}


template <>
std::string
m(std::string s, bool t)
{
  std::ostringstream ss;
  ss << s << ", default: ";
  if (t)
    ss << "yes";
  else
    ss << "no";
  return ss.str();
}

void
print_command_line(int argc, char* argv[], FILE* f)
{
  fprintf(f, "%s", argv[0]);
  for (int i=1; i < argc; ++i)
    fprintf(f, " %s", argv[i]);
  fprintf(f, "\n");
}

void
print_command_line(int argc, char* argv[])
{
  print_command_line(argc, argv, stdout);
}




matrix<double> positional_background(0,0);

// template <typename Value, typename Key>
// bool
// key_sort_functor_asc(Value x, Value y)
// {
//   Key k;
//   return k(x) < k(y); 
// }

std::vector<int>
count_frequencies(const std::string& s)
{
  std::vector<int> result(4, 0);
  for (int i=0; i < s.size(); ++i)
    ++result[to_int(s[i])];
  return result;
}

std::string
to_string(std::string format, ...)
{
  char* tmpstr;

  va_list args;
  va_start( args, format );
  int count = vasprintf(&tmpstr, format.c_str(), args);
  va_end( args );

  assert(count != -1);
  std::string result = tmpstr;
  free(tmpstr);

  return result;
}

std::string
vto_string(const std::string& format, va_list args)
{
  char* tmpstr;

  int count = vasprintf(&tmpstr, format.c_str(), args);

  assert(count != -1);
  std::string result = tmpstr;
  free(tmpstr);

  return result;
}

void
underline(std::string format, ...)
{
  va_list args;
  va_start( args, format );
  std::string temp = vto_string(format, args);
  va_end( args );

  std::cout << temp << std::endl;
  std::cout << std::string(temp.length(), '-') << std::endl;

}

inline
void
add_to_vector(std::vector<int>& to, const std::vector<int>& from)
{
  assert(to.size() == from.size());
  for (size_t i=0; i < to.size(); ++i)
    to[i] += from[i];
  
}

std::vector<int>
count_frequencies(const std::vector<std::string>& sequences)
{
  std::vector<int> result(4, 0);
  for (int i=0; i < sequences.size(); ++i)
    add_to_vector(result, count_frequencies(sequences[i]));
  return result;
}

std::string
itoa(int i)
{
  std::ostringstream s;
  s << i;
  return s.str();
}

int
atoi(const std::string& s)
{
  return atoi(s.c_str());
}

double
atof(const std::string& s)
{
  return atof(s.c_str());
}


void
normalize_vector(std::vector<double>& v)
{
  double s = sum(v);

  for (int i=0; i < v.size(); ++i)
    v[i] = v[i]/s;

  return;
}



std::vector<double>
normalize_vector_copy(const std::vector<int>& v)
{
  double s = sum(v);

  std::vector<double> result(v.size());

  for (size_t i=0; i < v.size(); ++i)
    result[i] = v[i]/s;

  return result;
}

std::vector<double>
to_double_vector(const std::vector<int>& v)
{
  std::vector<double> result(v.size());

  for (size_t i=0; i < v.size(); ++i)
    result[i] = v[i];

  return result;
}

void
normalize_map(std::map<big_int, double>& v)
{
  double s = sum(v);

  typedef std::map<big_int, double>::iterator iterator;
  for (iterator i=v.begin(); i != v.end(); ++i)
    i->second = i->second/s;

  return;
}

bool
is_valid_string(const std::string& str, const character_to_values<bool>& is_valid)
{
  for (int i=0; i < str.length(); ++i) {
    if (not is_valid(str[i]))
      return false;
  }
  return true;
}

void
check_data(const std::vector<std::string>& seqs, const std::string& valid_chars)
{
  assert(seqs.size() != 0);
  character_to_values<bool> is_valid(valid_chars, true);
  for (int i=0; i < seqs.size(); ++i) {
    const std::string& line = seqs[i];
    assert(is_valid_string(line, is_valid));
  }

}

bool
is_nucleotide_string(const std::string& str)
{
  //std::string nucs = "ACGT";
  static character_to_values<bool> isnuc("ACGTU", true);
  for (int i=0; i < str.length(); ++i) {
    if (not isnuc(str[i]))
      return false;
  }
  return true;
}

std::string
reverse_complement(const std::string& s) 
{
  std::string t(s.size(), '-');
  int j=s.size()-1;
  for (int i=0; i < s.size(); ++i, --j)
    t[j] = complement(s[i]);
  return t;
}

std::string
reverse_complement_rna(const std::string& s) 
{
  std::string t(s.size(), '-');
  int j=s.size()-1;
  for (int i=0; i < s.size(); ++i, --j)
    t[j] = complement_rna(s[i]);
  return t;
}

std::string
reverse(const std::string& s) 
{
  std::string t(s.size(), ' ');
  int j=s.size()-1;
  for (int i=0; i < s.size(); ++i, --j)
    t[j] = s[i];
  return t;
}

bool
is_palindromic(const std::string& s)
{
  int len =s.length();
  if (len % 2 == 1)  // odd string cannot be palindromic
    return false;
  int middle = len/2;
  for (int i=0; i<middle; ++i) {
    if (s[i] != complement(s[len-i-1]))
      return false;
  } 
  return true;
}

int
palindromic_index(const std::string& s)
{
  return hamming_distance(s, reverse_complement(s));
}

matrix<double>
reverse(const matrix<double>& m)
{
  assert( m.get_rows() == 4 );
  
  int c = m.get_columns();
  matrix<double> result(4, c);
  for (int i = 0; i < 4; ++i)
    for (int j = 0; j < c; ++j) {
      result(i, j) = m(i, c-j-1);
    }

  return result;
}

// extend the matrix by k-1 columns of even distribution to both sides of the matrix
dmatrix
extend_matrix(const dmatrix& orig, int k)
{
  int width = orig.get_columns();
  int new_width = 2*(k-1) + width;
  dmatrix result(4, new_width);
  result.fill_with(0.25);
  result.inject(orig, 0, k-1);

  return result;
}


int
hamming_distance(const std::string& s, const std::string& t)
{
  assert(s.length() == t.length());

  int count = 0;
  for (int i=0; i < s.length(); ++i)
    if (s[i] != t[i])
      ++count;

//   if (count == 1)
//     std::cout << t << std::endl;

  return count;
}

std::vector<int>
iupac_hamming_mismatches(const std::string& s, const std::string& pattern)
{
  assert(s.length() == pattern.length());

  std::vector<int> result;
  for (int i=0; i < s.length(); ++i)
    if (not iupac_match(s[i], pattern[i]))
      result.push_back(i);

  return result;
}


int
iupac_hamming_dist(const std::string& str, const std::string& pattern, int max_hd)
{
  assert(str.length() == pattern.length());
  int hd=0;
  for (int i=0; i < str.length(); ++i) {
    if (not iupac_match(str[i], pattern[i])) {
      ++hd;
      if (hd > max_hd)
	return hd;
    }
  }
  
  return hd;
}

// finds the minimum Hamming distance between t and the substrings of s
int
min_hamming_distance(const std::string& s, const std::string& t)
{
  assert(s.length() >= t.length());
  int len = s.length();
  int k = t.length();
  int dist = k;
  for (int i = 0; i < len-k+1; ++i) {
    int temp = hamming_distance(s.substr(i,k), t);
    if (temp < dist)
      dist = temp; 
  }

  return dist;
}


std::string
join(const std::vector<std::string>& v, char c)
{
  int L = v[0].length();
  int n = v.size();
  std::string temp;
  temp.reserve((L+1)*n);
  for (int i=0; i < n-1; ++i) {
    temp.append(v[i]);
    temp.push_back(c);
  }
  temp.append(v[n-1]);

  return temp;
}

std::string
join(const std::vector<std::string>& v, const std::string& sep)
{
  int l = v.size();
  std::string temp;
  for (int i=0; i < l-1; ++i) {
    temp.append(v[i]);
    temp.append(sep);
  }
  temp.append(v.back());

  return temp;
}


// reverse complement of above catenation
// only separators are correctly placed
std::string
join_rev(const std::vector<std::string>& v, char c)
{
  int L = v[0].length();
  int lines = v.size();
  assert(lines>0);
  std::string temp;
  temp.reserve((L+1)*lines);
  for (int i=lines; i>1;) {
    --i;
    temp.append(reverse_complement(v[i]));
    temp.push_back(c);
  }
  temp.append(reverse_complement(v[0]));

  return temp;
}

std::vector<std::string>
split(const std::string& s, char c)
{
  size_t b=0;
  size_t e;
  std::vector<std::string> result;

  while (b != s.size() && (e=s.find(c, b)) != std::string::npos) {
    if (e>b)
      result.push_back(s.substr(b, e-b));
    b=e+1;
  }
  if (b != s.size())
    result.push_back(s.substr(b, s.size()-b));
    
  return result;
}


std::pair<int,int>
read_sequences(const std::string& filename, std::vector<std::string>& seqs, bool allow_iupac)
{
  assert(seqs.size() == 0);
  std::ifstream f;
  std::string line;
  int bad_lines = 0;
  f.open(filename.c_str(), std::ios_base::in);
  if (not f.is_open()) {
    std::cerr << "Couldn't open file " << filename << std::endl;
    exit(1);
  }

  while (getline(f, line)) {
    if (line.length() == 0) {
      ++bad_lines;
      continue;
    }
    while (not line.empty() and line.back() == '\r')
      line.pop_back();
	   
    if ((allow_iupac and is_iupac_string(line)) || is_nucleotide_string(line))
      seqs.push_back(line);
    else
      ++bad_lines;
  }
  f.close();

  error(seqs.size() == 0, "No valid sequences found! Exiting.\n");
  int lines = seqs.size();

  return std::make_pair(lines, bad_lines);
}

std::istream&
mygetline(std::istream& f, std::string& line)
{
  std::getline(f, line);
  while (not line.empty() and line.back() == '\r')
    line.pop_back();
  return f;
}

std::pair<int,int>
read_fastq_sequences(const std::string& filename, std::vector<std::string>& seqs, bool allow_iupac)
{
  assert(seqs.size() == 0);
  std::ifstream f;
  std::string line;
  int bad_lines = 0;
  f.open(filename.c_str(), std::ios_base::in);
  if (not f.is_open()) {
    std::cerr << "Couldn't open file " << filename << std::endl;
    exit(1);
  }

  std::string msg = "Does not appear to be a fastq file\n";
  std::string sequence;
  int block=0;
  int offset=0;
  bool e=false;
  while (mygetline(f, line)) {
    if (line.empty() or line[0] != '@') {
      e=true;
      offset=0;
      break;
    }

    if (not mygetline(f, sequence) or sequence.empty() or not is_iupac_string(sequence)) {
      e=true;
      offset=1;
      break;
    }
    
    if (not mygetline(f, line) or line.empty() or line[0] != '+') {
      e=true;
      offset=2;
      break;
    }

    if (not mygetline(f, line) or line.length() != sequence.length()) {
      e=true;
      offset=3;
      break;
    }

    ++block;
    if (sequence.length() == 0) {
      ++bad_lines;
      continue;
    }
	   
    if ((allow_iupac and is_iupac_string(sequence)) || is_nucleotide_string(sequence))
      seqs.push_back(sequence);
    else
      ++bad_lines;
  }
  f.close();
  if (e) {
    fprintf(stderr, "Read error on line %i\n", 4*block+offset+1);
    switch (offset) {
    case 0:
      fprintf(stderr, "Expected a line beginning with @\n");
      break;
    case 1:
      fprintf(stderr, "Expected a sequence of nucleotides\n");
      line = sequence;
      break;
    case 2:
      fprintf(stderr, "Expected a line beginning with +\n");
      break;
    case 3:
      fprintf(stderr, "Expected a sequence of quality codes, one for each nucleotide\n");
      break;
    }
    fprintf(stderr, "Instead got:\n%s\n", line.c_str());
    exit(1);
  }
  error(seqs.size() == 0, "No valid sequences found! Exiting.\n");
  int lines = seqs.size();

  return std::make_pair(lines, bad_lines);
}


std::pair<int,int>
read_fasta_sequences(const std::string& filename, std::vector<std::string>& seqs, bool allow_iupac)
{
  assert(seqs.size() == 0);
  std::ifstream f;
  std::string line;
  int bad_lines = 0;
  f.open(filename.c_str(), std::ios_base::in);
  if (not f.is_open()) {
    std::cerr << "Couldn't open file " << filename << std::endl;
    exit(1);
  }
  bool header_read=false;
  std::string current;
  while (getline(f, line)) {
    while (not line.empty() and line.back() == '\r')
      line.pop_back();

    if (not header_read) {                             // wait till we find the next header line
      if (line.length() > 0 and line[0] == '>')
	header_read = true;
      else
	++bad_lines;
      continue;
    }
	
     
    if (line.length() == 0) {
      ++bad_lines;
      current.clear();
      header_read=false;
      continue;
    }
    if (line[0] == '>') {
      if (current.length() == 0)
	++bad_lines;
      else {
	seqs.push_back(current);
      }
      current.clear();
      continue;
    }
    if ((allow_iupac and is_iupac_string(line)) || is_nucleotide_string(line))
      current += line;
    else {
      ++bad_lines;
      current.clear();
      header_read=false;
    }
  }
  f.close();
  if (current.length() != 0)
    seqs.push_back(current);

  error(seqs.size() == 0, "No valid sequences found! Exiting.\n");
  int lines = seqs.size();

  return std::make_pair(lines, bad_lines);
}

// x == base^3 * result[0] + base^2 * result[1] + base * result[2] + result[3]
std::vector<int>
decode_base(int base, int x)
{
  std::vector<int> result(4);
  result[3] = x % base;
  x /= base;
  result[2] = x % base;
  x /= base;
  result[1] = x % base;
  x /= base;
  result[0] = x;

  assert(result[0] < base);
  return result;
}



int
code_base(int base, int a, int b, int c, int d)
{
  assert(0 <= a && a < base);
  assert(0 <= b && b < base);
  assert(0 <= c && c < base);
  assert(0 <= d && c < base);
  return a*base*base*base + b*base*base + c*base + d;
}

int
code_base(int base, const std::vector<int>& v)
{
  assert(v.size() == 4);
  assert(0 <= v[0] && v[0] < base);
  assert(0 <= v[1] && v[1] < base);
  assert(0 <= v[2] && v[2] < base);
  assert(0 <= v[3] && v[3] < base);
  return v[0]*base*base*base + v[1]*base*base + v[2]*base + v[3];
}




std::vector<std::string>
integer_range(int begin, int end)
{
  std::vector<std::string> result;
  for (int i=begin; i < end; ++i)
    result.push_back(to_string("%i", i));

  return result;
}

// This is a helper function to be used inside gdb debugger
std::string&
SSS(const char* s)
{
  return *(new std::string(s));
}



std::vector<std::string>
get_n_neighbourhood(const std::string&seed, int n, bool expand)
{
  const int k = seed.length();
  assert(n >= 0);
  assert(n <= k);
  //  char nucs[] = "ACGT";

  std::vector<std::string> result;

  std::vector<std::string> complements(k);  // These are set complements, not nucleotide complements
  std::vector<int> bases(k);
  unsigned long long N_mask=0;          // bitmask for positions that contain 'N'. Those positions cannot contain an error
  for (int i=0; i < k; ++i) {
    complements[i] = complement_set(seed[i]);
    bases[i]=complements[i].length() - 1;
    if (seed[i]=='N')
      N_mask |=  (1 << (k-1-i));
  }

  result.push_back(seed);                   // hamming distance 0
  for (int error=1; error <= n; ++error) { // handles hamming distances 1 <= hd <= n
    // bitvector c has 1-bit for each member of the subset, rightmost bit is bit number k-1
    unsigned long long c = (1ull<<error)-1;  // initially rightmost 'error' bits are 1
    // iterate through all subsets c of {0, ..., k-1} that have size 'error'
    while (c < (1ull<<k)) {   // Superset has only k elements
      assert(__builtin_popcountll(c) == error);
      if ((c & N_mask) == 0)  { // subset positions don't contain 'N'
	std::vector<int> P;  // positions that contain error
	int number_of_combinations = 1;
	std::string temp = seed;
	for (int pos=0; pos < k; ++pos) {
	  if ((c & (1ull << (k-1-pos))) != 0) {   // pos belongs to set c
	    P.push_back(pos);
	    temp[pos] = complements[pos][0];  // initialize to first string that has mismatches at these positions
	    number_of_combinations *= complements[pos].length();
	  }
	}
	  
	std::vector<int> y(error, 0);
	y[error-1]=-1;  // Initialize
	for (int r=0; r < number_of_combinations; ++r) {
      
	  int i;
	  for (i=error-1; y[i] == bases[P[i]]; --i) {
	    y[i]=0;
	    temp[P[i]] = complements[P[i]][y[i]];
	  }
	  y[i]++;
	  temp[P[i]] = complements[P[i]][y[i]];

	  result.push_back(temp);
	}  // end for r
	
      }    // end if subset doesn't contain 'N'

      unsigned long long a = c&-c;
      unsigned long long b = c+a;   // update bitvector c. This is "Gosper's hack"
      c = (c^b)/4/a|b;
      
    } // end foreach subset c
  }  // end for error

  if (not is_nucleotide_string(seed) and is_iupac_string(seed) and expand) {                     // expand iupac strings into nucleic strings
    std::set<std::string> expanded;                           // Tämä on tod. näk. turha, koska eri patternien laajennokset ovat pistevieraita
    BOOST_FOREACH(std::string pattern, result) {
      sequences_of_iupac_string x(pattern);
      int size = x.number_of_strings();
      for (int i=0; i < size; ++i)
	expanded.insert(x.next());
      
    }
    return std::vector<std::string>(expanded.begin(), expanded.end());
  }
  else
    return result;
}
