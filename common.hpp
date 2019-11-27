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
#ifndef COMMON_HPP
#define COMMON_HPP

#include "matrix.hpp"
#include "data.hpp"
#include "cob.hpp"
#include "orientation.hpp"

#include <cfloat>
#include <vector>
#include <map>
#include <string>
#include <cstring>
#include <boost/multi_array.hpp>
#include <boost/unordered_map.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/shared_ptr.hpp>

extern std::string header;

typedef std::vector<int> ivector;

double
cut(double x);


void
error(bool b, const std::string& message);

std::string
mybasename(const std::string& path);

int
bits_in_range(const std::bitset<64>& occ, int begin, int end);


std::string
print_bitvector(unsigned int x, int number_of_bits=sizeof(unsigned int)*8);


const char*
yesno(bool b);


template <typename T>
std::string
m(std::string s, T t)
{
  std::ostringstream ss;
  ss << s << ", default: " << t;
  return ss.str();
}

void
print_command_line(int argc, char* argv[]);

void
print_command_line(int argc, char* argv[], FILE* f);

std::vector<int>
count_frequencies(const std::string& s);


std::string
itoa(int i);

int
atoi(const std::string& s);

double
atof(const std::string& s);


std::string
to_string(std::string format, ...);

void
underline(std::string format, ...);

extern matrix<double> positional_background;

template <typename Iter, typename Key>
void
key_sort(Iter begin, Iter end, Key key);

template <typename Value, typename Key>
class key_sort_functor
{
public:

  key_sort_functor(Key k) : k_(k) {} 

  bool
  operator()(Value x, Value y) const
  {
    return k_(x) < k_(y); 
  }

private:
  Key k_;
};

template <typename Iter, typename Key>
void
key_sort(Iter begin, Iter end, Key key)
{
  typedef typename Iter::value_type Value;
  std::sort(begin, end, key_sort_functor<Value, Key>(key));
}

class division_functor
{
public:
  
  division_functor(double d_) : d(d_) {}

  double
  operator()(double x) const {
    return x/d;
  }

private:
  double d;
};

class multiplication_functor
{
public:
  
  multiplication_functor(double d_) : d(d_) {}

  double
  operator()(double x) const {
    return x*d;
  }

private:
  double d;
};

std::string
reverse_complement(const std::string& s);

std::string
reverse_complement_rna(const std::string& s);

std::string
reverse(const std::string& s);

bool
is_palindromic(const std::string& s);

int
palindromic_index(const std::string& s);

bool
is_nucleotide_string(const std::string& str);



// Reflects over both diagonals, that is rotate 180 degrees. This is NOT the transpose of the matrix
template<typename T>
matrix<T>
reverse_complement(const matrix<T>& m)
{
  assert( m.get_rows() == 4 );
  
  int c = m.get_columns();
  matrix<T> result(4, c);
  for (int i = 0; i < 4; ++i)
    for (int j = 0; j < c; ++j) {
      result(i, j) = m(4-i-1, c-j-1);
    }

  return result;
}


matrix<double>
reverse(const matrix<double>& m);

dmatrix
extend_matrix(const dmatrix& orig, int k);



//char complement(char c);

/*
// convert nucleotides ACGT to 0123
// Other characters are undefined
int
to_int(char c);

void
initialize_to_int();
*/


template <typename T>
boost::tuple<int,int,int,int>
get_ranges(const boost::multi_array<T, 2>& a)
{
  int row_begin, row_end, col_begin, col_end;
  row_begin = a.index_bases()[0];
  col_begin = a.index_bases()[1];
  row_end = row_begin + a.shape()[0];
  col_end = col_begin + a.shape()[1];
  return boost::make_tuple(row_begin, row_end, col_begin, col_end);
}




void
normalize_vector(std::vector<double>& v);

// std::vector<double>
// normalize_vector_copy(const std::vector<double>& v);

std::vector<double>
normalize_vector_copy(const std::vector<int>& v);

std::vector<double>
to_double_vector(const std::vector<int>& v);

void
normalize_map(std::map<big_int, double>& v);




std::pair<int,int>
read_sequences(const std::string& filename, std::vector<std::string>& seqs, bool allow_iupac=false);

std::pair<int,int>
read_fasta_sequences(const std::string& filename, std::vector<std::string>& seqs, bool allow_iupac=false);

std::pair<int,int>
read_fastq_sequences(const std::string& filename, std::vector<std::string>& seqs, bool allow_iupac=false);

std::string
join(const std::vector<std::string>& v, char c);

std::string
join(const std::vector<std::string>& v, const std::string& sep);

std::string
join_rev(const std::vector<std::string>& v, char c);

std::vector<std::string>
split(const std::string& s, char c);

template <typename T,
	  typename F,
	  typename R>
std::vector<R>
transform(const std::vector<T>& orig, F f)
{
  std::vector<R> result;
  std::transform(orig.begin(), orig.end(), std::back_inserter(result), f);
  return result;
}


class character_translation
{
public:

  character_translation(const std::string& from, const std::string& to) 
    : table(256)
  {
    assert(from.length() == to.length());

    for (int i=0; i < 256; ++i)
      table[i]=i;

    for (int i=0; i < from.length(); ++i)
      table[(unsigned char)from[i]] = to[i];
  }

  char
  operator()(char c) const
  {
    return table[(unsigned char) c];
  }

private:
  std::vector<char> table;
};

template <typename T>
class character_to_values
{
public:

  character_to_values(const std::string& from, const std::vector<T>& to, T default_value=T()) 
    : table(256)
  {
    assert(from.length() == to.size());

    for (int i=0; i < 256; ++i)
      table[i] = default_value;

    for (int i=0; i < from.length(); ++i)
      table[(unsigned char)from[i]] = to[i];
  }

  character_to_values(const std::string& from, T to, T default_value=T()) 
    : table(256)
  {

    for (int i=0; i < 256; ++i)
      table[i] = default_value;

    for (int i=0; i < from.length(); ++i)
      table[(unsigned char)from[i]] = to;
  }

  T
  operator()(char c) const
  {
    return table[(unsigned char) c];
  }

private:
  std::vector<T> table;
};

int
hamming_distance(const std::string& s, const std::string& t);

std::vector<int>
iupac_hamming_mismatches(const std::string& s, const std::string& pattern);

int
iupac_hamming_dist(const std::string& str, const std::string& pattern, int max_hd);

template <typename T>
std::string
iupac_string_giving_max_probability(const matrix<T>& dm, bool use_rna)
{

  const char* nucs = use_rna ? "ACGU" : "ACGT";
  int k = dm.get_columns();
  std::string result(k, '-');

  for (int i=0; i < k; ++i) {
    const std::vector<T>& c = dm.column(i);
    std::vector<int> v = {0,1,2,3};  // sort indexes for the column
    std::sort(v.begin(), v.end(), [c,v](int a, int b) { return c[v[a]] >= c[v[b]]; });
    if (c[v[0]] > 0.5 and c[v[0]] >= 2*c[v[1]])
      result[i] = nucs[v[0]];   // single nucleotide
    else if (c[v[0]]+c[v[1]] > 0.75)  // two nucleotides
      result[i] = iupac_class.bits_to_char(iupac_class.char_to_bits(nucs[v[0]]) |
					   iupac_class.char_to_bits(nucs[v[1]]));
    else if (c[v[3]] < 0.01)    // three nucleotides
      result[i] = iupac_class.bits_to_char(iupac_class.char_to_bits(nucs[v[0]]) |
					   iupac_class.char_to_bits(nucs[v[1]]) |
					   iupac_class.char_to_bits(nucs[v[2]]));
    else
      result[i] = 'N';
  }
  return result;
}

// finds the minimum Hamming distance between t and the substrings of s
int
min_hamming_distance(const std::string& s, const std::string& t);


void
check_data(const std::vector<std::string>& seqs, const std::string& valid_chars="ACGT");

// recalculate matrix m using frequencies of motif nucleotides 
// in subset [begin,end] of sequences
template <typename Iter>
matrix<double>
calculate_frequencies_in_subset(const std::vector<std::string>& sequences, Iter begin, Iter end, 
				const std::vector<int>& alignments, 
				const std::vector<int>& direction, 
				int k)
{
  matrix<double> m_new(4,k);
  m_new.fill_with(0);

  for (Iter i=begin;  i < end; ++i) {
    const std::string& line = sequences[*i];
    int n = line.length();
    assert(n>=k);
    int l = alignments[*i];
    std::string tmp = line.substr(l,k);
    if (direction[*i]==-1)
      tmp=reverse_complement(tmp);
    for (int j=0; j < k; ++j) {
      ++m_new(to_int(tmp[j]), j);
    }

  }

  // add pseudo counts
//   int pseudo_count=0;
//   const int pad = 1;
//   for (int i=0; i < 4; ++i)
//     for (int j=0; j < k; ++j)
//       if (m_new(i,j)==0) {
// 	m_new(i,j)=pad;
// 	pseudo_count += pad;
//       }

  //printf("pseudocount =%i    charcount=%i\n", pseudo_count, char_count);

  return m_new;

}

template <typename T, size_t N>
T
sum(T (&array)[N])
{
  T sum=0;
  for (size_t i=0; i < N; ++i)
    sum += array[i];

  return sum;
}

template <typename T>
T
sum(const std::vector<T>& array)
{
  T sum=0;
  size_t N = array.size();
  for (size_t i=0; i < N; ++i)
    sum += array[i];

  return sum;
}

template <typename T>
T
sum(const boost::multi_array<T, 2>& a) 
{
  int row_begin, row_end;
  int col_begin, col_end;
  boost::tie(row_begin, row_end, col_begin, col_end) = get_ranges(a);

  T sum = 0.0;
  for (int row=row_begin; row < row_end; ++row) {
    for (int col=col_begin; col < col_end; ++col) {
      sum += a[row][col];
    }
  }
  return sum;
}

/*

template <typename T>
T
sum(const matrix<T>& a) 
{
  int row_end;
  int col_end;
  boost::tie(row_end, col_end) = a.dim();

  T sum = 0.0;
  for (int row=0; row < row_end; ++row) {
    for (int col=0; col < col_end; ++col) {
      sum += a(row,col);
    }
  }
  return sum;
}
*/

template <typename K, typename T>
T
sum(const std::map<K, T>& m)
{
  T sum=0;

  typedef typename std::map<K,T>::const_iterator iterator;
  for (iterator i=m.begin(); i != m.end(); ++i)
    sum += i->second;

  return sum;
}

template <typename K, typename T, typename H>
T
sum(const boost::unordered_map<K, T, H>& m)
{
  T sum=0;

  typedef typename boost::unordered_map<K,T,H>::const_iterator iterator;
  for (iterator i=m.begin(); i != m.end(); ++i)
    sum += i->second;

  return sum;
}


template <typename T>
std::vector<T>
normalize_vector_copy(const std::vector<T>& v)
{
  T s = sum(v);

  std::vector<T> result(v.size());

  for (size_t i=0; i < v.size(); ++i)
    result[i] = v[i]/s;

  return result;
}


template <typename T>
std::string
print_vector(const std::vector<T>& v)
{
  if (v.size() == 0)
    return std::string("[]");
  std::ostringstream result;
  result << "[";
  for (int i=0; i < v.size()-1; ++i)
    result << v[i] << ", ";
  result << v[v.size()-1] << "]";

  return result.str();
}

template <typename T>
std::string
print_vector(const std::vector<boost::shared_ptr<T> >& v)
{
  if (v.size() == 0)
    return std::string("[]");
  std::ostringstream result;
  result << "[";
  for (int i=0; i < v.size()-1; ++i)
    result << *v[i] << ", ";
  result << *v[v.size()-1] << "]";

  return result.str();
}

template <typename T>
std::string
print_vector(const std::vector<T>& v, const std::string& sep, int precision)
{
  if (v.size() == 0)
    return std::string("");
  std::ostringstream result;
  result.precision(precision);
  for (int i=0; i < v.size()-1; ++i)
    result << v[i] << sep;
  result << v[v.size()-1];

  return result.str();
}


template <>
std::string
print_vector<bool>(const std::vector<bool>& v);

template <typename T>
T
max_element(const std::vector<T>& array)
{
  T max=-std::numeric_limits<T>::max();
  size_t N = array.size();
  for (size_t i=0; i < N; ++i)
    max = std::max(max, array[i]);

  return max;
}

template <typename T>
T
max_element(const boost::multi_array<T, 2>& a)
{
  int row_begin, row_end, col_begin, col_end;
  boost::tie(row_begin, row_end, col_begin, col_end) = get_ranges(a);
  T max = -std::numeric_limits<T>::max();
  for (int r = row_begin; r < row_end; ++r)
    for (int c = col_begin; c < col_end; ++c)
      max = std::max(a[r][c], max);
  return max;
}


template <typename T>
T
min_element(const std::vector<T>& array)
{
  T min=array[0];
  size_t N = array.size();
  for (size_t i=1; i < N; ++i)
    min = std::min(min, array[i]);

  return min;
}

/*
template <typename K, typename V>
V
max_element(const my_unordered_map<K,V>& hash)
{
  return std::max_element(hash.begin(), hash.end())->second;
}
*/

template <typename K, typename V, typename H>
V
max_element(const boost::unordered_map<K, V, H>& m)
{
  if (m.size() == 0)
    throw "sika";

  V max=m.begin()->second;

  typedef typename boost::unordered_map<K,V, H>::const_iterator iterator;
  for (iterator i=m.begin(); i != m.end(); ++i)
    max = std::max(max, i->second);

  return max;
}


template <typename K, typename V, typename H>
K
arg_max(const boost::unordered_map<K, V, H>& m)
{
  if (m.size() == 0)
    throw "Cannot take maximum of an empty container";

  //  K arg_max=m.begin()->first;
  //  V max=m.begin()->second;
  //  std::cout << "Tassa" << std::endl;
  //  std::cout << m.begin()->first << std::endl;
  //  std::cout << m.begin()->second << std::endl;
  K arg_max=m.begin()->first;
  V max=m.begin()->second;

  typedef typename boost::unordered_map<K,V, H>::const_iterator iterator;
  for (iterator i=m.begin(); i != m.end(); ++i)
    if (i->second > max) {
      arg_max = i->first;
      max = i->second;
    }

  return arg_max;
}


/*
template <typename K, typename V>
V
arg_max(const boost::unordered_map<K,V>& hash)
{
  return std::max_element(hash.begin(), hash.end())->first;
}
*/

template <typename T>
typename std::vector<T>::size_type
arg_max(const std::vector<T>& array)
{
  assert(array.size() > 0);
  size_t maxarg=0;
  T max=array[maxarg];
  size_t N = array.size();
  for (size_t i=1; i < N; ++i) {
    if (array[i]>max) {
      max=array[i];
      maxarg=i;
    }
  }
  return maxarg;
}

template <typename T>
double
add_pseudo_counts(matrix<T>& m, const double pad = 1.0)
{
  int rows = m.get_rows();
  int columns = m.get_columns();
  double count = 0;
  for (int i=0; i < rows; ++i)
    for (int j=0; j < columns; ++j) {
	m(i,j)+=pad;
	count+=m(i,j);
    }
  return count;
}

template <typename T>
class prior
{
public:

  prior() : pseudo_counts(4,0) {}
  
  void
  use_add_one(double p = 1.0) 
  {
    for (int i=0; i < 4; ++i)
      pseudo_counts[i] = p;
  }

  void
  use_dirichlet(double b, const std::vector<double>& q)
  {
    assert( q.size() == 4 );
    for (int i=0; i < 4; ++i)
      pseudo_counts[i] = q[i] * b;
  }

  void
  add(matrix<T>& m) const
  {
    int rows = m.get_rows();
    int columns = m.get_columns();
    assert(rows == 4);

    for (int i=0; i < rows; ++i)
      for (int j=0; j < columns; ++j)
	m(i,j) += pseudo_counts[i];
  }

  void
  add(std::vector<T>& v) const 
  {
    int size = v.size();
    
    for (int i=0; i < size; ++i)
      v[i] += pseudo_counts[i];
  }

  std::vector<double>
  get() const
  {
    return pseudo_counts;
  }

private:
  std::vector<double> pseudo_counts;
};

template <typename T>
class dinucleotide_prior
{
public:

  dinucleotide_prior() : pseudo_counts(16,0) {}
  
  void
  use_add_one(double p = 1.0) 
  {
    for (int i=0; i < 16; ++i)
      pseudo_counts[i] = p;
  }

  void
  use_dirichlet(double b, const std::vector<double>& q)
  {
    background_probabilities = q;
    assert( q.size() == 4 );
    for (int i=0; i < 4; ++i)
      for (int j=0; j < 4; ++j)
        pseudo_counts[(i<<2) + j] = q[i] * q[j] * b;
  }

  void
  add(matrix<T>& m) const
  {
    int rows = m.get_rows();
    int columns = m.get_columns();
    assert(rows == 16);

    for (int i=0; i < rows; ++i)
      for (int j=0; j < columns; ++j)
        m(i,j) += pseudo_counts[i];
  }

  void
  add(std::vector<T>& v) const
  {
    int size = v.size();
    assert(size == 16);
    
    for (int i=0; i < size; ++i)
      v[i] += pseudo_counts[i];
  }

  std::vector<double>
  get() const
  {
    return pseudo_counts;
  }
  std::vector<double> background_probabilities;

private:
  std::vector<double> pseudo_counts;
};


// x == base^2 * result[0] + base * result[1] + result[2]
std::vector<int>
decode_base(int base, int x);


int
code_base(int base, int a, int b, int c, int d);


int
code_base(int base, const std::vector<int>& v);


namespace {

  template <typename T>
  void
  print_cell(FILE* fp, const std::string& format, const T& cell)
  {
    fprintf(fp, format.c_str(), cell);
  }

  template <>
  void
  print_cell(FILE* fp, const std::string& format, const std::string& cell)
  {
    fprintf(fp, format.c_str(), cell.c_str());
  }

}

/*
template <typename T, unsigned long N, unsigned long I>
typename boost::multi_array<T, N>::extent_gen::template gen_type<N>::type
xget_shape_helper(const boost::multi_array<T, N>& a, typename boost::multi_array<T, N>::extent_gen::template gen_type<I>::type ext)
{
  typedef typename boost::multi_array<T, N>::extent_range range;
  int row_begin, row_end, col_begin, col_end;
  row_begin = a.index_bases()[I];
  row_end = row_begin + a.shape()[I];
  return xget_shape_helper<T, N, I+1>(a, ext[range(row_begin, row_end)]);
}

template <typename T, unsigned long N>
typename boost::multi_array<T, N>::extent_gen::template gen_type<N>::type
xget_shape_helper<T,N,N>(const boost::multi_array<T, N>& a, typename boost::multi_array<T, N>::extent_gen::template gen_type<N>::type ext)
{
  return ext;
}


template <typename T, unsigned long N>
typename boost::multi_array<T, N>::extent_gen::template gen_type<N>::type
xget_shape(const boost::multi_array<T, N>& a)
{
  typedef typename boost::multi_array<T, N>::extent_range range;
  int row_begin, row_end, col_begin, col_end;
  row_begin = a.index_bases()[0];
  row_end = row_begin + a.shape()[0];
  return xget_shape_helper<T, N, 1>(a, boost::extents[range(row_begin, row_end)]);
}
*/

template <typename T>
typename boost::multi_array<T, 2>::extent_gen::template gen_type<2>::type
get_shape(const boost::multi_array<T, 2>& a)
{
  typedef typename boost::multi_array<T, 1>::extent_range range;
  int row_begin, row_end, col_begin, col_end;
  row_begin = a.index_bases()[0];
  row_end = row_begin + a.shape()[0];
  col_begin = a.index_bases()[1];
  col_end = col_begin + a.shape()[1];
  return boost::extents[range(row_begin, row_end)][range(col_begin,col_end)];
}

template <typename T>
typename boost::multi_array<T, 1>::extent_gen::template gen_type<1>::type
get_shape(const boost::multi_array<T, 1>& a)
{
  typedef typename boost::multi_array<T, 1>::extent_range range;
  int row_begin, row_end;
  row_begin = a.index_bases()[0];
  row_end = row_begin + a.shape()[0];
  return boost::extents[range(row_begin, row_end)];
}

template <typename T>
void
fill_with(boost::multi_array<T, 2>& a, const T& value)
{
  int row_begin, row_end, col_begin, col_end;
  boost::tie(row_begin, row_end, col_begin, col_end) = get_ranges(a);

  for (int r = row_begin; r < row_end; ++r)
    for (int c = col_begin; c < col_end; ++c)
      a[r][c] = value;
}


std::vector<std::string>
integer_range(int begin, int end);

template <typename T>
void
print_matrix(FILE* fp, const matrix<T>& a, 
             const std::vector<std::string>& row_headers = std::vector<std::string>(), 
             const std::vector<std::string>& col_headers = std::vector<std::string>(),
             const std::string& format = "%.2f")
{
  //  int row_begin;
  int row_end;
  //int col_begin;
  int col_end;
  //col_begin = 0;
  boost::tie(row_end, col_end) = a.dim();

  bool use_row_headers = row_headers.size() == row_end;
  bool use_col_headers = col_headers.size() == col_end;
  if (use_col_headers) {
    if (use_row_headers)
      fprintf(fp, "\t");
    for (int col=0; col < col_end-1; ++col)
      fprintf(fp, "%s\t", col_headers[col].c_str());
    fprintf(fp, "%s\n", col_headers[col_end-1].c_str());
  }

  
  for (int row=0; row < row_end; ++row) {
    if (use_row_headers)
      fprintf(fp, "%s\t", row_headers[row].c_str());
    for (int col=0; col < col_end-1; ++col) {
      print_cell(fp, format, a(row, col));
      fprintf(fp, "\t");
    }
    print_cell(fp, format, a(row, col_end-1));
    fprintf(fp, "\n");
  }
}



template <typename T>
void
print_array(FILE* fp, const boost::multi_array<T, 2>& a, 
	    const std::vector<std::string>& row_headers = std::vector<std::string>(), 
	    const std::vector<std::string>& col_headers = std::vector<std::string>(),
	    const std::string& format = "%.2f")
{
  int row_begin, row_end;
  int col_begin, col_end;
  boost::tie(row_begin, row_end, col_begin, col_end) = get_ranges(a);

  bool use_row_headers = row_headers.size() == (row_end - row_begin);
  bool use_col_headers = col_headers.size() == (col_end - col_begin);
  if (use_col_headers) {
    if (use_row_headers)
      fprintf(fp, "\t");
    for (int col=col_begin; col < col_end-1; ++col)
      fprintf(fp, "%s\t", col_headers[col-col_begin].c_str());
    fprintf(fp, "%s\n", col_headers[col_end-col_begin-1].c_str());
  }

  
  for (int row=row_begin; row < row_end; ++row) {
    if (use_row_headers)
      fprintf(fp, "%s\t", row_headers[row-row_begin].c_str());
    for (int col=col_begin; col < col_end-1; ++col) {
      print_cell(fp, format, a[row][col]);
      //      printf(format.c_str(), a[row][col]);
      fprintf(fp, "\t");
    }
    print_cell(fp, format, a[row][col_end-1]);
    //    printf(format.c_str(), a[row][col_end-1]);
    fprintf(fp, "\n");
  }
}

template <typename T>
boost::multi_array<T, 2>
read_array(FILE* fp)
{
  typedef typename boost::multi_array<T, 2> cob_type;

  typedef typename cob_type::extent_range range;
  typename cob_type::extent_gen extents;


  char* line = NULL;
  size_t dummy;
  int ret_val = getline(&line, &dummy, fp);
  assert(ret_val > 0);
  std::vector<std::string> parts = split(line, '\t');
  //assert(parts[0].length() == 0);
  free(line);
  int dmin = atoi(parts[0]);
  int dmax = atoi(parts[parts.size()-1]);
  int n=0;
  std::vector<std::string> lines;
  line = NULL;
  ret_val = getline(&line, &dummy, fp);
  while (ret_val > 0 and n < 4) {
    lines.push_back(line);
    free(line);
    line = NULL;
    ret_val = getline(&line, &dummy, fp);
    ++n;
  }
  free(line);

  assert(n == 3 or n == 4);
  cob_type result(extents[n][range(dmin, dmax+1)]);


  for (int i=0; i < n; ++i) {
    parts = split(lines[i], '\t');
    assert(orients[i] == parts[0]);
    for (int j=0; j < (dmax-dmin+1); ++j) {
      result[i][dmin+j] = atof(parts[j+1]);
    }
  }
  return result;
}

template <typename T>
boost::multi_array<T, 2>
read_array_file(const std::string& arrayfile)
{
  FILE* fp = fopen(arrayfile.c_str(), "r");
  if (fp == NULL) {
    fprintf(stderr, "Couldn't open file %s\n", arrayfile.c_str());
    exit(1);
  }
  boost::multi_array<T, 2> M = read_array<T>(fp);
  fclose(fp);
  return M;
}


template <typename T>
void
print_cob(FILE* fp, const boost::multi_array<T, 2>& a, const std::string& title = "", const std::string& format = "%.2f")
{
  int row_begin, row_end, col_begin, col_end;
  boost::tie(row_begin, row_end, col_begin, col_end) = get_ranges(a);  

  assert(row_begin == 0 and (row_end == 3 or row_end == 4 or row_end == 1 or row_end == 2));
  std::vector<std::string> col_headers = integer_range(col_begin, col_end);
  std::vector<std::string> row_headers;
  row_headers.push_back("HT");
  if (row_end == 2)
    row_headers.push_back("TH");
  else if (row_end >= 3) {
    row_headers.push_back("HH");
    row_headers.push_back("TT");
    if (row_end == 4)
      row_headers.push_back("TH");
  }

  if (title.length() > 0)
    fprintf(fp, "%s", title.c_str());
  print_array(fp, a, row_headers, col_headers, format);
}

template <typename T>
void
write_cob_file(const std::string& filename, const boost::multi_array<T, 2>& a)
{
  FILE* fp = fopen(filename.c_str(), "w");
  if (fp == NULL) {
    fprintf(stderr, "Couldn't open file %s for writing\n", filename.c_str());
    exit(1);
  }
  print_cob(fp, a, "", "%.5f");
  fclose(fp);
}

std::string&
SSS(const char* s);

std::vector<std::string>
get_n_neighbourhood(const std::string&seed, int n);

#endif // COMMON_HPP

