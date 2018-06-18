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
#ifndef MATRIX_TOOLS_HPP
#define MATRIX_TOOLS_HPP

#include "matrix.hpp"

#include <string>
#include <limits>

template <typename T>
class accumulate
{
public:
  accumulate() : sum(0) {}

  T 
  operator()(T v)
  {
    sum += v;
    return sum;
  }

  T get() { return sum; }
private:
  T sum;
};

template <typename T>
class gather_maximum
{
public:
  gather_maximum(T t) : max(t) {}

  void 
  operator()(T v)
  {
    max = std::max(max, v);
  }

  T get() const { return max; }

private:
  T max;
};

template <typename T>
class gather_minimum
{
public:
  gather_minimum(T t) : min(t) {}

  void 
  operator()(T v)
  {
    min = std::min(min, v);
  }

  T get() const { return min; }

private:
  T min;
};

template <typename From, typename To>
matrix<To>
convert(const matrix<From>& f)
{
  matrix<To> result(f.dim());
  for (int i=0; i < f.get_rows(); ++i)
    for (int j=0; j < f.get_columns(); ++j)
      result(i, j) = f(i, j);

  return result;
}

template <typename T, typename F>
matrix<T>
log2(const matrix<F>& orig)
{
  matrix<T> result(orig.dim());
  int rows=orig.get_rows();
  int cols=orig.get_columns();
  for (int i=0; i < rows; ++i)
    for (int j=0; j < cols; ++j)
      result(i, j) = log2l(orig(i, j));
  //      result(i, j) = orig(i,j)==0.0 ? 0.0 : log2l(orig(i, j));

  return result;
}

// This version assumes that on the first column only first four values are non-zero.
template <typename T, typename F>
matrix<T>
log2_special(const matrix<F>& orig)
{
  matrix<T> result(orig.dim());
  int rows=orig.get_rows();
  int cols=orig.get_columns();
  for (int i=0; i < 4; ++i)
    for (int j=0; j < cols; ++j)
      result(i, j) = log2l(orig(i, j));
  //result(i, j) = orig(i,j)==0.0 ? 0.0 : log2l(orig(i, j));
  for (int i=4; i < rows; ++i)
    for (int j=1; j < cols; ++j)
      result(i, j) = log2l(orig(i, j));
  //result(i, j) = orig(i,j)==0.0 ? 0.0 : log2l(orig(i, j));

  return result;
}


template <typename T>
T
sum(const matrix<T>& m)
{
  accumulate<T> b;
  b=m.iterate(b);
  
  return b.get();
}

template <typename T>
T
max(const matrix<T>& m)
{
  gather_maximum<T> b(std::numeric_limits<T>::min());
  b=m.iterate(b);
  
  return b.get();
}

template <typename T>
T
min(const matrix<T>& m)
{
  gather_minimum<T> b(std::numeric_limits<T>::max());
  b=m.iterate(b);
  
  return b.get();
}

void
normalize_matrix_rows(matrix<double>& m);

// void
// normalize_matrix_columns(matrix<double>& m);

template <typename T>
void
normalize_matrix_columns(matrix<T>& m)
{
  // sums of cols should be 1
  for (int i=0; i < m.get_columns(); ++i) {
    T sum = 0;
    for (int j=0; j < m.get_rows(); ++j)
      sum += m(j,i);
    
    for (int j=0; j < m.get_rows(); ++j)
      m(j,i) /= sum;
  }
}


void
normalize_whole_matrix(matrix<double>& m);


dmatrix
normalize_matrix_columns_copy(matrix<double> m);


dmatrix
normalize_matrix_rows_copy(matrix<double> m);

// bool
// is_column_stochastic_matrix(const matrix<double>& m);

// bool
// is_row_stochastic_matrix(const matrix<double>& m);

template <typename T>
bool
is_stochastic_matrix(const matrix<T>& m)
{
  //const double delta = 0.000001;
  const T delta = 0.00001;
  for (int i=0; i < m.get_rows(); ++i) {
    T sum = 0;
    for (int j=0; j < m.get_columns(); ++j)
      sum += m(i,j);
    if (fabs(1.0 - sum) >= delta)  // isn't close enough to 1
      return false;
  }
  return true;
}

template <typename T>
bool
is_column_stochastic_matrix(const matrix<T>& m)
{
  return is_stochastic_matrix(transpose(m));
}

template <typename T>
bool
is_row_stochastic_matrix(const matrix<T>& m)
{
  return is_stochastic_matrix(m);
}


bool
is_palindromic_matrix(const matrix<double>& m);

#include <stdexcept>

class fileformat_error : public std::runtime_error
{
public:

  fileformat_error(const std::string& msg) : runtime_error(msg) {}

  //  const char* what() const throw() { return message.c_str(); }
 
  virtual ~fileformat_error() throw() {}

private:
  std::string message;
};

matrix<double>
read_matrix(FILE* fp);


template<typename T>
void
write_matrix(FILE* fp, const matrix<T>& m, const std::string& tag, 
	     std::string format = "", bool dimensions=true)
{
  fprintf(fp, "%s", tag.c_str());

  int r = m.get_rows();
  int c = m.get_columns();
  /*
  if (dimensions)
    fprintf(fp, "%ix%i\n", r, c);
  */
  
  if (format != "") 
    m.print3(fp, format);
  else
    m.print3(fp, "%10lf", "\t");

}

void
write_matrix_file(const std::string& matrixfile, const dmatrix& M, std::string format= "%.6f");

dmatrix
read_matrix_file(const std::string& matrixfile);

dmatrix
matrix_sum(const dmatrix& m1, const dmatrix& m2, int d);

dmatrix
matrix_product(const dmatrix& m1, const dmatrix& m2, int d);

#endif //MATRIX_TOOLS_HPP
