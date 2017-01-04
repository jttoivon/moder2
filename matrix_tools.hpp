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

void
normalize_matrix_columns(matrix<double>& m);

void
normalize_whole_matrix(matrix<double>& m);


dmatrix
normalize_matrix_columns_copy(matrix<double> m);


dmatrix
normalize_matrix_rows_copy(matrix<double> m);

bool
is_column_stochastic_matrix(const matrix<double>& m);

bool
is_row_stochastic_matrix(const matrix<double>& m);

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


void
write_matrix(FILE* fp, const matrix<double>& m, const std::string& tag, 
	     std::string format = "", bool dimensions=true);

void
write_matrix_file(const std::string& matrixfile, const dmatrix& M, std::string format= "%.6f");

dmatrix
read_matrix_file(const std::string& matrixfile);

dmatrix
matrix_sum(const dmatrix& m1, const dmatrix& m2, int d);

dmatrix
matrix_product(const dmatrix& m1, const dmatrix& m2, int d);

#endif //MATRIX_TOOLS_HPP
