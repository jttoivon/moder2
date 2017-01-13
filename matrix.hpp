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
#ifndef MATRIX_HPP
#define MATRIX_HPP

#include <cassert>
#include <cstdio>
#include <cmath>
#include <algorithm>
#include <iomanip>
#include <iostream>
#include <vector>

template <typename T>
class matrix
{
public:

  matrix(int n=0, int m=0) : rows(n), columns(m), p(n*m) {}

  matrix(const std::pair<int,int>& pair) : rows(pair.first), columns(pair.second), p(pair.first * pair.second) {}

  matrix(const matrix& other) {
    rows=other.rows;
    columns=other.columns;
    p=other.p;
  }

  void
  assign(int n, int m, T t = T()) 
  {
    rows=n; 
    columns=m;
    p.assign(n*m, t);
  }

  matrix&
  operator=(const matrix& other) {
    rows=other.rows;
    columns=other.columns;
    p=other.p;
    return *this;
  }

  // inject another matrix to coordinates (n,m) of *this matrix
  void
  inject(const matrix& other, int n, int m) {
    int orows=other.rows;
    int ocolumns=other.columns;
    assert(n >= 0 && m >= 0);
    assert(n+orows <= rows && m+ocolumns <= columns);

    for (int i=0;i<orows;++i) {
      for (int j=0;j<ocolumns;++j) 
	(*this)(n+i,m+j)=other(i,j);
    }
  }

  matrix
  cut(int a, int b, int n, int m) const
  {
    matrix result(n, m);
    for (int i=0;i<n;++i) {
      for (int j=0;j<m;++j) 
	result(i,j)=(*this)(a+i,b+j);
    }
    return result;
  }

  int get_rows() const { return rows; }
  int get_columns() const { return columns; }

  std::pair<int,int>
  dim() const { return std::make_pair(rows, columns); }

  std::vector<T> row(int l) const {
    assert ( l >= 0 && l < rows );
    std::vector<T> r;
    for (int i=0; i < columns; ++i)
      r.push_back((*this)(l, i));
    return r;
  }

  void set_row(int l, const std::vector<T>& r) {
    assert( l >= 0 && l < rows );
    assert( r.size() == columns );
    for (int i=0; i < columns; ++i)
      (*this)(l, i) = r[i];
  }

  void set_column(int l, const std::vector<T>& c) {
    assert( l >= 0 && l < columns );
    assert( c.size() == rows );
    for (int i=0; i < rows; ++i)
      (*this)(i, l) = c[i];
  }

  std::vector<T> column(int l) const {
    assert( l >= 0 && l < columns );
    std::vector<T> c;
    for (int i=0; i < rows; ++i)
      c.push_back((*this)(i, l));
    return c;
  }

  void
  insert_column(int pos, const std::vector<T>& column)
  {
    assert(0 <= pos);
    assert(pos <= columns);
    assert(column.size() == rows);

    matrix result(rows, columns+1);
    result.inject(this->cut(0, 0, rows, pos), 0, 0);
    result.set_column(pos, column);
    result.inject(this->cut(0, pos, rows, columns-pos), 0, pos+1);
    this->swap(result);
  }

  void
  swap(matrix& other)
  {
    std::swap(rows, other.rows);
    std::swap(columns, other.columns);
    std::swap(p, other.p);
  }

  std::vector<T> to_vector() const {
    std::vector<T> c;
    for (int i=0; i < rows; ++i)
      for (int j=0; j < columns; ++j)
      c.push_back((*this)(i, j));
    return c;
    //return p;
  }


  T&
  operator()(int i, int j) {
    assert(i>=0 && j >=0);
    assert(i<rows);
    assert(j<columns);

    return p[columns*i+j];
  }

  T
  operator()(int i, int j) const {
    assert(i>=0 && j >=0);
    assert(i<rows);
    assert(j<columns);

    return p[columns*i+j];
  }

  void
  fill_with(T c)
  {
    /*
    for (int i=0;i<rows;++i) {
      for (int j=0;j<columns;++j) 
	(*this)(i,j)=c;
    }
    */
    const int sz = p.size();
    for (int i=0; i < sz; ++i)
      p[i] = c;
  }

  template <typename Functor>
  Functor
  apply(Functor f)
  {
    for (int i=0;i<rows;++i) {
      for (int j=0;j<columns;++j) 
	(*this)(i,j) = f((*this)(i,j));
    }
    return f;
  }

  template <typename Functor>
  Functor
  iterate(Functor f) const
  {
    for (int i=0;i<rows;++i) {
      for (int j=0;j<columns;++j) 
	f((*this)(i,j));
    }
    return f;
  }

  void
  print(int field_width = 2, int precision=6) const
  {
    for (int i=0;i<rows;++i) {
      for (int j=0;j<columns;++j) 
	std::cout << std::setw(field_width) << std::setprecision(precision) << std::fixed << (*this)(i,j) << "\t";
      //std::cout << std::setw(field_width) << std::setprecision(precision) << std::scientific << (*this)(i,j) << "\t";
      std::cout << std::endl;
	/*
	printf("%*i\t", field_width, (*this)(i,j));
      printf("\n");
	*/
    }
       
  }

  matrix<T>&
  operator+=(const matrix<T>& m2);

  
  // void
  // print2(FILE* fp, std::string format) const
  // {
  //   // format == "%*.*lf\t" , field_width, precision
  //   for (int i=0;i<rows;++i) {
  //     for (int j=0;j<columns;++j) 
  // 	fprintf(fp, format.c_str(), (*this)(i,j));
  //     fprintf(fp, "\n");
  //   }
       
  // }

  void
  print3(FILE* fp, std::string format, const char* sep = "\t") const
  {
    for (int i=0;i<rows;++i) {
      fprintf(fp, format.c_str(), (*this)(i,0));
      for (int j=1;j<columns;++j) {
	fprintf(fp, "%s", sep);
	fprintf(fp, format.c_str(), (*this)(i,j));
      }
      fprintf(fp, "\n");
    }
  }

  // void
  // print_line(int row, int field_width = 2) const
  // {
  //   for (int j=0;j<columns;++j) 
  //     printf("%*i", field_width, (*this)(row,j));
  //   printf("\n");
       
  // }


  ~matrix()
  {
  }

private:
  int rows;
  int columns;
  std::vector<T> p;
};

typedef matrix<double> dmatrix;
typedef matrix<int> imatrix;


// HUOM ei toimi oikein
template <typename T>
class sub_matrix
  : public matrix<T>
{
public:

  sub_matrix(matrix<T>& base_, int i, int j, int m, int n) : 
    matrix<T>(base_), startx(i), starty(j), 
    rows(n), columns(m) 
  {

  }

  T&
  operator()(int i, int j) {
    assert(i>=0 && j >=0);
    assert(i<rows);
    assert(j<columns);

    return matrix<T>(startx+i, starty+j);
  }

  T
  operator()(int i, int j) const {
    assert(i>=0 && j >=0);
    assert(i<rows);
    assert(j<columns);

    return matrix<T>(startx+i, starty+j);
  }

private:
  //matrix<T>& base;
  int startx;
  int starty;
  int rows;
  int columns;

};

template <typename T>
bool
operator==(const matrix<T>& m1, const matrix<T>& m2)
{
  int rows1 = m1.get_rows();
  int cols1 = m1.get_columns();
  int rows2 = m2.get_rows();
  int cols2 = m2.get_columns();
  if (rows1 != rows2)
    return false;
  if (cols1 != cols2)
    return false;

  for (int i=0; i < rows1; ++i)
    for (int j=0; j < cols1; ++j)
      if (m1(i,j) != m2(i,j))
	return false;

  return true;
}

template <typename T>
matrix<T>
operator*(double d, const matrix<T>& m)
{
  int rows = m.get_rows();
  int cols = m.get_columns();
  matrix<T> ret(rows, cols);

  for (int i=0; i < rows; ++i)
    for (int j=0; j < cols; ++j)
      ret(i,j) = d * m(i,j);

  return ret;
}
template <typename T>
matrix<T>
operator+(double d, const matrix<T>& m)
{
  int rows = m.get_rows();
  int cols = m.get_columns();
  matrix<T> ret(rows, cols);

  for (int i=0; i < rows; ++i)
    for (int j=0; j < cols; ++j)
      ret(i,j) = d + m(i,j);

  return ret;
}


// nämä olisi parempi toteuttaa geneerisella binaarisella operaatiolla
// jolle annetaan parametriksi esim. + tai /
template <typename T>
matrix<T>
operator-(const matrix<T>& m1, const matrix<T>& m2)
{
  int cols1 = m1.get_columns();
  int rows1 = m1.get_rows();
  int cols2 = m2.get_columns();
  int rows2 = m2.get_rows();
  assert(rows1==rows2 && cols1==cols2);

  matrix<T> ret(rows1, cols1);

  for (int i=0; i < rows1; ++i)
    for (int j=0; j < cols1; ++j)
      ret(i,j) = m1(i,j) - m2(i,j);

  return ret;
}

template <typename T>
matrix<T>&
matrix<T>::operator+=(const matrix<T>& m2)
{
  int cols1 = this->get_columns();
  int rows1 = this->get_rows();

  int cols2 = m2.get_columns();
  int rows2 = m2.get_rows();

  assert(rows1==rows2 && cols1==cols2);

  for (int i=0; i < rows1; ++i)
    for (int j=0; j < cols1; ++j)
      (*this)(i,j) += m2(i,j);

  return *this;
}

template <typename T>
matrix<T>
operator+(const matrix<T>& m1, const matrix<T>& m2)
{
  matrix<T> result = m1;
  result += m2;
  return result;
}

template <typename T>
matrix<T>
operator/(const matrix<T>& m1, const matrix<T>& m2)
{
  int cols1 = m1.get_columns();
  int rows1 = m1.get_rows();
  int cols2 = m2.get_columns();
  int rows2 = m2.get_rows();
  assert(rows1==rows2 && cols1==cols2);

  matrix<T> ret(rows1, cols1);

  for (int i=0; i < rows1; ++i)
    for (int j=0; j < cols1; ++j)
      ret(i,j) = (m1(i,j) == m2(i,j)) ? 1 : m1(i,j) / m2(i,j);   // 0/0 -> 1
 
  return ret;
}

template <typename T>
matrix<T>
transpose(const matrix<T>& m)
{
  int new_rows = m.get_columns();
  int new_cols = m.get_rows();
  matrix<T> ret(new_rows, new_cols);

  for (int i=0; i < new_rows; ++i)
    for (int j=0; j < new_cols; ++j)
      ret(i,j)=m(j,i);

  return ret;
}


// computes the maximum element-wise distance between two matrices
template <typename T>
T
distance(const matrix<T>& s, const matrix<T>& t)
{
  T max_dist=0;
  int rows = s.get_rows();
  int columns = s.get_columns();
  assert(s.get_rows() == t.get_rows());
  assert(s.get_columns() == t.get_columns());
 
 for (int i=0;i<rows;++i) {
    for (int j=0;j<columns;++j) 
      max_dist = std::max( fabs(s(i,j)-t(i,j)) , max_dist);
  }

 return max_dist;
}


template <typename T>
T
relative_distance(const matrix<T>& s, const matrix<T>& t)
{
  T max_dist=0;
  int rows = s.get_rows();
  int columns = s.get_columns();
  assert(s.get_rows() == t.get_rows());
  assert(s.get_columns() == t.get_columns());
 
 for (int i=0;i<rows;++i) {
    for (int j=0;j<columns;++j) 
      max_dist = std::max( fabs(s(i,j)-t(i,j))/s(i,j) , max_dist);
  }

 return max_dist;
}

#endif // MATRIX_HPP
