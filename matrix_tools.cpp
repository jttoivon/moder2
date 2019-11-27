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
#include "matrix_tools.hpp"
#include "common.hpp"

std::vector<double>
parse_row(char* tmp)
{
  std::vector<double> v;
  int m;
  double d;
  while (true) {
    if (sscanf(tmp, "%lg%n", &d, &m) != 1)
      break;
    v.push_back(d);
    tmp += m;
  }
  return v;
}

std::vector<double>
read_row(FILE* fp)
{
  char* line=NULL;
  size_t n;
  if (getline(&line, &n, fp) == -1)
    return std::vector<double>();
  std::vector<double> v = parse_row(line);
  free(line);

  return v;
}

dmatrix
dimensionless_read(FILE* fp, char* line)
{
  std::vector<double> v = parse_row(line);
  free(line);

  int cols = v.size();
  if (cols == 0)
    throw fileformat_error("Couldn't read matrix elements");

  std::vector<std::vector<double> > rows;
  rows.push_back(v);

  while (true) {
    v = read_row(fp);
    if (v.size() == cols)
      rows.push_back(v);
    else 
      break;
  }
  
  matrix<double> m(rows.size(), cols);
  for (int i=0; i < rows.size(); ++i)
    m.set_row(i, rows[i]);
  return m;

}

matrix<double>
read_matrix(FILE* fp)
{
  assert(fp != NULL);
  int rows, cols;
  double d;
  char* line=NULL;
  size_t n;
  if (getline(&line, &n, fp) == -1)
    throw fileformat_error("Couldn't read matrix elements");

  int count=sscanf(line, "%ix%i\n", &rows, &cols);
  if (count != 2)
    return dimensionless_read(fp, line);
  free(line);

  assert(count==2);
  //printf("%ix%i\n", rows, cols);

  matrix<double> m(rows, cols);
  for (int i=0; i<rows; ++i)
    for (int j=0; j<cols; ++j) {
      if (fscanf(fp, "%lg", &d) != 1)
	throw fileformat_error("Couldn't read matrix elements");
      m(i,j)=d;
    }
  return m;
}




void
write_matrix_file(const std::string& matrixfile, const dmatrix& M, std::string format)
{
  FILE* fp=fopen(matrixfile.c_str(), "w");
  if (fp == NULL) {
    fprintf(stderr, "Couldn't create matrixfile %s\n", matrixfile.c_str());
    perror("write_matrix_file");
    exit(1);
  }
  write_matrix(fp, M, "", format, false); // write used to file
  fclose(fp);
  
  // This should be conditional to some log level request
  //printf("Wrote resulting matrix to file %s\n", matrixfile.c_str());
}

dmatrix
read_matrix_file(const std::string& matrixfile)
{
  FILE* fp = fopen(matrixfile.c_str(), "r");
  if (fp == NULL) {
    fprintf(stderr, "Couldn't open file %s\n", matrixfile.c_str());
    exit(1);
  }
  matrix<double> M = read_matrix(fp);
  fclose(fp);
  return M;
}

    




// void
// normalize_matrix(matrix<double>& m)
// {
//   // sums of rows should be 1
//   for (int i=0; i < m.get_rows(); ++i) {
//     double sum = 0;
//     for (int j=0; j < m.get_columns(); ++j)
//       sum += m(i,j);
    
//     for (int j=0; j < m.get_columns(); ++j)
//       m(i,j) /= sum;
//   }

// }

// void
// normalize_matrix2(matrix<double>& m)
// {
//   // sums of cols should be 1
//   for (int i=0; i < m.get_columns(); ++i) {
//     double sum = 0;
//     for (int j=0; j < m.get_rows(); ++j)
//       sum += m(j,i);
    
//     for (int j=0; j < m.get_rows(); ++j)
//       m(j,i) /= sum;
//   }

// }

void
normalize_whole_matrix(matrix<double>& m)
{
  // sums of cols should be 1
  double sum = 0;
  for (int i=0; i < m.get_columns(); ++i) {
    for (int j=0; j < m.get_rows(); ++j)
      sum += m(j,i);
  }
    
  for (int i=0; i < m.get_columns(); ++i) {
    for (int j=0; j < m.get_rows(); ++j)
      m(j,i) /= sum;
  }

}

void
normalize_matrix_rows(matrix<double>& m)
{
  // sums of rows should be 1
  for (int i=0; i < m.get_rows(); ++i) {
    double sum = 0;
    for (int j=0; j < m.get_columns(); ++j)
      sum += m(i,j);

    if (sum != 0.0) {
      for (int j=0; j < m.get_columns(); ++j)
	m(i,j) /= sum;
    }
  }
}


dmatrix
normalize_matrix_columns_copy(matrix<double> m)
{
  normalize_matrix_columns(m);
  return m;
}

dmatrix
normalize_matrix_rows_copy(matrix<double> m)
{
  normalize_matrix_rows(m);
  return m;
}


bool
is_palindromic_matrix(const matrix<double>& m)
{
  int rows=m.get_rows();
  int cols=m.get_columns();
  assert(rows==4);

  int ret=1;

  int middle=static_cast<int>(ceil(cols/2.0));
  for (int j=0; j < middle; ++j) {
    ret *= (m(0, j) == m(3, cols-j-1));
    ret *= (m(1, j) == m(2, cols-j-1));
  }

  return ret;
}

dmatrix
matrix_sum(const dmatrix& m1, const dmatrix& m2, int d)
{
  int width = m1.get_columns() + d + m2.get_columns();

  dmatrix result(4, width);

  result.inject(m1, 0, 0);

  for (int i=0; i < m2.get_columns(); ++i)
    for (int c=0; c < 4; ++c)
      result(c, m1.get_columns()+d+i) += m2(c, i);

  return result;
}

dmatrix
matrix_product(const dmatrix& m1, const dmatrix& m2, int d)
{
  int width = m1.get_columns() + d + m2.get_columns();

  dmatrix result(4, width);
  result.fill_with(1.0);
  result.inject(m1, 0, 0);

  for (int i=0; i < m2.get_columns(); ++i)
    for (int c=0; c < 4; ++c)
      result(c, m1.get_columns()+d+i) *= m2(c, i);

  return result;
}
