/*

    MODER is a program to learn DNA binding motifs from SELEX datasets.
    Copyright (C) 2016, 2017  Jarkko Toivonen

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
#include <vector>
#include <string>

bool
is_smaller_freq(const std::vector<int>& freq1, const std::vector<int>& freq2);

bool
is_non_negative(const std::vector<int>& freq);


bool
equal_vectors(const std::vector<int>& v1, const std::vector<int>& v2);


std::vector<int>
subtract_frequencies(const std::vector<int>& freq1, const std::vector<int>& freq2);

std::vector<int>
add_frequencies(const std::vector<int>& freq1, const std::vector<int>& freq2);

template<typename T>
std::vector<T>
add_vectors(const std::vector<T>& v1, const std::vector<T>& v2)
{
  int len = v1.size();
  assert(len == v2.size());

  std::vector<T> result(len);
  for (int i=0; i < len; ++i)
    result[i] = v1[i] + v2[i];

  return result;

}

template<typename T>
std::vector<T>&
operator+=(std::vector<T>& v1, const std::vector<T>& v2)
{
  int len = v1.size();
  assert(len == v2.size());

  for (int i=0; i < len; ++i)
    v1[i] += v2[i];

  return v1;

}

template <typename T>
std::vector<T>&
operator/=(std::vector<T>& v, double d)
{
  for (int i=0; i < v.size(); ++i)
    v[i] /= d;
  return v;
}


// compute frequency table for vectors containing values 0, 1, 2, 3
std::vector<int>
get_freqs(const std::vector<int>& s);


std::vector<int>
get_freqs(const std::string& s);
