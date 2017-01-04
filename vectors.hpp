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
