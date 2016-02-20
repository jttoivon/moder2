#include "orientation.hpp"
#include "matrix.hpp"
#include "common.hpp"

#include <boost/tuple/tuple.hpp>
#include <cassert>

boost::tuple<dmatrix,dmatrix>
get_matrices_according_to_hetero_orientation(int o, const dmatrix& m1, const dmatrix& m2)
{
  assert(o >= 0 and o <=3);

  dmatrix x, y;
  switch (o) {
  case HT: x = m1; y = m2; break;                      // HT
  case HH: x = m1; y = reverse_complement(m2); break;  // HH
  case TT: x = reverse_complement(m1); y = m2; break;  // TT
  case TH: x = reverse_complement(m1); y = reverse_complement(m2); break;  // TH
  }
  return boost::make_tuple(x, y);
}

boost::tuple<std::string,std::string>
get_seeds_according_to_hetero_orientation(int o, const std::string& s1, const std::string& s2)
{
  assert(o >= 0 and o <=3);

  std::string x, y;
  switch (o) {
  case HT: x = s1; y = s2; break;                      // HT
  case HH: x = s1; y = reverse_complement(s2); break;  // HH
  case TT: x = reverse_complement(s1); y = s2; break;  // TT
  case TH: x = reverse_complement(s1); y = reverse_complement(s2); break;  // TH
  }
  return boost::make_tuple(x, y);
}

int
orientation(int o1, int o2)
{
  assert(o1 == 0 || o1 == 1);
  assert(o2 == 0 || o2 == 1);

  if (o1 == o2)
    return HT;
  if (o1 == 0 && o2 == 1)
    return HH;
  if (o1 == 1 && o2 == 0)
    return TT;

  assert(false);
}

int
orientation2(int o1, int o2)
{
  assert(o1 == 0 || o1 == 1);
  assert(o2 == 0 || o2 == 1);

  if (o1 == 0 && o2 == 0)
    return HT;
  if (o1 == 0 && o2 == 1)
    return HH;
  if (o1 == 1 && o2 == 0)
    return TT;
  if (o1 == 1 && o2 == 1)
    return TH;

  assert(false);
}

const char* orients[4] = {"HT", "HH", "TT", "TH"};

const int hetero_orientation_class::a[4][2] = { {128+2,16+4}, {128+1,32+4}, {64+2,16+8}, {64+1,32+8}};

const int homo_orientation_class::a[4] = { 8+2, 8+1, 4+2, 4+1};

hetero_orientation_class hetero_orientation;
homo_orientation_class get_homo_orientation;


