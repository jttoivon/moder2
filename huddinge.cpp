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
#include "huddinge.hpp"
#include "data.hpp"
#include "common.hpp"
#include "matrix.hpp"
#include "suffix_array_wrapper.hpp"

#include <functional>

bool require_contiguous_gap=false;
bool use_counts=false;
bool quiet = true;
bool print_maximum_only = false;
bool get_paths = false;


// These are for bit operations
enum {up   = 1,
      left = 2,
      diag = 4};


std::string
mypad(const std::string& s, int pad, char c, int total_length)
{
  int right = total_length - pad - s.length();
  return std::string(pad, c) + s + std::string(right, c);
}


int
min3(int a, int b, int c)
{
  return std::min(a, std::min(b, c));
}

int
max3(int a, int b, int c)
{
  return std::max(a, std::max(b, c));
}

boost::tuple<dmatrix, imatrix>
huddinge_distance_helper(const std::string& x, const std::string& y)
{
  int xlen = x.length();
  int ylen = y.length();

  dmatrix d(ylen+1, xlen+1);
  imatrix path(ylen+1, xlen+1);
  
  for (int i=1; i < ylen; ++i) {
    for (int j=1; j < xlen; ++j) {
      d(i, j) = d(i-1, j-1) + (y[i-1] == x[j-1] and y[i-1] != 'n');
    }
  }

  // Lowest row
  for (int j=1; j < xlen; ++j) {
    d(ylen, j) = std::max(d(ylen, j-1), d(ylen-1, j-1) + (y[ylen-1] == x[j-1] and y[ylen-1] != 'n') );
    if (d(ylen, j) == d(ylen, j-1))
      path(ylen, j) |= left;
    if (d(ylen, j) == d(ylen-1, j-1) + (y[ylen-1] == x[j-1] and y[ylen-1] != 'n'))
      path(ylen, j) |= diag;
  }
  
  // rightmost column
  for (int i=1; i < ylen; ++i) {
    d(i, xlen) = std::max(d(i-1, xlen), d(i-1, xlen-1) + (y[i-1] == x[xlen-1] and x[xlen-1] != 'n') );
    if (d(i, xlen) == d(i-1, xlen))
      path(i, xlen) |= up;
    if (d(i, xlen) == d(i-1, xlen-1) + (y[i-1] == x[xlen-1] and x[xlen-1] != 'n'))
      path(i, xlen) |= diag;
  }
  
  d(ylen, xlen) = max3(d(ylen-1, xlen), d(ylen, xlen-1), d(ylen-1, xlen-1) + (y[ylen-1] == x[xlen-1] and x[xlen-1] != 'n'));

  // Which routes gave maximum value
  if (d(ylen, xlen) == d(ylen-1, xlen))
    path(ylen, xlen) |= up;
  if (d(ylen, xlen) == d(ylen, xlen-1))
    path(ylen, xlen) |= left;
  if (d(ylen, xlen) == d(ylen-1, xlen-1) + (y[ylen-1] == x[xlen-1] and x[xlen-1] != 'n'))
    path(ylen, xlen) |= diag;

  return boost::make_tuple(d, path);

}

character_to_values<bool> is_defined_base("ACGTacgt", true);


int
defined_bases(const std::string& x)
{
  int count = 0;
  int len = x.length();
  for (int i=0; i < len; ++i)
    count += is_defined_base((unsigned char)x[i]);
  return count;
}

bool
one_contiguous_gap(std::string s)
{
  assert(s[0] != 'n' and s[s.length()-1] != 'n');

  int number_of_switches=0;
  bool was_previous_defined = true;
  for (int i=0; i < s.length(); ++i) {
    char c = s[i];
    if (was_previous_defined != is_defined_base(c))
      number_of_switches += 1;
    was_previous_defined = is_defined_base(c);
  }
  return number_of_switches <= 2;

}


int
huddinge_distance(const std::string& x, const std::string& y)
{
  int xlen = x.length();
  int ylen = y.length();

  dmatrix d = huddinge_distance_helper(x, y).get<0>();
  return std::max(defined_bases(x), defined_bases(y)) - d(ylen, xlen);
}

void
check(int i, int j, std::vector<int>& result, const imatrix& path)
{
  if (path(i, j) & diag) {
    int r = std::min(i, j);
    if (i >= j)   // Last row, y starts before x
      result.push_back(-(i-r));
    else       // Last column, x starts before y
      result.push_back((j-r));
  }
}

std::vector<int>
huddinge_alignment(const std::string& x, const std::string& y)
{
  int xlen = x.length();
  int ylen = y.length();
  dmatrix d;
  imatrix path;
  boost::tie(d, path) = huddinge_distance_helper(x, y);
  std::vector<int> result;  // In alignment a \in result, y starts after 'a' positions after x
  check(ylen, xlen, result, path);  // Check the lower right corner

  // Check the last column
  int i = ylen;
  while (path(i, xlen) & up and i >= 0) {
    i -= 1;
    check(i, xlen, result, path);
  }
		// Check the last row
  int j = xlen;
  while (path(ylen, j) & left and j >= 0) {
    j -= 1;
    check(ylen, j, result, path);
  }

  return result;

}


// Align consequent strings in the input list
std::vector<std::string>
huddinge_align_all(const std::vector<std::string>& l)
{
  int n=l.size();
  std::vector<int> aligns(1, 0);
  std::vector<int> lengths;

  BOOST_FOREACH(const std::string& x, l) {
    lengths.push_back(x.length());
  }
  int current_align=0;
  for (int i=0; i < n-1; ++i) {
    std::vector<int> a = huddinge_alignment(l[i], l[i+1]);
        
    current_align=current_align+a[0];
    aligns.push_back(current_align);
  }
  int max_align=-min_element(aligns);
  int total_length = 0;
  for (int i=0; i < n; ++i) {
    total_length = std::max(total_length, max_align + aligns[i] + lengths[i]);
  }
  std::vector<std::string> result;
  for (int i=0; i < n; ++i) {
    result.push_back(mypad(l[i], max_align+aligns[i], 'n', total_length));
  }
  return result;
}


// The next algorithm is from page
// https://alistairisrael.wordpress.com/2009/09/22/simple-efficient-pnk-algorithm/
void
next_k_permutation(std::vector<int>& a, int k)
{
  int n = a.size();
  assert(n >= k);
  int edge = k - 1;
  // find j in (k to n-1) where a_j > a_edge
  int j = k;
  while (j < n and a[edge] >= a[j])
    j += 1;

  if (j < n)
    std::swap(a[edge], a[j]);
  //a[edge], a[j] = a[j], a[edge];
  else {
    // reverse a[k] to a[n-1]
    std::reverse(a.begin()+k , a.end());
    
    // find rightmost ascent to left of edge
    int i = edge - 1;
    while (i >= 0 and a[i] >= a[i+1])
      i -= 1;

    if (i < 0) //no more permutations
      return;

    // find j in (n-1 to i+1) where a[j] > a[i]
    j = n - 1;
    while (j > i and a[i] >= a[j])
      j -= 1;
    
    std::swap(a[i], a[j]);
    //a[i], a[j] = a[j], a[i]  # swap
    //reverse a[i+1] to a[n-1]
    std::reverse(a.begin()+(i+1), a.end());
    //    a[i+1:n] = reversed(a[i+1:n])
  }
  //    return a
}


class operation
{
public:

  virtual
  void
  execute(std::string& s, std::string& t) const=0;

  virtual
  std::ostream&
  print(std::ostream& f) const=0;

};

std::ostream&
operator<<(std::ostream& f, const operation& o)
{
  return o.print(f);
}

class swap : public operation
{
public:

  swap(int i_, int j_) : i(i_), j(j_) {}
  
  void
  execute(std::string& s, std::string& t) const
  {
    s[i] = t[i];
    s[j] = 'n';
  }

  std::ostream&
  print(std::ostream& f) const
  {
    return f << to_string("Swap(%i,%i)", i, j);
  }
  
private:
  int i;
  int j;
};

class doswitch : public operation
{
public:

  doswitch(int i_) : i(i_) {}
  
  void
  execute(std::string& s, std::string& t) const
  {
    s[i] = t[i];
  }

  std::ostream&
  print(std::ostream& f) const
  {
    return f << to_string("Switch(%i)", i);
  }
  
private:
  int i;
};

class change : public operation
{
public:

  change(int i_) : i(i_) {}
  
  void
  execute(std::string& s, std::string& t) const
  {
    s[i] = t[i];
  }

  std::ostream&
  print(std::ostream& f) const
  {
    return f << to_string("Change(%i)", i);
  }
  
private:
  int i;
};

suffix_array sa;

int
get_count(const std::string& x)
{
  int count = sa.count_iupac(x);
  if (is_palindromic(x))
    return count / 2;
  else 
    return count;
}

double
get_scaled_count(const std::string& x, int reference_length)
{
  return get_count(x) * pow(4.0, defined_bases(x) - reference_length);
}

typedef std::pair<std::string, double> my_comp_param;

bool
mycomp(const my_comp_param& a, const my_comp_param& b)
{
  return a.second < b.second;
}

std::string
make_header(int l)
{
  std::string result(l, '-');
  for (int i=0; i < l; ++i)
    result[i] = i % 10;
  return result;
}

counted_path_type
get_counted_path(const path_type& path, int reference_length)
{
  counted_path_type result;
  BOOST_FOREACH(std::string w, path) {
    double count = use_counts ? get_scaled_count(w, reference_length) : 0;
    result.push_back(boost::make_tuple(w, count));
  }
  return result;
}


// This is more or less Dijkstra's algorithm except vertices are added to the priority queue lazily
std::vector<counted_path_type>
check_for_longer_paths(std::string x, std::string y, int dmax, int wmin, int L, int reference_length)
{
  if (not quiet)
    printf("Finding route with highest weight, but at least weight %i, and length at most %d\n", wmin, dmax);
  typedef std::map<std::string, int> int_map_type;
  typedef std::map<std::string, double> double_map_type;
  typedef std::map<std::string, std::vector<std::string> > list_map_type;
  int_map_type d;
  double_map_type w;
  list_map_type prev;
  prev[x]=std::vector<std::string>();
  d[x]=0;
  w[x]=get_scaled_count(x, reference_length);  // non-final weights
  double_map_type S;        // final weights

  while (w.size() > 0) {
    std::string u;
    double w1;
    boost::tie(u, w1) = *std::max_element(w.begin(), w.end(), mycomp);
    //u, w1 = sorted(w.items(), key = lambda x: x[1])[-1];  // get the maximum
    S[u] = w1;
    if (not quiet) {
      printf("Current string is %s with weight %f", u.c_str(), w1);
      printf(" and count  %f and path length %i\n", get_scaled_count(u, reference_length), d[u]);
    }
    w.erase(u);
    if (d[u] == dmax)          // Path is already too long to be extended
      continue;
    huddinge_neighbourhood neighbourhood(u, 1, L, 1, L);
    std::vector<boost::tuple<std::string, int> > neigh=neighbourhood.compute(require_contiguous_gap); // The neighbourhood contains u as well
    std::string v;
    BOOST_FOREACH(boost::tie(v, boost::tuples::ignore), neigh) {
      if (S.count(v))
	continue;
      double c = get_scaled_count(v, reference_length);
      if (c < wmin)
	continue;
      if (not w.count(v)) {
	w[v] = 0;
	d[v] = 1000000;
	prev[v] = std::vector<std::string>();
      }
      double temp = std::min(S[u], c);
      if (w[v] < temp or (w[v] == temp and d[v] >= d[u] + 1)) {   // update weight of a node
	if (w[v] < temp or d[v] > d[u] + 1)
	  prev[v] = std::vector<std::string>(1, u);
	else
	  prev[v].push_back(u);
	w[v] = temp;
	d[v] = d[u] + 1;
      }
    }
  }

  std::vector<counted_path_type> result;
  if (not quiet)
    printf("Weight of y is %f\n", S[y]);
  if (S.count(y)) {
    std::string temp=y;
    int depth=d[y];
    std::vector<std::string> path(depth+1, "");
    std::vector<std::string> nodes(1, y);
    while (nodes.size() > 0) {
      temp=nodes.back();
      nodes.pop_back();
      path[d[temp]]=temp;
      BOOST_FOREACH(std::string child, prev[temp])
	nodes.push_back(child);
      if (temp == x)
	result.push_back(get_counted_path(path, reference_length));
    }
  }

  if (not quiet)
    printf("Size of result set from find_longer_paths is %zu\n", result.size());

  return result;
}

std::string
strip(const std::string& s, char c)
{
  int start, end;
  for (start=0; start < s.length() and s[start] == c; ++start)
    ;
  for (end=s.length(); end > 0 and s[end-1] == c; --end)
    ;
  if (start >= end)
    return std::string();
  else
    return std::string(start, end - start);
}

bool
is_valid_path(const std::vector<std::string>& path)
{
  if (not require_contiguous_gap)
    return true;
  
  BOOST_FOREACH(std::string s, path) {
    if (not one_contiguous_gap(strip(s, 'n'))) {
      return false;
    }
  }
  return true;
}

bool
has_positive_count(const counted_path_type& counted_path, int reference_length)
{
  std::string s;
  double count;
  BOOST_FOREACH(boost::tie(s, count), counted_path) {
    if (count == 0.0)
      return false;
  }

  return true;
}

// Finds paths of length Huddinge(x,y)
std::vector<counted_path_type>
find_paths(const std::string& x, const std::string& y)
{
  int xlen = x.length();
  int ylen = y.length();
  int xd = defined_bases(x);
  int yd = defined_bases(y);
  int reference_length=std::max(xd,yd);
  std::vector<int> alignments=huddinge_alignment(x, y);
  if (not quiet)
    printf ("Alignments are %s\n", print_vector(alignments).c_str());

  std::vector<counted_path_type> result;
  BOOST_FOREACH(int h, alignments) {          // Iterate over all alignments
    int l=0;
    if (not quiet)
      printf("%s\n", std::string(40, '-').c_str());
    std::string s;
    std::string t;
    if (h >= 0) {
      l = std::max(xlen, h+ylen);
      s = mypad(x, 0, 'n', l);
      t = mypad(y, h, 'n', l);
    }
    else{
      int l = std::max(xlen - h, ylen);
      s = mypad(x, -h, 'n', l);
      t = mypad(y, 0, 'n', l);
    }
    std::string header = make_header(l);
    if (not quiet) {
      printf("%s\n", header.c_str());
      printf("%s\n", s.c_str());
      printf("%s\n", t.c_str());
      printf("\n");
    }
    std::vector<int> Jzero;
    std::vector<int> Jplus;
    std::vector<int> Jminus;
    for (int i=0; i < l; ++i) {
      if (s[i] != t[i]) {
	if (s[i] == 'n')
	  Jplus.push_back(i);
	else if (t[i] == 'n')
	  Jminus.push_back(i);
	else
	  Jzero.push_back(i);
      }
    }
    if (not quiet) {
      printf("Jzero %s\n", print_vector(Jzero).c_str());
      printf("Jplus %s\n", print_vector(Jplus).c_str());
      printf("Jminus %s\n", print_vector(Jminus).c_str());
    }
    int k = std::min(Jplus.size(), Jminus.size());
    int n = std::max(Jplus.size(), Jminus.size());
    int number_of_k_permutations=1;   //reduce(lambda x, y: x*y, xrange(n, n-k, -1), 1);
    for (int i=n; i > n-k; --i)
      number_of_k_permutations *= i;
  
    if (not quiet) {
      printf("n is %i\n", n);
      printf("k is %i\n", k);
      printf("Number of k-permutations is %i\n", number_of_k_permutations);
    }
    for (int counter = 0; counter < number_of_k_permutations; ++counter) {  // go through all sets of operations
      std::vector<boost::shared_ptr<operation> > operations;
      for (int i=0; i < k; ++i)
	operations.push_back(boost::shared_ptr<operation>(new swap(Jplus[i], Jminus[i])));
      if (xd >= yd)    
	for (int i=k; i < Jminus.size(); ++i)   // Deletions from x
	  operations.push_back(boost::shared_ptr<operation>(new doswitch(Jminus[i])));
      else
	for (int i=k; i < Jplus.size(); ++i)    // Insertions to x
	  operations.push_back(boost::shared_ptr<operation>(new doswitch(Jplus[i])));
      BOOST_FOREACH(int i, Jzero)
	operations.push_back(boost::shared_ptr<operation>(new change(i)));

      if (not quiet)
	printf("Operations: %s\n", print_vector(operations).c_str());
      int no = operations.size();
      std::vector<int> order(no);
      for (int i=0; i < no; ++i)
	order[i] = i;
	
      int number_of_orders = 1; // reduce(lambda x, y: x*y, xrange(1, no+1), 1);
      for (int i=1; i < no + 1; ++i)
	number_of_orders *= i;
      for (int order_counter = 0; order_counter < number_of_orders; ++order_counter) {  // go through orderings of operations
	std::vector<boost::shared_ptr<operation> > sorted_operations;
	BOOST_FOREACH(int i, order)
	  sorted_operations.push_back(operations[i]);
	std::string a = s;
	std::string b = t;
	path_type path;
	
	path.push_back(a);   // first sequence of the path
	BOOST_FOREACH(boost::shared_ptr<operation> op, sorted_operations) {
	  op->execute(a,b);
	  path.push_back(a);
	}
	bool valid_path = is_valid_path(path);

	if (valid_path) { // all strings in the path have at most one gap, if require_contiguous_gap
	  counted_path_type counted_path = get_counted_path(path, reference_length);
	  if (not use_counts or has_positive_count(counted_path, reference_length)) {
	    result.push_back(counted_path);
	    if (not quiet) {
	      printf("\n");
	      printf("%s\n", print_vector(sorted_operations).c_str());
	      printf("%s\n", header.c_str());   // header of indices for easier reading
	      std::string s;
	      double count;
	      BOOST_FOREACH(boost::tie(s, count), counted_path)
		printf("%s\t%f\n", s.c_str(), count);
	      printf("\n");
	    }

	  }
	} // if is_valid_path
	next_k_permutation(order, no);
      }  // for order_counter

      if (xd >= yd)    
	next_k_permutation(Jminus, k);
      else
	next_k_permutation(Jplus, k);
    }
  } // FOREACH h in alignments

  if (not quiet)
    printf("%s\n", std::string(40, '=').c_str());

  return result;
}



std::vector<boost::tuple<std::string,int> >
huddinge_neighbourhood::compute(bool allow_only_one_gap)
{
  for (int j=1; j < x.length() + 1; ++j)
    U[0][j] = insert(0, j);

  // Compute array U

  for (int i=1; i < L+1; ++i) {
    for (int j=0; j < x.length()+1; ++j) {
      str_dist_container ins = insert(i, j);
      str_dist_container de = del(i, j);
      str_dist_container ma = match(i, j);

      std::set<std::string> keys;
      std::string s;
      BOOST_FOREACH(boost::tie(s, boost::tuples::ignore), ins)
	keys.insert(s);
      BOOST_FOREACH(boost::tie(s, boost::tuples::ignore), de)
	keys.insert(s);
      BOOST_FOREACH(boost::tie(s, boost::tuples::ignore), ma)
	keys.insert(s);

      BOOST_FOREACH(s, keys) {
	int xerr = min3(my_get(ins, s).get<0>(), my_get(de, s).get<0>(), my_get(ma, s).get<0>());
	int yerr = min3(my_get(ins, s).get<1>(), my_get(de, s).get<1>(), my_get(ma, s).get<1>());
	assert(xerr != INT_MAX and yerr != INT_MAX);
	U[i][j][s] = boost::make_tuple(xerr, yerr);
	
	/*	
	int xpad=0, ypad = 0;
	if (i > j)
	  xpad = i - j;
	else
	  ypad = j - i;
	printf("j=%i i=%i xerr=%i yerr=%i:\n", j, i, xerr, yerr);
	std::string pad(xpad, ' ');
	printf("%s%s\n", pad.c_str(), x.substr(0, j).c_str());
	pad = std::string(ypad, ' ');
	printf("%s%s\n", pad.c_str(), s.c_str());
	*/
      }

    }
  }

  // Compute array T

  boost::multi_array<str_dist_container, 2> T(boost::extents[L+1][x.length()+1]);
  for (int i = 1; i < L+1; ++i) {
    std::string s;
    err_t e;
    BOOST_FOREACH(boost::tie(s, e), U[i][0]) {
      if (s[i-1] != 'N')
	T[i][0][s] = e;
    }
    for (int j=1; j < x.length()+1; ++j) {
      std::set<std::string> keys;
      BOOST_FOREACH(boost::tie(s, boost::tuples::ignore), U[i][j])
	keys.insert(s);
      BOOST_FOREACH(boost::tie(s, boost::tuples::ignore), T[i][j-1])
	keys.insert(s);

      BOOST_FOREACH(s, keys) {
	if (s[i-1] != 'N') {
	  err_t e1 = my_get(U[i][j], s);
	  err_t e2 = my_get(T[i][j-1], s);
	  int xerr_temp = e2.get<0>();
	  if (xerr_temp != INT_MAX)   // This makes sure that the value doesn't overflow
	    xerr_temp += (x[j-1] != 'N');
	  err_t temp = boost::make_tuple(std::min(e1.get<0>(), xerr_temp),
					 std::min(e1.get<1>(), e2.get<1>()));
	  if (std::max(temp.get<0>(), temp.get<1>()) <= h)
	    T[i][j][s] = temp;
	}
      }
    }
  }

  std::vector<boost::tuple<std::string, int> > result;
  for (int i=1; i < L+1; ++i) {
    std::string s;
    err_t err;
    BOOST_FOREACH(boost::tie(s, err), T[i][x.length()]) {
      if (non_N_count(s) >= min_kmer_len and non_N_count(s) <= max_kmer_len and (not allow_only_one_gap or has_atmost_one_gap(s))) {
	int distance = std::max(err.get<0>(), err.get<1>());
	//printf("%%%s %i %i\n", s.c_str(), err.get<0>(), err.get<1>());
	//printf("%%%s %i\n", s.c_str(), std::max(err.get<0>(), err.get<1>()));
	result.push_back(boost::make_tuple(s, distance));
      }
    }
  }
  return result;
}


str_dist_container
huddinge_neighbourhood::insert(int i, int j)
{
  str_dist_container result;
  if (i==0 and j > 0) {
    std::string s;
    err_t err;
    //printf("Size of U[%i][%i] is %zu\n", i, j-1, U[i][j-1].size());
    BOOST_FOREACH(boost::tie(s, err), U[i][j-1]) {
      //	if (std::max(err.get<0>(), err.get<1>()) == INT_MAX)
      //	  continue;

      err_t temp = boost::make_tuple(err.get<0>() + (x[j-1] != 'N'), err.get<1>());
      if (std::max(temp.get<0>(), temp.get<1>()) <= h) {
	result[s] = temp;
	// printf("j=%i i=%i xerr=%i yerr=%i:\n", j, 0, temp.get<0>(), temp.get<1>());
	// printf("%s\n", x.substr(0, j).c_str());
	// printf("\n");
      }
    }
  }

  return result;
}

str_dist_container
huddinge_neighbourhood::del(int i, int j)
{
  str_dist_container result;
  std::string sigma = (i==1 or i==L) ? "ACGT" : "ACGTN";

  if (i > 0 and (j==0 or j == x.length())) {
    std::string s;
    err_t err;
    BOOST_FOREACH(boost::tie(s, err), U[i-1][j]) {
      //      if (std::max(err.get<0>(), err.get<1>()) == INT_MAX)
      //	continue;

      BOOST_FOREACH(char w, sigma) {
	err_t temp = boost::make_tuple(err.get<0>(), err.get<1>() + (w != 'N'));
	if (std::max(temp.get<0>(), temp.get<1>()) <= h)
	  result[s+w] = temp;
      }
    }
  }

  return result;
}

bool
cerror(char a, char b)
{
  return a != 'N' and a != b;
}

str_dist_container
huddinge_neighbourhood::match(int i, int j)
{
  str_dist_container result;
  std::string sigma = (i==1 or i==L) ? "ACGT" : "ACGTN";

  if (i > 0 and j > 0) {
    std::string s;
    err_t err;
    BOOST_FOREACH(boost::tie(s, err), U[i-1][j-1]) {
      BOOST_FOREACH(char w, sigma) {
	// if (std::max(err.get<0>(), err.get<1>()) == INT_MAX)
	//   continue;
	err_t temp = boost::make_tuple(err.get<0>() + cerror(x[j-1], w), 
				       err.get<1>() + cerror(w, x[j-1]));
	if (std::max(temp.get<0>(), temp.get<1>()) <= h)
	  result[s+w] = temp;
      }
    }
  }

  return result;
}


err_t
huddinge_neighbourhood::my_get(const str_dist_container& h, const std::string& key)
{
  typedef str_dist_container::const_iterator iterator;
  iterator i = h.find(key);
  if (i == h.end())
    return boost::make_tuple(INT_MAX, INT_MAX);
  else
    return i->second;
}


