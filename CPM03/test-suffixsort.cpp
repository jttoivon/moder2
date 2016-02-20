/*
 * Copyright 2003 Stefan Burkhardt, Juha K"arkk"ainen
 *
 * This file contains a test program for dc-sort.
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version. 
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *
 */

#include "suffixsort.hpp"

#include <iostream>
#include <vector>
#include <cstdio>
#include <cstdlib>


int main(int argc, char* argv[])
{
  // Check Arguments
  if (argc != 2) {
    std::cout << "ERROR: Wrong number of arguments!\n";
    std::cout << "Usage: " << argv[0] << " textfile" << std::endl;
    std::exit(EXIT_FAILURE);
  }
  
  std::FILE *file1;
  file1 = std::fopen(argv[1], "r");
  if (!file1) {
    std::cout << "ERROR: Unable to open text file " << argv[1] << std::endl;
    std::exit(EXIT_FAILURE);
  }
  
  // Read Data
  std::fseek(file1, 0, SEEK_END);
  int length = std::ftell(file1);
  std::rewind(file1);
  std::vector<unsigned char> s(length);
  length = std::fread((void *)&s[0], sizeof(char), length, file1);
  std::fclose(file1);
  if (length != static_cast<int>( s.size() )) {
    std::cout << "ERROR: wrong number of bytes read" << std::endl;
    exit(EXIT_FAILURE);
  }

#ifndef SCRIPT
  std::cout << "Read " << length << " characters from " << argv[1]<< "\n"; 
  std::cout << "  Step 1pre     Step 1  Step 2pre     Step 2"
	    << "     Step 3     Step 4      Total   Checking" 
#ifdef LCP
            << "     Avg. LCP     Max. LCP"
#endif
	    << std::endl;
#endif
  
  // Construct the suffix array
  typedef std::vector<unsigned char>::difference_type offset_type;
  std::vector<offset_type> sa(length);
  try {
    construct_suffix_array(s.begin(), s.end(), sa.begin(), sa.end());
  } catch (const std::exception& error) {
    std::cout << std::endl << error.what() << std::endl;
    exit(EXIT_FAILURE);
  }

#ifdef LCP
  // Compute inverse suffix array
  std::vector<offset_type> isa(length);
  for (std::vector<offset_type>::iterator i = sa.begin(); i!=sa.end(); ++i) {
    isa[*i] = i-sa.begin();
  }	
  // Compute lcp array H (from Kasai et al. CPM 2001)
  std::vector<offset_type> H(length, 0);
  offset_type h=0;
  offset_type a,k;
  for(a=0; a<length; a++) {
    if(isa[a] > 1) {
      k = sa[isa[a]-1];
      while(s[a+h] == s[k+h]) {
	h++;
      }
      H[isa[a]] = h;
      h -= (h>0);
    }
  }
  // Analyze lcp array
  int max_lcp = 0;
  int curr_lcp = 0;
  double lcp_sum = 0;
  for(a=0; a<length; a++) {
    curr_lcp = H[a];
    lcp_sum += curr_lcp;
    if(curr_lcp > max_lcp) {
      max_lcp = curr_lcp;
    }
  }
  std::fprintf(stdout, " %12.2f %12d ", lcp_sum/(length-1), max_lcp);
#endif

#ifndef SCRIPT
  std::cout << std::endl;
#endif
}
