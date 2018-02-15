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
#include<stdio.h>
#include<stdlib.h>

#ifndef NDEBUG
#include "my_assert.hpp"

extern const char *__progname;

namespace my_assert_detail
{


void
__my_assert(const char *assertion, const char *x, const char *y, const char *file, 
	    unsigned int line, const char *function)
{
  fprintf(stderr, "%s: %s:%u: %s: Assertion %s failed.\n", __progname, file, line, 
	 function,assertion);
  fprintf(stderr, "Evaluated operand values are: %s and %s\n", x, y);
  fflush(stderr);
  abort();
}
	    
void
__my_assert2(const char *assertion, const char *x, const char *file, 
	     unsigned int line, const char *function)
{
  fprintf(stderr, "%s: %s:%u: %s: Assertion %s failed.\n", __progname, file, line, 
	 function,assertion);
  fprintf(stderr, "Operand value was: %s\n", x);
  fflush(stderr);
  abort();
}
	    
} // end namespace my_assert_detail

#endif
