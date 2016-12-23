/*

    MODER is a program to learn DNA binding motifs from SELEX datasets.
    Copyright (C) 2016  Jarkko Toivonen

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
#ifndef MY_ASSERT_HPP
#define MY_ASSERT_HPP

#ifndef NDEBUG
#include<stdio.h>
#include<stdlib.h>
#include<assert.h>

#include <string>
#include <sstream>




namespace my_assert_detail
{

void
__my_assert(const char *assertion, const char *x, const char *y, const char *file, 
	    unsigned int line, const char *function);

void
__my_assert2(const char *assertion, const char *x, const char *file, 
	    unsigned int line, const char *function);

	    

template <typename T>
std::string
my_to_string(T t)
{
  std::ostringstream os;
  os << t;
  return os.str();
}

} // end namespace my_assert_detail

#define my_assert(x,y) \
  ((x)==(y))? (void)0 :	   \
  my_assert_detail::__my_assert(__STRING((x)==(y)),	\
	      my_assert_detail::my_to_string(x).c_str(),		\
	      my_assert_detail::my_to_string(y).c_str(),		\
	      __FILE__,				\
	      __LINE__,				\
	      __ASSERT_FUNCTION)

#define my_assert2(x,y,op)   \
  ((x)op(y))? (void)0 :	   \
  my_assert_detail::__my_assert(__STRING((x)op(y)), \
				my_assert_detail::my_to_string(x).c_str(), \
				my_assert_detail::my_to_string(y).c_str(), \
				__FILE__, __LINE__, __ASSERT_FUNCTION)

#define my_assert3(x,cond)   \
  (cond)? (void)0 :	   \
  my_assert_detail::__my_assert2(__STRING(cond), \
				 my_assert_detail::my_to_string(x).c_str(), \
				 __FILE__, __LINE__, __ASSERT_FUNCTION)

#else
#define my_assert(x,y)
#define my_assert2(x,y,z)    
#define my_assert3(x,cond)    
#endif // NDEBUG

#endif // MY_ASSERT_HPP
