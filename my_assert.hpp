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
