#include<stdio.h>
#include<stdlib.h>

#undef NDEBUG
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
