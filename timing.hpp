//
// 2003 Juha K"arkk"ainen
//
// Timer
//

#ifndef _TIMING_HPP
#define _TIMING_HPP

#include <sys/time.h>
#include <iostream>
#include <iomanip>

inline double get_cpu_time()
{
  static const double in_seconds = 1.0/static_cast<double>(CLOCKS_PER_SEC);
  return clock() * in_seconds;
}

#ifdef TIMING

#define TIME_START(t) \
  double timer_##t##_current_, \
    timer_##t##_last_ = get_cpu_time()

#define TIME_CHECK(t) \
  timer_##t##_current_ = get_cpu_time(); \
  std::cout << std::setw(10) << std::fixed << std::setprecision(2) \
            << timer_##t##_current_ - timer_##t##_last_ << std::endl; \
  timer_##t##_last_ = timer_##t##_current_

#define TIME_PRINT(format,t)		 \
  timer_##t##_current_ = get_cpu_time(); \
  printf((format), timer_##t##_current_ - timer_##t##_last_);	\
  timer_##t##_last_ = timer_##t##_current_

#define TIME_GET(t)		 \
  ((timer_##t##_current_ = get_cpu_time()), (timer_##t##_current_ - timer_##t##_last_))


#else

#define TIME_START(t)
#define TIME_CHECK(t)
#define TIME_PRINT(format,t)
#define TIME_GET(t)

#endif

#endif // _TIMING_HPP
