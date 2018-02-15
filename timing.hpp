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
//
// 2003 Juha K"arkk"ainen
//
// Timer
//

#ifndef _TIMING_HPP
#define _TIMING_HPP

#include <time.h>
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
  printf(format, (timer_##t##_current_ - timer_##t##_last_));	\
  timer_##t##_last_ = timer_##t##_current_

#define TIME_GET(t)		 \
  ((timer_##t##_current_ = get_cpu_time()), (timer_##t##_current_ - timer_##t##_last_))


#ifdef __APPLE__ // apple machines don't define clock_gettime. Not at least the older ones.

inline
timeval
get_wall_time()
{
  struct timeval tv;
  gettimeofday(&tv, NULL);
  return tv;
}

inline
double
subtract_time(timeval& tp_new, timeval& tp_old)
{
  timeval temp = tp_new;
  temp.tv_sec -= tp_old.tv_sec;
  temp.tv_usec -= tp_old.tv_usec;
  if (temp.tv_usec < 0) {
    --temp.tv_sec;
    temp.tv_usec += 1000000;
  }
  double result = temp.tv_sec;
  result += temp.tv_usec / 1000000.0;

  return result; 
}

#define WALL_TIME_START(t) \
  timeval timer_##t##_current_, \
  timer_##t##_last_ = get_wall_time()

#define WALL_TIME_PRINT(format,t)		 \
  timer_##t##_current_ = get_wall_time(); \
  printf(format, subtract_time(timer_##t##_current_, timer_##t##_last_)); \
  timer_##t##_last_ = timer_##t##_current_

#else

inline
timespec
get_wall_time()
{
  struct timespec tp;
  clock_gettime(CLOCK_REALTIME, &tp);
  return tp;
}

inline
double
subtract_time(timespec& tp_new, timespec& tp_old)
{
  timespec temp = tp_new;
  temp.tv_sec -= tp_old.tv_sec;
  temp.tv_nsec -= tp_old.tv_nsec;
  if (temp.tv_nsec < 0) {
    --temp.tv_sec;
    temp.tv_nsec += 1000000000;
  }
  double result = temp.tv_sec;
  result += temp.tv_nsec / 1000000000.0;

  return result; 
}

#define WALL_TIME_START(t) \
  timespec timer_##t##_current_, \
  timer_##t##_last_ = get_wall_time()

#define WALL_TIME_PRINT(format,t)		 \
  timer_##t##_current_ = get_wall_time(); \
  printf(format, subtract_time(timer_##t##_current_, timer_##t##_last_)); \
  timer_##t##_last_ = timer_##t##_current_
#endif

#else

#define TIME_START(t)
#define TIME_CHECK(t)
#define TIME_PRINT(format,t)
#define TIME_GET(t) 0.0
#define WALL_TIME_START(t)
#define WALL_TIME_PRINT(format,t)

#endif

#endif // _TIMING_HPP
