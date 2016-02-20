/*
 * Copyright 2003 Juha K"arkk"ainen
 *
 * This file is part of DC-sort
 *
 * DC-sort is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version. 
 *
 * DC-sort is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *
 */


//======================================================================
// Simple macros for timing
//======================================================================

#ifndef TIMING2_HPP
#define TIMING2_HPP

#include <sys/time.h>
#include <iostream>
#include <iomanip>

inline double get_cpu_time2()
{
  static const double in_seconds = 1.0/static_cast<double>(CLOCKS_PER_SEC);
  return clock() * in_seconds;
}

#ifdef TIMING2

#define TIME_START2(t) \
  double timer_##t##_current_, \
         timer_##t##_last_ = get_cpu_time2()

#define TIME_CHECK2(t) \
  timer_##t##_current_ = get_cpu_time(); \
  std::fprintf(stdout, " %10.2f",  timer_##t##_current_ - timer_##t##_last_); \
  std::cout << std::flush; \
  timer_##t##_last_ = timer_##t##_current_

// this didn't work with g++ 2.95.4
//  std::cout << std::setw(10) << std::fixed << std::setprecision(2) 
//          << timer_##t##_current_ - timer_##t##_last_ << std::flush; 

#else

#define TIME_START2(t)
#define TIME_CHECK2(t)

#endif

#endif // TIMING_HPP
