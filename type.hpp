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
//typedef unsigned long long int big_int;

#include <cstdint>

// The types uint32_t and uint64_t are defined in c++11 so the below definitions are commented out:
//typedef unsigned uint32_t;
//typedef unsigned long long uint64_t;

#ifdef __int128
typedef unsigned __int128 myuint128;
typedef __int128 myint128;
#else
typedef __uint128_t myuint128;   // for old gcc compilers
typedef __int128_t myint128;   // for old gcc compilers
#endif

typedef unsigned long long int big_int;
//typedef myuint128 big_int;

