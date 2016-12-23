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
#include "suffix_array_wrapper.hpp"

const suffix_array::char_classes_t suffix_array::char_classes2[16] =
    {{'A', "A"},
     {'C', "C"},
     {'G', "G"},
     {'T', "T"},
     {'W', "AT"},
     {'S', "CG"},
     {'M', "AC"},
     {'K', "GT"},
     {'R', "AG"},
     {'Y', "CT"},
     {'B', "CGT"},
     {'D', "AGT"},
     {'H', "ACT"},
     {'V', "ACG"},
     {'N', "ACGT"},
     {'n', "ACGT"}
    };

