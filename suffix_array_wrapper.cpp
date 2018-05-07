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
#include "suffix_array_wrapper.hpp"

const suffix_array::char_classes_t suffix_array::char_classes_dna[16] =
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

const suffix_array::char_classes_t suffix_array::char_classes_rna[16] =
    {{'A', "A"},
     {'C', "C"},
     {'G', "G"},
     {'U', "U"},
     {'W', "AU"},
     {'S', "CG"},
     {'M', "AC"},
     {'K', "GU"},
     {'R', "AG"},
     {'Y', "CU"},
     {'B', "CGU"},
     {'D', "AGU"},
     {'H', "ACU"},
     {'V', "ACG"},
     {'N', "ACGU"},
     {'n', "ACGU"}
    };
