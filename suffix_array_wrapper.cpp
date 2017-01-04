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

