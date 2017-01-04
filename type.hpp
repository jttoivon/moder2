//typedef unsigned long long int big_int;

typedef unsigned uint32_t;
typedef unsigned long uint64_t;

#ifdef __int128
typedef unsigned __int128 myuint128;
typedef __int128 myint128;
#else
typedef __uint128_t myuint128;   // for old gcc compilers
typedef __int128_t myint128;   // for old gcc compilers
#endif

typedef unsigned long long int big_int;
//typedef myuint128 big_int;

