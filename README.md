Computing pairwise huddinge distances
=====================================

This is a program for computing all pairwise huddinge distances for
input sequences. The output is streamed to stderr with 1st line giving
number N of sequences studied, next N lines give the sequences
themselves, and the last N(N-1)/2 lines give i and j and huddinge
distance between sequence i and sequence j.


Installing and pre-requisities
==============================

As the only dependency, moder requires Boost library, which is normally installed on Linux machines. Relatively recent version
should be used. At least version 1.49 is known to work.
To get full advantage of parallellism the compiler should support openmp 4.0 which, for instance, gcc 4.9.0 and later support.
Running make in the directory of the distribution should compile 'moder'.
If wanted, you can install by running the command

	sudo make install

or if you want to install to a non-standard location use, for example

	make prefix=$HOME/usr install

which installs the binary to $HOME/usr/bin.
Running command 'moder' should give brief instructions on the command line parameters.
The distribution also includes, for internal use, an implementation of suffix array by Juha Kärkkäinen in directory CPM03.

Running
=======

The input file to moder should consists of sequences separated by new lines. If you give a fasta file as
parameter, it ignores the header lines (ones that begin with '>'). If a fasta sequence is split on
several lines, then moder considers these as separate sequences. So, using a fasta file at the moment
is not recommended. Also, currently sequences containing non-base characters, such as 'N', are ignored.

