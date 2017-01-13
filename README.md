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

By default, no pseudo counts are used. Option '--prior addone' uses pseudo count 1. Option
'--prior dirichlet' uses pseudo count 0.01 times the initial background frequency of the corresponding nucleotide.

Use the option '--cob' to specify which cob tables you want to compute. For instance, if you have given
seeds for five monomeric binding motifs, then these are referred to with indices 0, 1, 2, 3, and 4.
If you are interested only in the interaction of motifs 0 and 2, and homodimeric cases of motif 4, then
you can give the following option:
    --cob 0-2,4-4

By default, MODER learns the gap area for dimers having positive gap length (d>=0). It
can however be restricted with option '--max-gap-learned' that gaps are learned only
for gaps of length at most the given limit. With option '--max-gap-learned -1'
the gaps are never learned.

If option '--unbound filename' is given, all the sequences X_i that had background mixing parameter \lambda_i0 larger
than for any other model's parameter \lambda_ik are dumped to this file, where \lambda_ik = \sum_j z_ikj. In essence, sequences for which the background
model is the most likely explanation are stored to this file. The purpose of this option is to be able to verify
later whether the supposed background sequences contain any interesting signal or not.

If option '--output' is given, the program writes monomer motifs, cob tables, deviation matrices and monomer weights to separate files.
By default these files are written to current directory. This behaviour can be changed with '--outputdir directory' option.
If option '--names name1,name2,...' is given, these names are used to construct the above filenames instead of TF0, TF1,...

Example of running moder:
./moder --names TFAP2A --outputdir TFAP2A_models --prior addone --cob 0-0 data/TFAP2A-head-1000.seq GGGCA > result.txt

The data file data/TFAP2A-head-1000.seq contains the first 1000 reads from ENA experiment ERX168813 (http://www.ebi.ac.uk/ena/data/view/ERX168813).

After the program has run, the full, unparsed result, is in file 'result.txt'.
In the directory 'TFAP2A_models' the following files are stored:
   *.pfm	 	 The pfm model of monomer motifs
   *.cob	     	 The cob table of a transcription factor pair, or of a pair of binding profiles
   *.deviation	     	 Correction table kappa for each detected overlapping dimeric case
   monomer_weights.txt	 Weight for each monomer model, separated by commas

If the R package 'pheatmap' is installed, the cob tables can be visualized as follows:
./heatmap.R -c TFAP2A_models/TFAP2A-TFAP2A.cob
creates file TFAP2A_models/TFAP2A-TFAP2A.png
./heatmap.R -c -s TFAP2A_models/TFAP2A-TFAP2A.cob
creates file TFAP2A_models/TFAP2A-TFAP2A.svg
Run ./heatmap.R without parameters to get brief instructions.
