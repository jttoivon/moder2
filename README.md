Installing and pre-requisities
==============================

MODER2 is implemented in C++ and has been tested in Linux platform.
As the only external dependency, the MODER2 requires the Boost library, which is normally installed on Linux machines. Relatively recent version
should be used. At least version 1.49 is known to work.
To get full advantage of parallellism the compiler should support openmp 4.0 which, for instance, gcc 4.9.0 and later supports.
Running `make` in the directory of the distribution should compile MODER2.
If on macOS/OSX you have problems compiling MODER2, it may be due to lack of openmp support. In this case,
first run `make clean`, followed by `make NOOPENMP=1`. This option unfortunately prevents running
MODER2 with multiple threads simultaneously.

You can also install MODER2 by running the command

	sudo make install

If you want to install to a non-standard location, use for example

	make prefix=$HOME/usr install

which installs the binary to $HOME/usr/bin.
Running command `moder2` should give brief instructions on the command line parameters.
The packaged also includes, for internal use, an implementation of suffix array by Juha Kärkkäinen in directory CPM03.
For visualization of PPM and ADM models, the package
includes a modified version of a program called spacek40 by Jussi Taipale.

Running
=======

The generic form of running MODER2 is as follows

    moder2 [ options ] inputfile monomer0,monomer1,...

If the inputfile has extension `fasta` or `fa`, then it is read as a fasta file.
Otherwise, the input file to MODER2 should consist of sequences separated by 'new line' characters.
That is, each sequence should appear on its own line.
Currently, sequences containing non-base characters, such as `N`, are ignored.

By default, MODER2 learns PPM models (order-zero inhomogeneous Markov chains).
To learn ADM models (order-one inhomogeneous Markov chains) use option `--model adm`.

The second parameter is a comma separated list of initial values for monomer models. These
can be given either as IUPAC sequences or as matrices. In the latter case, the option
`--matrices` should be given.

By default, dirichlet pseudo counts are used, before the count matrices are normalized. Option
`--prior dirichlet` uses pseudo count 0.01 times the initial background frequency of the corresponding nucleotide.
Option `--prior addone` uses pseudo count 0.000001.
The option `--prior none` disables the use of pseudo counts.

MODER2 learns PPM models by default. If you want to learn ADM models, give parameter `--model adm-fixed` to `moder2`.

Use the option `--cob` to specify which cob tables you want to compute. For instance, if you have given
seeds for five monomeric binding motifs, then these are referred to with indices 0, 1, 2, 3, and 4.
If you are interested only in the interaction of motifs 0 and 2, and homodimeric cases of motif 4, then
you can give the following option:
    --cob 0-2,4-4
The option `--cob all` can be given to generate all possible cob combinations. For example,
in the case of two monomeric binding motifs, this creates three cob tables: 0-0,1-1,0-1.

The dimeric cases with non-negative distance between the half-sites, contain position
that do not belong to either monomer model
By default, MODER2 learns also this gap area from the data. It
can however be restricted with option `--max-gap-learned` that gaps are learned only
for gaps of length at most the given limit. With option `--max-gap-learned -1`
the gaps are not learned. When gap positions are not learned from the data, they are filled
with uniform distribution.

If option `--unbound filename` is given, all the sequences X_i for which the posterior probability of the background model \theta_0 is higher
than the posterior probability of any other model \theta_k are dumped to this file. In essence, sequences for which the background
model is the most likely explanation are stored to this file. The purpose of this option is to be able to verify
later whether the supposed background sequences contain any remaining signal or not.

If option `--output` is given, the program writes monomer motifs, cob tables, deviation matrices and monomer weights to separate files.
By default these files are written to current directory. This behaviour can be changed with `--outputdir directory` option.
If option `--names name1,name2,...` is given, these names are used to construct the above filenames instead of TF0, TF1,...

If option `--flanks` is given, the program also computes a model including the flanking positions of the actual motif(s).
This option can be used to check that motif is long enough not to leave out any informational positions.
However, the matrices which include the flanks are not dumped to files, even if the option --output is given, but the
user has to parse these from the program output (or use the `to_html.py` utility described below).

The option `--number-of-threads n` instructs MODER2 to use 'n' parallel OpenMP threads to speed up execution time.

**Example of running MODER2:**

	./moder2 --model adm-fixed --names TFAP2A --outputdir TFAP2A_models --cob 0-0 data/TFAP2A-head-1000.seq GGGCA > result.txt

The data file data/TFAP2A-head-1000.seq included in this package contains the first 1000 reads from ENA experiment
ERX168813 (http://www.ebi.ac.uk/ena/data/view/ERX168813).

After the program has run, the full, unparsed result, will be in file 'result.txt'.
In the directory `TFAP2A_models` the following files are stored:

* \*.pfm	 	 The PPM model of motifs
* \*.adm                 The ADM model of motifs
* \*.cob	     	 The COB table of a transcription factor pair, or of a pair of binding profiles
* \*.dev	     	 Correction table kappa for each detected overlapping/gapped dimeric case
* monomer_weights.txt	 Weight for each monomer model, separated by commas

Run `moder2` without parameters to get description of all possible command line parameters.

Visualizing pfms and cob tables
===============================

Use the program `myspacek40` to visualize a pfm file to an image.

	./myspacek40 -paths -noname --logo monomer.0.pfm monomer.0.svg
	or
        ./myspacek40 -paths -noname --logo monomer.0.adm monomer.0.svg
	
	
The `myspacek40` program only supports svg output, but this format can
easily be converted to other formats using external tools.

If the python packages `numpy` and `matplotlib` are installed, the cob tables can be visualized as follows:

	./heatmap.py TFAP2A_models/TFAP2A-TFAP2A.cob

shows the visualization on display, and

	./heatmap.py TFAP2A_models/TFAP2A-TFAP2A.cob TFAP2A-TFAP2A.png
	
creates file `TFAP2A_models/TFAP2A-TFAP2A.png`. The extension of the output file
determines the format.
Run `./heatmap.py` without parameters to get more instructions.

Creating an html report of the results
======================================

If output of the above example of running MODER2 is stored in file `result.txt`, then
an html report can be generated using the command

        ./to_html.py TFAP2A result.txt

This creates a directory named `result.report` that contains all the model component (\*.pfm, \*.adm \*.cob, \*.dev)
in numerical form and also visualized form. The directory also contains the html report
`index.html` that can be viewed using a web browser.

The first parameter to `to_html.py` should be a comma separated list of factor names,
with as many factor names as there were seeds given to MODER2. An exception to this rule
is when all the seeds given to MODER2 are just different profiles of the same transcription factor.
Then, if the first parameter is, for example, 'TF', the report generator will automatically create profile names
TFa,TFb,TFc, etc, for each seed given to MODER2.

All the images in the html page are clickable and reveal more information.

Computing pairwise huddinge distances
=====================================

This part is not closely related to MODER2, but can be useful to
cluster or otherwise visualize the set of k-mers of a data set.

The program `all_pairs_huddinge` computes all pairwise huddinge distances for
input sequences. The output is by default streamed to file `huddinge.dists` with 1st line giving
number N of sequences studied and the next N lines give the sequences
themselves. Finally there is N(N-1)/2 bytes of huddinge distances
between pairs of sequences. The distance between sequences i > j is
provided by kth (unsigned) byte where k=i*(i-1)/2 + j.
