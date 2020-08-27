# NucEnergy

Nucleosome energy predictor based on DNA sequence.

This software is a port of the python script by [van der Heijden et al.](https://www.pnas.org/content/109/38/E2514/1) with minor bug fixes. This software is written in C++ to enable faster analysis of large numbers of sequences.

## Manual

To compile the NucEnergy software, navigate to the "Release" directory and use the command `make all` to build the binary.

Usage instructions: `NucEnergy -i [input file] -o [output file] -w [window size (74)] -b [amplitude (0.2)] -p [period (10)]`

Default values (in parentheses) are taken from the default optimal values in van der Heijden et al.

The input file should be a tab delimited file with the following fields: `chr start stop DNA_sequence'.

## Acknowledgements and Contact

This algorithm was originally published and written in Python by [van der Heijden et al.](https://www.pnas.org/content/109/38/E2514/1). It was ported to C++ with minor bug fixes by Rohan Shah (rohanshah@uchicago.edu).

For questions, comments, or issues with the software, please contact Rohan Shah (rohanshah@uchicago.edu) or Alex Ruthenburg (aruthenburg@uchicago.edu).
