# Massey-pari

This folder contains a program written by Eric Ahlqvist which computes Massey products in étale cohomology of the ring of integers of a number field plus some more which is explained in the output. The program is written in C using the library [PARI](http://pari.math.u-bordeaux.fr/) which makes it very fast for number fields of low degree. 

The program originates from the mathematics article <a href="https://www.degruyter.com/document/doi/10.1515/crelle-2025-0006/html">Massey products in the étale cohomology of number fields</a> written by Eric Ahlqvist and Magnus Carlson. 

To run the main program: ./massey p "pol(s)"
where you replace pol(s) by a defining polynomial of the number field K in the variable s, and replace p by a prime number dividing the class number. 

This repository is based on Eric Ahlqvist's Massey-pari program, developed
for the computations in Ahlqvist--Carlson, “Massey products in the étale
cohomology of number fields.” The present version adds prescribed-character
secondary norm computations, rank-three tensor and relation calculations, and
the verification routines used in the accompanying work.

See `REPRODUCING.md` for build and reproduction commands.
