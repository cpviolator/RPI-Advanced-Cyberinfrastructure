# Serial 2D HMC+SwendsenWang phi^4 theory

## Features

This code contains many of the structures one might expect to find in a finite element
field theory computation. We list the features here and place more elaborate
explanations in the code itself.

1. Hybrid Monte-Carlo
2. Swendsen-Wang cluster decomposition
3. Moment of phi measurements

## Using the code.

All of the variables are hardcoded for ease of development. Physical and algorithmic
variables such as \lambda, \musqr are defined the paramater structure, and can be
changed in the *main* file. The square lattice size L is a #define at line 17

To compile, simply type make in the source directory, and run with ./phi4Serial

Happy Hacking!