# OpenACC 2D HMC+SwendsenWang phi^4 theory

## Features

This code is a series of evolutions of `phi4Serial.cpp`:

1. v0 The same phi4Serial code, with timings.
2. v1 All of the HMC routines are parallelised with OpenACC
3. v2 Most of the Swendsen-Wang cluster routines are parallelised
      with OpenACC
v3. Up to you!

## Using the code.

The make file has four different compiler variations:

1. g++
2. pgi++
3. pgi++ with multicore
4. pgi++ with tesla

You must edit the makefile (uncomment the desired line) to invoke the compiler

There are four targets; `v0`, `v1`, `v2`, `v3`. Entering `make all` will 
make all four. Entering `make v0` will make just `v0`, and so on. Each 
compiler and each version will give a different exectuable.

Happy Hacking!
