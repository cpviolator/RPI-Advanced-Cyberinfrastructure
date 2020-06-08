# 2D-Schwinger Model

C++ code for a dynamic or quenched 2D Schwinger Model with Wilson fermions

## Features

The code is 'built for comfort, not for speed'. No attempts have been made to
optimise the code for performance. 

### Actions

We use the Wilson gauge and fermion actions throughout.

### Simulation

We use Leapfrog HMC to simulate dynamic fermions. You have the option to perform dynamic
or quenched simulations, dictated from the command line.

### Utilities

There are several utilities at your disposal. Please refer to the comments in
the source for further information. We list here the features as a synopsis:

   1. Gauge field saving/loading.
   2. Gauge field smearing (APE)
   3. Linear algebra suite (linAlgHelpers.h)

### Measurements

At the preset time, one may measure

   1. Average plaquette
   2. Nx(N-1), (N-1)xN, and NxN Wilson loops (for Creutz ratios)
   3. Polyakov loops
   4. Topological charge
   5. Vacuum trace
   6. Pion correlation function

### Usage

In the `wilson` directory we provide two files,

   Makefile
   main.cpp

and in the `scripts` directory there is a launching script,

   launcher.sh

We have deliberatly left some parameters of the code as preprocessor defines.
This is so that any attempt to parallelise the code with, say OpenACC, will
be easier. In order to constrcut an exectuable, edit the `main.cpp` file so
that it has the size that you desire, and type `make`. Then, edit the 
`launcher.sh` script to suit your run-time paramters.

Notice that the launcher script creates directories. The executable will expect
these directories to exist, and will hard bork if they are not there.

Make sure you make the `launcher.sh` script executable:

`chmod +x launcher.sh`

Happy Hacking!
