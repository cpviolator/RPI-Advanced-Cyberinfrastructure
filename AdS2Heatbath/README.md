# Phi4 theory on AdS2 with Heatbath

## Getting started

The repository should compile out of the box. Simply type

    make

in the source directory. There are several variables you can pass at runtime,
as seen in `run.sh`.

## Details

The main purpose of this code is to implement a heatbath algorithm on a q triangulated
AdS2 lattice. The usual phi4 parameters are passed (mass^2, lambda) to the
simulation, as well as a parameter `deltaPhi`. This parameter governs how much a
candidate phi value varies from the original one. You will find that a larger value of
deltaPhi will give a lower acceptance rate, and a smaller value gives a larger one.
However, high acceptance rates usually mean the algorithm is moving through phase space
very slowly. The best choice for heathbath is anything from 50%-70% acceptance.

The author advises that one errs on the side of caution. It is better to use a 50%
acceptance rate and collect more data, than to have a higher acceptance rate with
possibly autocorrelated data.

Happy Hacking!