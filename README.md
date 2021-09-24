# Sufficient Perturbation

A finite-element based numerical simulation program to investigate the required perturbation of a bubble enclosed in a Hele-Shaw channel flow for transient and long-term dynamics.

Adapted from the suplementary material for the article "The influence of invariant solutions on the transient behaviour of an air bubble in a Hele-Shaw channel." published in Proc. R. Soc. A, 2019 by J. S. Keeler et al.

Academic correspondence email: johnpaul.alexander@manchester.ac.uk.

## Installation

1. Clone or download oomph-lib from the Git repository: https://github.com/oomph-lib/oomph-lib
2. Download trillinos. See http://oomph-lib.maths.man.ac.uk/doc/html/index.html
3. Place trillinos tar file in `oomph-lib/external_distributions/trillinos/`
6. Run `./autogen.sh` in parent directory and follow prompts
7. Update `OOMPH-LIB_INSTALL_LOCATION` in Makefile of this repository
8. Run `make`

## Testing

Testing is easily run via,
```
make test
```

## Usage

* `./bubble_unsteady` - Starts from an arbitrary initial condition and time-steps. Choose the parameters for a given initial condition when prompted.

Note: This requires command line inputs. Add a `--help` to display the options.
