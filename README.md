# pmclib

Population Monte Carlo (PMC) library

## Information

### Authors

Karim Benabed

Olivier Cappé, Jean-François Cardoso, Gersende Fort, Martin Kilbinger, Simon Prunet, Christian P. Robert, Darren Wraith

### Version

1.1

### Installation

Download the library from the github repository:

```bash
git clone https://github.com/martinkilbinger/pmclib
```
A directory `pmclib` is created automatically. Change into that directory, and install the code using the `python `script `waf`.
First, configure the software:

```bash
./waf configure
```

You might need to set paths to required libraries and other options. `waf` can install libraries for you. Note that some libraries are not essential for the use of `pmclib` (but optional e.g. `lapack`, `lua`).

Next, compile the code:

```bash
./waf build
```

On success, install the library:

```bash
./waf install
```

## References

`pmclib` has been used for cosmological inference in [Wraith, Kilbinger, Benabed et al. (2009)](https://arxiv.org/abs/0903.0837).
Please cite at least this paper when using `pmclib`. Cosmological model comparison using the Bayesian evidence estimated with PMC is discussed in [Kilbinger, Wraith, Robert et al. (2010)](https://arxiv.org/abs/0912.1614). 


