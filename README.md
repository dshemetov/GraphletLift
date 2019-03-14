# GraphletLift Comparisons
GraphletLift is an algorithm to estimate the graphlet counts by a method we call lifting. The method is described and analyzed in [this paper](pdf/graphlet_lift.pdf), with some theoretical results in the [supplementary materials](pdf/graphlet_lift_supp.pdf). We also have an introductory [Jupyter notebook](Introduction.ipynb).

## Instructions

To install, clone the repo and run 

```
git clone https://github.com/dshemetov/LiftSRW
cd LiftSRW
pipenv install
pipenv shell
cd pynauty-0.6.0
cd nauty
make clean
make 
cd ..
make clean
make pynauty
make user-ins
make tests
exit
```

To verify that you have installed everything correctly run

```
pipenv run python tests.py
```

which should conclude with an average time report.

The code has been tested on 
* OS Mojave 10.14.3, Python 3.5.2, gcc Developer tools
* Linux Ubunti 16.04, Python 3.5.2, gcc

See the "Introduction.ipynb" for a tutorial on usage.

The bulk of the code is contained in "lift.py". The unit tests are in "tests.py". Various conversion and batch experiment scripts are in "/scripts". Some graphs we tested were too large to fit in repo; you can find the missing graphs by searching on http://networkrepository.com/ (you may have to convert them from the .mtx format into a .edges or .edgelist format; converter script in "/scripts"). 

## Dependencies

The code uses [pynauty-0.6.0](https://web.cs.dal.ca/~peter/software/pynauty/html/), a wrapper for [nauty](http://pallini.di.uniroma1.it/) by Peter Dobcsányi, and [networkx](https://networkx.github.io/). 

We compared our method to two other modern methods: [PGD](https://github.com/nkahmed/pgd#input-file-format) and [ESCAPE](https://bitbucket.org/seshadhri/escape).
