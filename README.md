# GraphletLift Comparisons
GraphletLift is an algorithm to estimate the graphlet counts by a method we call lifting. The method is described and analyzed in [this paper](pdf/graphlet_lift.pdf), with some theoretical results in the [supplementary materials](pdf/graphlet_lift_supp.pdf). We also have a [Jupyter notebook](Lift.ipynb) used for experiments.

The code uses [pynauty-0.6.0](https://web.cs.dal.ca/~peter/software/pynauty/html/), a wrapper for [nauty](http://pallini.di.uniroma1.it/), by Peter Dobcsányi.

We compare our method to two other methods: [PGD](https://github.com/nkahmed/pgd#input-file-format) and [ESCAPE](https://bitbucket.org/seshadhri/escape).

