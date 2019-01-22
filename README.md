# LiftSRW
An algorithm to estimate the graphlet counts by a method we call lifting. The method is described and analyzed in [this paper](pdf/graphlet_lift.pdf), with some theoretical results in the [supplementary materials](pdf/graphlet_lift_supp.pdf). We also have a [Jupyter notebook](Lift.ipynb) used for experiments.

Uses [pynauty-0.6.0](https://web.cs.dal.ca/~peter/software/pynauty/html/), a wrapper for [nauty](http://pallini.di.uniroma1.it/), by Peter Dobcsányi.

# TODO:
- [x] Add tests via small test graphs for which frequency counts are known and hardcode those.
- [ ] Remove dependence on networkx.
- [ ] Speed up isomorphism.
- [ ] Rethink symbolic computation.
