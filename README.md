# LiftSRW
An algorithm to estimate the graphlet counts by a method we call lifting. The method is described and analyzed in [this paper](pdf/graphlet_lift.pdf), with some theoretical results in the [supplementary materials](pdf/graphlet_lift_supp.pdf). We also have a [Jupyter notebook](Lift.ipynb) used for experiments.

# TODO:
- [x] Add tests via small test graphs for which frequency counts are known and hardcode those.
- [ ] Remove dependence on networkx (Why? Slow. Also the atlas only goes up to 7 node graphlets.)
- [ ] Profile the code. https://docs.python.org/2/library/profile.html
- [ ] Speed up by reducing the isomorphism computations (match at the end instead of constantly; memory trade for time).
