# LiftSRW
An algorithm to estimate the motif (or graphlet) counts by "lifting" a graph node to a connected induced subgraph of size k.

The repository includes [the paper](graphlet_lift.pdf) describing the algorithm with theoretical analysis and experiments, includes [supplementary material](graphlet_lift_supp.pdf) for theoretical bound of the variation of the estimate, and includes the [Python code](Lift.ipynb) used for experiments.

# TODO:
- [x] Add tests via small test graphs for which frequency counts are known and hardcode those.
- [ ] Remove dependence on networkx (Why? Slow. Also the atlas only goes up to 7 node graphlets.)
- [ ] Profile the code. https://docs.python.org/2/library/profile.html
- [ ] Speed up by reducing the isomorphism computations (match at the end instead of constantly; memory trade for time).
