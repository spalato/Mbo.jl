# Cross peak dynamics for the EMA model in QDs

The calculation shows the cross-peak dynamics for a more realistic example. The
system is a simple model for a QD described at the EMA level, as described in
the paper<sup>[1](#f1)</sup>. The system contains a total of 6 states: the ground state *G*, the exciton states $X_A$ and *X_B* and the biexciton states *XX_AA*, *XX_BB* and *XX_AB*. The 2D spectrum thus has 2 diagonal peaks and two cross peaks. There is an ESA contribution to each of these features, with a constant binding energy. The full model, including lineshapes, is described in the paper and its supplement.

Two configuration files are available: `ema.yaml` and `ema_BX.yaml`. The first excludes the biexciton states from the calculation by setting the transition dipole moment to 0 for all X->XX transitions.

<a id="f1">[1]</a> "S. Palato, H. Seiler, P. Nijjar, O. Prezhdo, P. Kambhampat, submitted (2019)"