# Schrödinger’s Bridge toolbox
Genji Kawakita and Masafumi Oizumi (The University of Tokyo) \
Email: gkawakita (at) g.ecc.u-tokyo.ac.jp, c-oizumi (at) g.ecc.u-tokyo.ac.jp

This repository contains a MATLAB code for computing the optimal control cost in stochastic systems by solving the Schrödinger’s bridge problem (solveSBP.m).
Please run "SB_demo.m" first for understanding how to use the code.
For details, see: Quantifying brain state transition cost via Schrödinger’s bridge \
Genji Kawakita, Shunsuke Kamiya, Shuntaro Sasai, Jun Kitazono, Masafumi Oizumi \
bioRxiv: https://doi.org/10.1101/2021.05.24.445394

For solving the Schrödinger’s bridge problem (mathematically equivalent to the entropy-reguralized optimal transport problem), Sinkhorn's algorithm implemented by Marco Cuturi (sinkhornTransport.m) is used. \
https://marcocuturi.net/SI.html
