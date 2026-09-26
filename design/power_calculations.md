# Power calculations

Section 6.1 uses two-sided normal tests at nominal 5%, before Holm adjustment.
The approximate 80%-power minimum detectable absolute difference is
(z_0.975 + z_0.80) times the realised holdout HAC standard error.

For a main-OOS difference mu and holdout standard error s, the rejection
probability is Phi(−z_0.975 − mu/s) + Phi(mu/s − z_0.975). The output includes
both full-precision standard errors and the three-decimal standard errors
used in the prose; this accounts for the reported DAX probability of about
0.33. The required sample size assumes the effect remains mu, long-run
variance is constant, and standard errors scale with the inverse square root
of sample size. It is n times ((z_0.975 + z_0.80) s / abs(mu)) squared.
The seven-year description uses roughly 250 trading days per year.

These are conditional approximations, not assurances about future market
conditions or estimates of the probability that the observed conclusion is
correct. `code/appendix_outputs.R` computes the quantities from Table 3 data.
