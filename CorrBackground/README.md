This folder section have codes that uses combinatorial background to get background (both combinatorial and correlated background). Initial commit by Sijan Regmi on Sept 1, 2024.

The second commit on Dec 26 2024 is for "HistBothRange.C" to get new update to account for newer update for correlated background. As of today this will be still in testing phase and may not work.

Update April, 26, 2025

HistBothRangeV1.C is working as expected but does not produce best result. We are working on that to optimize this (we may modify our algorithm in code a bit).

Update April 27, 2025
HistBothRangeV2.C is working as expected but requires fine tune. This is different from V1 such that we directly project as per initial pairpT ranges. We are working on that to optimize this (we may modify our algorithm in code a bit).

The main issue seem to be related to choice of fit function based on error (and remove over fitting and underfitting) and may be the range as well.

Update May 1, 2025

We have newer and more roboust version for HistBothRangeV2.C named as "HistBothRangeV3.C" and semi-stable version with some level of human interferance to make it work but is better than any of it's predecessor.
