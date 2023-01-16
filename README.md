# Make fishery-related parameter files for Atlantis GOA

This code is adapted from R code written by Javier Porobic and Hem Nalini Morzaria-Luna.

It writes the harvest.prm file for Atlantis GOA. 

Currently, the code builds a file that results in a background F of 25% FMSY for species where we have FMSY, and of 25% M for all other exploited species. The purpose of this background F is to account for some fishing pressure while calibrating the model, even though calibration is done to reach unfished state. The reason for this is that if we are not actually calibrating to B0 but to B1990, which was not unfished. So, applying catch removals to a 1990 system calibrated to "equilibrium"" is likely to wipe out several species if we do not build in some tolerance to fishing. The code is equipped to deal with the real fleets once we introduce those.

The code needs the following components:

1. A template with names and values of the fishery parameters.
2. Catch data, turned to F in `make_dummy_catch.R` and then converted to daily harvest rate mFC.
3. The names of the fishing fleets.
4. Some ancillary files containing some parameter vectors created in `make_ancillary_prm.R`.
5. The BGM file and Groups.csv
6. Values of FMSY from GOA stock assessments and M from parametrization.

The script `tune_mFC.R` can be used to calibrate the values of mFC by running the model repeatedly for a year changing values of mFC until the realized catch approximately matches the catch that is in the data.