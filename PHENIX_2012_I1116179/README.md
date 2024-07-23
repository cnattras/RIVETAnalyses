Direct photon spectra

Written by Zhandong Sun.  Contributions from Adam Tilley, Sean Grace, Nik Nelson, Joesph Beller.

This one is in great shape!  Needs bin shift correction and to be set up so that it can read beams in either as options set on the command line or by reading them from the HEPMC.

Update: Added bin shift and beam options

To do:
This does not actually look at prompt photons.  It looks at all prompt particles.  That is helpful for debugging because it helps us see that histograms are filled, but it needs to be updated so it selects the CORRECT particles.  Please see PHENIX_2019_I1672476 for correct selection of direct photons.

Update:
Changed the selection to just be prompt photons.
