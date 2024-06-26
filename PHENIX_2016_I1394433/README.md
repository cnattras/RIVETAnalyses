This is a transverse energy analysis with several different collision systems and energies.  Originally written by Ejiro.  Contributions from Sean Grace, Nik Nelson, Joesph Beller.

Most histograms are not implemented.

The .info file points to the wrong article on the arxiv.

The local centrality definitions/files should be removed to point to something which is in the main repository.  I won't clean them up now to avoid breaking the analysis.

The logic in init needs to be changed.  The histograms should always be declared, even if they're not filled.

The logic in analyze should also be changed.  Each event has one system/energy so you should always add the energy.  However, you fill the histogams in different places depending on the system.

Needs beams set up to take either a command line argument or to read it from the HEPMC.

Double check output - I don't know if this is calculating the average per event in a centrality bin or if it needs to be normalized in finalize.

Update: Fixed the .info file.
Update: Fixed the calibration location.
Update: Fixed beams set up.
Update: Fixed logic in init.