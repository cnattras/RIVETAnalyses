This is a measurement of neutral pions in Au+Au collisions at sqrt(sNN)=200 GeV.

It should be nearly done but it needs to be updated for YODA2, which probably means going through and updating the indexing of the histograms.  Bin shift correction needs checks for division by zero.  Double check that all histograms are labeled.
5/27: R_AA vs Centrality is not binned right. Properly scaling said  R_AA calculation will need to be done over a loop of N^pp_(cent)/N^AA_(cent) and 1/N_(binary). Christine Plz help

This appears to have been initially written by Christal Martin, with updates by Adam Tilley, Olivia Bartoshesky, Aidan Hill.