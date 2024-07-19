PHENIX_2003_I619987

Done by Joesph Beller and Olivia Bartoshky

This one appears to be done but I get a lot of errors when making plots

Update: The errors occur from NaN values in the .dat files that are created when runing rivet-mkhtml. Updated .cc file to add exceptions to these errors in non histogram division. Needs checked by someone with more experience.

To Do:

Ratio plots are not filled correctly.  It looks like the ratio is calculated but it's not ending up in the right histogram.  Rcp is getting filled correctly.  This looks like logic errors in finalize - it should be possible to fix it.

Please also suppress output which does not include data points.
