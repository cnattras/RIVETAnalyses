Written by Austin Schmier.  Contributions by Adam Tilley, Antonio Da Silva, and Sean Grace.

This is a spectra analysis from 62.4 GeV pp.  It should be pretty easy to wrap up.

It's unclear if the spectra are normalized correctly because our test file is completely wrong.  We need to test with higher statistics!  This analysis should be prioritized because it would be a really interesting analysis to have in Rivet.

The counters ("sow_*") should be turned off in the output.
Please add protection to prevent dividing by zero in finalize.
This should be written so it can either take a beam option set on the command line or read from the HEPMC.

This needs the bin shift correction applied.

Update: Added bin shift correction and divide by zero protection. Added beam options.
    - histograms turned off

This has some bug with the spectra.  The normalization is off by orders of magnitudes.

Update: I believe this is off for the same reason as some of the other analyses and is not being divided by 2pi in the normalization. I have added it and we will see if it fixes it in our next higher stats test.
