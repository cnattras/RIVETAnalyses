PHENIX_2004_I623413

Primarily done by Ralph Davenport.  Contributions by Adam Tilley, Nik Nelson, Sean Grace, and Joesph Beller.

This one is almost done but needs plot clean up.  Needs beams added as options and needs bin shift correction applied.

Line 28 uses the final state particles projection.  That is not the right projection!  We are copying the ALICE projection, which we need to mention in the approvals, but that gets the gist correctly while this does not.
This doesn't implement ANY beam check; it needs to have the beam check implemented.

CN 5/31/24: It is not plotting nicely - it hangs on certain plots.  I narrowed down the problem, but I didn't solve it.  I commented several lines out.  When certain histograms are declared, the histograms no longer plot.  It isn't clear to me why it's doing this.  There do not appear to be any overlapping bins.  (I know overlapping bins cause similar crashes.)  The histograms do not even have to be filled - as long as they are declared (at which point Rivet tries to plot them), they cause rivet-mkhtml to freeze while making plots.  I therefore think that this is something in the PHENIX_2004_I623413.yoda files.  I think it may be time to ask the Rivet developers for help understanding what's wrong, but I thought I'd check with Niseem and Christal first.