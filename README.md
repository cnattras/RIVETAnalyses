# RIVETAnalyses
This repository is for RIVET analyses developed at UTK and in the 2020 Rivetizing Heavy Ion Collisions at RHIC workshop.  The goal is that these analyses will eventually be included in the official version of RIVET.  After analyses have been officially committed, they may be removed from this repository.

Please see individual READMEs for status on PHENIX papers.  Please make sure to check the .info file before use as well.

Analyses which need work (see README for individual analyses): PHENIX_2024_IPPG252, PHENIX_2019_I1672476, PHENIX_2013_I1227971

Analysis ready for higher statistics sorted into collision system and energy: 
AUAU39:  PHENIX_2012_I1107625, PHENIX_2016_I1394433, PHENIX_2019_I1672476

AUAU62: PHENIX_2012_I1107625, PHENIX_2014_I1273625, PHENIX_2016_I1394433, PHENIX_2019_I1672476

AUAU130: PHENIX_2001_I562409, PHENIX_2001_I555603, PHENIX_2014_I1273625, PHENIX_2016_I1394433
 
AUAU200: PHENIX_2003_I619987, PHENIX_2004_I624474, PHENIX_2006_I711951, PHENIX_2007_I731133, PHENIX_2008_I777211, PHENIX_2008_I778168,  PHENIX_2009_I815824, PHENIX_2010_I856259, PHENIX_2011_I900703, PHENIX_2012_I1116179, PHENIX_2013_I1227971, PHENIX_2014_I1273625, PHENIX_2016_I1394433

DAU200: PHENIX_2006_I711951, PHENIX_2007_I731133, PHENIX_2013_I1227971, PHENIX_2014_I1273625, PHENIX_2016_I1394433

PP39: PHENIX_2012_I1107625

PP62: PHENIX_2011_I886590, PHENIX_2012_I1107625 

PP200: PHENIX_2006_I711951, PHENIX_2007_I731133, PHENIX_2008_I777211, PHENIX_2008_I776624, PHENIX_2009_I815824 ,PHENIX_2010_I856259, PHENIX_2011_I886590, PHENIX_2014_I1273625, PHENIX_2024_IPPG252

CUCU200: PHENIX_2008_I776624, PHENIX_2016_I1394433

AUAU7: PHENIX_2016_I1394433

AUAU14: PHENIX_2016_I1394433

AUAU19: PHENIX_2016_I1394433

AUAU27: PHENIX_2016_I1394433

CUCU62: PHENIX_2016_I1394433

CUAU200: PHENIX_2016_I1394433


Analyses implemented at UTK
STAR Analyses updates: 

0604018/STAR_2006_I715470: Amber Threlkeld
- Correlation analysis
- This one needs more work

1004.2377/STAR_2010_I851937: Eden Ross
- Correlation Analysis
- Some histograms to be implemented
- Bkg subtraction with v2 needs to be tested and debugged
- "Segmentation fault"

1107.2955/STAR_2012_I918779: Christal Martin
- Spectra and ratio analysis
- Some histograms to be implemented (Ratios)

1110.0121/ALICE_2012_I930312: Shelby
- Not implemented

1110.0579/STAR_2012_I930463 - Implemented by Christal Martin at UTK. Charged pi, K, p and K0S and rho spectra in Au+Au and p+p at 200 GeV.
Outstanding questions/issues:
Need to implement RAA Ratios (RAA/RAA)
Need to normalize RAA

-"Error in STAR_2012_I930463:beam=AUAU:cent=GEN::analyze method: Missing implementation for a Cuts::Quantity."

1110.5800/STAR_2012_I943192 - Started by Christine Nattrass.  Dihadron correlations at 62.4 and 200 GeV in Cu+Cu and Au+Au.  Not complete.  Needs a lot of work.  Probably should be treated virtually like an unstarted analysis.
-"Error in PHENIX_2012_I1107625:beam=AUAU:cent=GEN::analyze method: Missing implementation for a Cuts::Quantity."

1603.05477/STAR_2016_I1429700 - First attempt by Shelby Richmond Fall 2019.  Identified K0-h, Lambda-h, h-Lambda, and h-K0S correlations in d+Au, Cu+Cu, and Au+Au collisions at 200 GeV.  Christine's thesis analysis.  I think this one needs more serious work.  The HEPData may not be rivet-compatible.

1604.01117/STAR_2016_I1442357 - First attempt by Erica Irwin Fall 2019.  Progress Brittney Contreras Summer 2020.  Gamma-hadron correlations in Au+Au and p+p collisions at 200 GeV.  I believe this may mainly need a higher stats test.  Can ask Brittney for a reminder.

1812.10224/STAR_2019_I1711377 - Implemented by Christal Martin at UTK.  D0 meson spectra in Au+Au at 200 GeV.
Outstanding questions/issues:
Need to implement cross-sections
Need to normalize RAA
Does RAA and RCP have D0 defined as D0+D0bar/2? Is it scaled correctly in the finalize section (yields)?

1910.04812/CMS_2020_I064906 - Implemented by Lauren Kasper at Vanderbilt.  Lambda, K0S, Omega, and Xi production in pp and pPb collisions at 5.02.  Request input from Lauren.

