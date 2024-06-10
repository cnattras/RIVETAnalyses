# RIVETAnalyses
This repository is for RIVET analyses developed at UTK and in the 2020 Rivetizing Heavy Ion Collisions at RHIC workshop.  The goal is that these analyses will eventually be included in the official version of RIVET.  After analyses have been officially committed, they may be removed from this repository.

Here is a list of analyses and a brief summary of their status.  Please make sure to check the .info file before use as well.

Analyses implemented at UTK
0603010/PHENIX_2006_I711951: Andi Mankolli
- Spectra and ratio analysis
- Need to implement ratio plots between particle species and R_dAu.
- Need to implment double ratios
- Remove secondary decays "by hand", need to implment primary Particles

0604018/STAR_2006_I715470: Amber Threlkeld
- Correlation analysis
- This one needs more work

PHENIX_2006_I0611006: Lauren Kasper
- Spectra and ratio Analysis
- Normalizations needs double check

0801.1665/PHENIX_2008_I777211 - Implemented by Christal Martin at UTK.  Pi0 RAA in 0-5% central Au+Au at 200 GeV.  Complete.
- "Error in PHENIX_2008_I777211:beam=AUAU200:cent=GEN::analyze method: Missing implementation for a Cuts::Quantity."

0801.4020/PHENIX_2008_I778168: Mani
- Spectra and ratio Analysis
- Check ratio normalizations
- Use bin center instead of particle pt for normalizations

- "Error in PHENIX_2008_I778168:beam=AUAU:cent=GEN::analyze method: Missing implementation for a Cuts::Quantity."

08014545/PHENIX_2008_I778396: Nora Bauer
- Correlation Analysis
- Need some double check (Many histograms and very long code)

0903.3399/PHENIX_2009_I815824: John Bridges
- Correlation Analysis
- Histograms booked but not properly implemented in analyze

0903.4886/PHENIX_2009_I816486 - Implemented by Christal Martin at UTK.  Pi0 RAA vs rxn plane Au+Au 200 GeV
Outstanding questions/issues:
Need to implement v_2 and phi vs RAA.
Need to normalize RAA
CN: Maybe we can get Takahito's help with this one, with the reaction plane?

1004.2377/STAR_2010_I851937: Eden Ross
- Correlation Analysis
- Some histograms to be implemented
- Bkg subtraction with v2 needs to be tested and debugged

- "Segmentation fault"

1006.1347/PHENIX_2010_I857187: Jacob
- Correlation Analysis
- Seems in good shape. Should run and plot histograms for debugging

1010.1521/PHENIX_2011_I872172: Andrew Bryant
- Correlation Analysis
- Some histograms to be implemented (I_AA)
- Bkg subtraction to be implemented

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


1204.1526/PHENIX_2012_I1107625 - Implemented by Christal Martin at UTK.  Pi0 RAA at 39, 62.4, 200 GeV.  Complete.

1208.2254/PHENIX_2013_I1127262 - Scarcely started.  Pi0 production vs reaction plane in Au+Au at 200 GeV.

1212.3323/PHENIX_2013_I1207323 - Implemented by Brandi Skipworth Fall 2019.  Gamma-hadron correlations in Au+Au collisions at 200 GeV.  I think it was close to done, but it needs checking and higher statistics tests.

1304.3410/PHENIX_2013_I1227971 - Implemented by Christal Martin at UTK.  Spectra and ratios of charged pi/k/p in Au+Au and d+Au collisions at sqrt(s_NN)=200 GeV
Outstanding questions/issues:
Check RCP
Check finalize section (Ratio of yields, RCP, Ratio of Spectra)

1603.05477/STAR_2016_I1429700 - First attempt by Shelby Richmond Fall 2019.  Identified K0-h, Lambda-h, h-Lambda, and h-K0S correlations in d+Au, Cu+Cu, and Au+Au collisions at 200 GeV.  Christine's thesis analysis.  I think this one needs more serious work.  The HEPData may not be rivet-compatible.

1604.01117/STAR_2016_I1442357 - First attempt by Erica Irwin Fall 2019.  Progress Brittney Contreras Summer 2020.  Gamma-hadron correlations in Au+Au and p+p collisions at 200 GeV.  I believe this may mainly need a higher stats test.  Can ask Brittney for a reminder.

1812.10224/STAR_2019_I1711377 - Implemented by Christal Martin at UTK.  D0 meson spectra in Au+Au at 200 GeV.
Outstanding questions/issues:
Need to implement cross-sections
Need to normalize RAA
Does RAA and RCP have D0 defined as D0+D0bar/2? Is it scaled correctly in the finalize section (yields)?

1910.04812/CMS_2020_I064906 - Implemented by Lauren Kasper at Vanderbilt.  Lambda, K0S, Omega, and Xi production in pp and pPb collisions at 5.02.  Request input from Lauren.

Analyses implemented as part of Rivetizing Heavy Ion Collisions at RHIC
PHENIX_2008_I776624/PHENIX_2008_I776624
PHENIX_2011_I886590/PHENIX_2011_I886590
PHENIX_2011_I900703/PHENIX_2011_I900703
PHENIX_2012_I1116179/PHENIX_2012_I1116179
PHENIX_2016_I1393529/PHENIX_2016_I1393529
PHENIX_2016_I1394433/PHENIX_2016_I1394433
PHENIX_2018_1672859/PHENIX_2018_1672859

PHENIX_2018_I1672476/PHENIX_2018_I1672476
Done by Zhangdong Sun.  HEPData needs to be uploaded!  This analysis was started but the analysis is empty.

PHENIX_2019_I1672133/PHENIX_2019_I1672133
PHENIX_2020_I1773662/PHENIX_2020_I1773662
STAR_2003_I619063/STAR_2003_I619063
STAR_2009_I793126/STAR_2009_I793126
STAR_2010_I837075/STAR_2010_I837075
STAR_2010_I840766/STAR_2010_I840766
STAR_2016_I1420183/STAR_2016_I1420183
STAR_2020_I1771348/STAR_2020_I1771348



Analyses with all problems fixed, just needs Histograms turned off: 
PHENIX_2012_I1107625
PHENIX_2008_I777211
PHENIX_2008_I778168
PHENIX_2011_I886590


Analyses ready for higher stats with collision system:
PHENIX_2001_I562409: pp130 and AuAu130
PHENIX_2006_I711951: pp200, AuAu200, and DAu200
PHENIX_2007_I731133: pp200, AuAu200, and DAu200
PHENIX_2004_I624474 : AuAu200