BEGIN PLOT /ALICE_2019_I1693308/d01-x01-y01
Title=Inclusive Charged Jet Cross Section for PP Collisions at $\sqrt{s}$ = 7 TeV
XLabel=$p^{ch\,jet}_{T}$(GeV/c)
YLabel=$d^2 \sigma^{ch\,jet}/dp_{T}d\eta$ (mb c/GeV)
RatioPlotYMax=2.5
RatioPlotYMin=-0.2
END PLOT
BEGIN PLOT /ALICE_2019_I1693308/d02-x01-y01
Title=Charged Jet Differential Cross Section Ratios for $p^{ch\,jet}_{T}$ 5-10 GeV/c
XLabel=$z^{ch}$
YLabel=$1/N_jets$ $dN/dz^{ch}$
RatioPlotYMin=0
END PLOT
BEGIN PLOT /ALICE_2019_I1693308/d03-x01-y01
Title=Charged Jet Differential Cross Section Ratios for $p^{ch\,jet}_{T}$ 10-15 GeV/c
XLabel=$z^{ch}$
YLabel=$1/N_jets$ $dN/dz^{ch}$
RatioPlotYMin=-0.1
END PLOT
BEGIN PLOT /ALICE_2019_I1693308/d04-x01-y01
Title=Charged Jet Differential Cross Section Ratios for $p^{ch\,jet}_{T}$ 15-20 GeV/c
XLabel=$z^{ch}$
YLabel=$1/N_jets$ $dN/dz^{ch}$
RatioPlotYMax=2.2
RatioPlotYMin=-0.1
# + any additional plot settings you might like, see make-plots documentation
END PLOT

# ... add more histograms as you need them ...
