BEGIN PLOT /ALICE_2013_I1210881/d01-x01-y01
Title=Jet spectra with R=0.2
XLabel=$p_T$
YLabel=$\frac{d^2\sigma}{dp_Td\eta}$
# + any additional plot settings you might like, see make-plots documentation
END PLOT

BEGIN PLOT /ALICE_2013_I1210881/d01-x01-y02
Title=Jet spectra with R=0.4
XLabel=$p_T$
YLabel=$\frac{d^2\sigma}{dp_Td\eta}$
# + any additional plot settings you might like, see make-plots documentation
END PLOT

BEGIN PLOT /ALICE_2013_I1210881/d02-x01-y01
Title=Ratios of R=0.2 to R=0.4 jets
XLabel=$p_T$
YLabel=$\sigma(R=0.2)/\sigma(R=0.4)$
LogY=<0>
# + any additional plot settings you might like, see make-plots documentation
END PLOT

# ... add more histograms as you need them ...
