.plot complete 
.cc: added books and beam, needs to be implemented in the analyze and finalize

Cross section has wrong units!  Off by six orders of magnitude.  Cross section is implemented in Rivet as picobarns, and in the paper it's millibarns.  Easy fix.

Xi, R, and Zg distributions are clearly not normalized correctly.  I think this is not taking the cross section into account.  I don't know off the cuff how to fix this but it should be straightforward when it's done.
