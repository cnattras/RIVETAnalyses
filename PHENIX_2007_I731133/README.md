CN 12/17/25:
This analysis appears to be done except for protection for dividing by zero and larger scale testing.

Analysis contains: Au+Au, d+Au, and p+p 200 GeV

Fixed beams checking

Needs option where you can set the centrality by hand

Some histograms need protection against division by zero in plots
d07-x01-y04.py
d07-x01-y01.py
d07-x01-y03.py
d07-x01-y05.py
d08-x01-y01.py
d07-x01-y02.py
d08-x01-y02.py
d08-x01-y03.py
Right now what's happening is that all histograms are always scaled by the cross section.  When there are no entries for a given collision system, this gives you the wrong answer and these histograms are getting filled with nan.  This is happening in finalize.  I added protections for the first set of histograms, but it needs to be added for all of them.

Merging needs to be checked (For now table because this is something I need to do)
Needs checking with high statistics
