#!/usr/bin/perl
@lines = `cat MakePlotsPretty.csv`;
$analysisname = "PHENIX_2018_I1672859";
#print "Hello\n";
for($i=0;$i<=$#lines;$i++){
	chomp($lines[$i]);
	($Name,$Title,$Xaxis,$Yaxis,$Options) = split(",",$lines[$i]);
	#print "$lines[$i]\n";
	print "BEGIN PLOT /$analysisname/$Name\n";
print "Title=$Title\n";
print "XLabel=$Xaxis\n";
print "YLabel=$Yaxis\n";
print "Options=$Options\n";
print "END PLOT\n";

}
