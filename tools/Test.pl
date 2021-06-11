#!/bin/perl
#@analyses = (619063,731133,711951,776624,777211,778168,793126,837075,851937,886590,900703,1107625,930463,1227971,1393529,1658594,1672859,1798493);
@analyses = (619063,711951,776624,777211,778168,793126,837075,886590,1107625,930463,1227971,1393529,1658594,1672859,1798493);
for($i=0;$i<=$#analyses;$i++){
	$test = `ls ../ | grep $analyses[$i]`;
	#print "$test";
	chomp($test);
	#system "ls $test/*.sh"
	chdir("../$test");
	#system "bash RunAnalysis.sh";
	print "Analysis $analyses[$1] $test :\n";
	system "rivet-merge --pwd -o /tmp/Rivet.yoda Rivet.yoda Rivet.yoda";
	print "\n\n\n";
}
