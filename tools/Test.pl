#!/bin/perl
#@analyses = (619063,731133,711951,776624,777211,778168,793126,837075,851937,886590,900703,1107625,930463,1227971,1393529,1658594,1672859,1798493);
@analyses = (619063,711951,776624,777211,778168, 793126,837075,886590,1107625,930463, 1227971,1393529,1658594,1672859,1798493);
@isAu200 =   (1,0,0,1,1,  0,0,0,0,1, 1,1,1,0,1);
@isAu130 =   (1,0,0,0,0,  1,0,0,0,0, 0,0,0,0,0);
@ispp200 =   (1,1,1,1,1,  0,1,1,0,1, 1,1,0,1,1);
@isCuCu200 = (0,0,1,0,0,  0,1,0,0,0, 0,0,0,0,0);
@isdAu200 =  (0,1,0,0,0,  1,0,0,0,0, 1,0,0,0,1);
@isAu62 =    (0,0,0,0,0,  1,0,0,1,0, 0,0,0,0,0);
@ispp62 =    (0,0,0,0,0,  0,0,1,1,0, 0,0,0,0,0);
@ispp39 =    (0,0,0,0,0,  0,0,0,1,0, 0,0,0,0,0);
@isAu39 =    (0,0,0,0,0,  0,0,0,1,0, 0,0,0,0,0);
@isCuAu200=  (0,0,0,0,0,  0,0,0,0,0, 0,0,0,1,0);


$workdir = "/phenix/scratch/cen/sampleHepData/LocalAnalyses/";
$AuAuCommand = "";
$ppCommand = "";
for($i=0;$i<=$#analyses;$i++){
	$test = `ls ../ | grep $analyses[$i]`;
	#print "$test";
	chomp($test);
	#system "ls $test/*.sh"
	chdir("../$test");
	#system "bash RunAnalysis.sh";
	#print "Analysis $analyses[$1] $test :\n";
	system "cp ~/RIVETAnalyses/$test/$test.* $workdir/.";
	print "rivet-build Rivet$test.so $test.cc\n";
	#system "rivet-merge --pwd -o /tmp/Rivet.yoda Rivet.yoda Rivet.yoda";i
	#print "\n\n\n";
	if($isAu200[$i]==1){$AuAu200=" -a $test:cent=GEN $AuAu200";}
        if($isAu130[$i]==1){$AuAu130=" -a $test:cent=GEN $AuAu130";}
        if($ispp200[$i]==1){$pp200=" -a $test:cent=GEN $pp200";}
        if($isCuCu200[$i]==1){$CuCu200=" -a $test:cent=GEN $CuCu200";}
        if($isdAu200[$i]==1){$dAu200=" -a $test:cent=GEN $dAu200";}
        if($isAu62[$i]==1){$AuAu62=" -a $test:cent=GEN $AuAu62";}
        if($ispp62[$i]==1){$pp62=" -a $test:cent=GEN $pp62";}
        if($ispp39[$i]==1){$pp39=" -a $test:cent=GEN $pp39";}
        if($isAu39[$i]==1){$AuAu39=" -a $test:cent=GEN $AuAu39";}
        if($isCuAu200[$i]==1){$CuAu200=" -a $test:cent=GEN $CuAu200";}
}
print "AuAu200:\n$AuAu200\n\n";
print "AuAu130:\n$AuAu130\n\n";
print "pp200:\n$pp200\n\n";
print "CuCu200:\n$CuCu200\n\n";
print "dAu200:\n$dAu200\n\n";
print "AuAu62:\n$Au62\n\n";
print "pp62:\n$pp62\n\n";
print "pp39:\n$pp39\n\n";
print "Au39:\n$Au39\n\n";
print "CuAu200:\n$CuAu200\n\n";

