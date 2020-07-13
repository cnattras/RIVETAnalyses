
$filename = $ARGV[0];

#path where call files are located 
$perlpath = "/nics/a/proj/UTK0019/Phys494Fa2019/JacobG/rivet";

#path to your working file
$filepath = "/nics/a/proj/UTK0019/Phys494Fa2019/JacobG/rivet/PHENIX/1006.1347";

#makes .cc .plot .info .yoda (if hepdata exists) files 
system "rivet-mkanalysis $filename";
system "rm $filename.plot";

#this will make histogram booking names
system "perl $perlpath/filemakers/WriteHistoBookingCode.pl $filename.yoda > histogram.txt"; 

#searches the yaml file to get table titles
system "perl $perlpath/filemakers/InputPlotLabels.pl $filepath/yaml/submission.yaml > dfile.txt"; 

#temp file that will be delete automatically
system "mkdir junk";

open(DATA, "<dfile.txt");
@d_lines = <DATA>;
close(DATA);
for ($i = 1; $i <= $#d_lines; $i++) {
	chomp($lines[$i]);
	($d[$i], $dfile[$i], $junk) = split("=", $d_lines[$i]);
	$dfile[$i]=~ s/\n//g;
	
	#searches all the table files for x and y axis labels (may need to be edited depending on yaml format)
	system "perl $perlpath/filemakers/MakePlotLabels.pl $filepath/yaml/$dfile[$i] > $filepath/junk/xyfile$i.txt";
	
	#updates .plot file
	system "perl $perlpath/filemakers/MakeMyAxesPretty.pl $d[$i] $dfile[$i] dfile.txt $filepath/junk/xyfile[$i].txt >> $filename.plot";
}

system "rm -r junk";

#updates .cc file (this will need to be edited depending on paper)
system "perl $perlpath/filemakers/Editcc.pl $filename > $filename.txt";
system "mv $filename.txt $filename.cc";

#copy and paste histogram books from histogram.txt to .cc file 
#make any neccessary corrections to .plot 
#edit .info file 
#then run Plotit.pl to get plots
