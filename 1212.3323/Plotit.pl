$filename = $ARGV[0];

#will compile and print out errors
system "rivet-buildplugin Rivet$filename.so $filename.cc";

#runs the assigned generator file
system "rivet --pwd -a $filename /nics/a/proj/UTK0019/cnattras/rivet/data/pyhi.hepmc";

#makes plots with a html link 
system "rivet-mkhtml --pwd Rivet.yoda";
