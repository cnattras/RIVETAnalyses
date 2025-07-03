#!/bin/bash
#This script is not a complete package for conversion but will make the most straightforward changes
sed -i 's/applyProjection/apply/g' *.cc
 sed -i 's/DECLARE\_RIVET\_PLUGIN/RIVET\_DECLARE\_PLUGIN/' *.cc
 sed -i 's/DEFAULT\_RIVET\_ANALYSIS\_CTOR/RIVET\_DEFAULT\_ANALYSIS\_CTOR/' *.cc

