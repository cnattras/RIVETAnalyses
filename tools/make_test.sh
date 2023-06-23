#!/bin/bash

#Thanks to ChatGPT for helping with writing!

# Iterate through each argument given in terminal
for dir in "$@"; do

    if [ -d ../"$dir" ]; then

        #Here, we copy the RunAnalysis.sh in variable directory ($dir) and rename it test.sh and place it into the same directory
        cp "../$dir/RunAnalysis.sh" "../$dir/test.sh"
        echo "rivet-merge --pwd -o /tmp/Rivet.yoda Rivet.yoda Rivet.yoda" >> ../$dir/test.sh
        echo "$dir"
        
        
    else
        echo "Directory $dir does not exist"
        
    fi
    
done
